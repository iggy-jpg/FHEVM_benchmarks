use std::error::Error;
use std::fs;
use std::path::PathBuf;
use std::process::Command;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LutMode {
    Single,
    Hybrid,
    Multi,
}

impl LutMode {
    pub fn parse(value: &str) -> Result<Self, Box<dyn Error>> {
        match value {
            "single" => Ok(Self::Single),
            "hybrid" => Ok(Self::Hybrid),
            "multi" => Ok(Self::Multi),
            _ => Err(format!("invalid LUT mode '{value}'; use single, hybrid, or multi").into()),
        }
    }

    fn execute_flag(self) -> &'static str {
        match self {
            Self::Single => "--single-lut-mode",
            Self::Hybrid => "--hybrid-lut-mode",
            Self::Multi => "--multi-lut-mode",
        }
    }

    fn needs_merge(self) -> bool {
        !matches!(self, Self::Single)
    }
}

#[derive(Clone, Debug)]
pub struct BenchmarkOptions {
    pub output_root: PathBuf,
    pub mode: LutMode,
    pub repeat: usize,
    pub seed: Option<u64>,
    pub emit_only: bool,
    pub prepare_only: bool,
}

#[derive(Clone, Debug)]
pub struct CircuitSource {
    pub name: String,
    pub blif: String,
}

/// Small gate-level BLIF builder used by the benchmarks. The generated Boolean
/// network is subsequently optimized and mapped to LUT15 by ABC before SALUS
/// preprocessing, avoiding any dependence on RTL arithmetic parsing.
pub struct BlifBuilder {
    model: String,
    inputs: Vec<String>,
    outputs: Vec<String>,
    nodes: Vec<String>,
    next_node: usize,
    zero: Option<String>,
    one: Option<String>,
}

impl BlifBuilder {
    pub fn new(model: impl Into<String>) -> Self {
        Self {
            model: model.into(),
            inputs: Vec::new(),
            outputs: Vec::new(),
            nodes: Vec::new(),
            next_node: 0,
            zero: None,
            one: None,
        }
    }

    pub fn input_word(&mut self, name: &str, width: usize) -> Vec<String> {
        let bits = bit_names(name, width);
        self.inputs.extend(bits.iter().cloned());
        bits
    }

    pub fn output_word(&mut self, name: &str, width: usize) -> Vec<String> {
        let bits = bit_names(name, width);
        self.outputs.extend(bits.iter().cloned());
        bits
    }

    pub fn output_bit(&mut self, name: &str) -> String {
        self.outputs.push(name.to_owned());
        name.to_owned()
    }

    pub fn constant(&mut self, value: bool) -> String {
        if let Some(existing) = if value { &self.one } else { &self.zero } {
            return existing.clone();
        }
        let out = self.fresh(if value { "one" } else { "zero" });
        let body = if value { "1\n" } else { "" };
        self.nodes.push(format!(".names {out}\n{body}"));
        if value {
            self.one = Some(out.clone());
        } else {
            self.zero = Some(out.clone());
        }
        out
    }

    pub fn constant_word(&mut self, value: u64, width: usize) -> Vec<String> {
        (0..width)
            .map(|bit| self.constant(((value >> bit) & 1) != 0))
            .collect()
    }

    pub fn not(&mut self, input: &str) -> String {
        self.gate("not", &[input], "0 1\n")
    }

    pub fn and(&mut self, left: &str, right: &str) -> String {
        self.gate("and", &[left, right], "11 1\n")
    }

    pub fn or(&mut self, left: &str, right: &str) -> String {
        self.gate("or", &[left, right], "1- 1\n-1 1\n")
    }

    pub fn xor(&mut self, left: &str, right: &str) -> String {
        self.gate("xor", &[left, right], "01 1\n10 1\n")
    }

    pub fn xnor(&mut self, left: &str, right: &str) -> String {
        self.gate("xnor", &[left, right], "00 1\n11 1\n")
    }

    pub fn mux(&mut self, select: &str, when_true: &str, when_false: &str) -> String {
        let selected_true = self.and(select, when_true);
        let not_select = self.not(select);
        let selected_false = self.and(&not_select, when_false);
        self.or(&selected_true, &selected_false)
    }

    pub fn and_many(&mut self, values: &[String]) -> String {
        assert!(!values.is_empty());
        let mut level = values.to_vec();
        while level.len() > 1 {
            let mut next = Vec::with_capacity((level.len() + 1) / 2);
            for pair in level.chunks(2) {
                next.push(if pair.len() == 2 {
                    self.and(&pair[0], &pair[1])
                } else {
                    pair[0].clone()
                });
            }
            level = next;
        }
        level.pop().unwrap()
    }

    pub fn add(&mut self, left: &[String], right: &[String]) -> Vec<String> {
        assert_eq!(left.len(), right.len());
        let mut carry = self.constant(false);
        let mut result = Vec::with_capacity(left.len());
        for (a, b) in left.iter().zip(right) {
            let a_xor_b = self.xor(a, b);
            result.push(self.xor(&a_xor_b, &carry));
            let both = self.and(a, b);
            let carry_and_one = self.and(&carry, &a_xor_b);
            carry = self.or(&both, &carry_and_one);
        }
        result
    }

    pub fn sub(&mut self, left: &[String], right: &[String]) -> Vec<String> {
        assert_eq!(left.len(), right.len());
        let inverted = right.iter().map(|bit| self.not(bit)).collect::<Vec<_>>();
        let mut carry = self.constant(true);
        let mut result = Vec::with_capacity(left.len());
        for (a, b) in left.iter().zip(&inverted) {
            let a_xor_b = self.xor(a, b);
            result.push(self.xor(&a_xor_b, &carry));
            let both = self.and(a, b);
            let carry_and_one = self.and(&carry, &a_xor_b);
            carry = self.or(&both, &carry_and_one);
        }
        result
    }

    pub fn eq(&mut self, left: &[String], right: &[String]) -> String {
        assert_eq!(left.len(), right.len());
        let equal_bits = left
            .iter()
            .zip(right)
            .map(|(a, b)| self.xnor(a, b))
            .collect::<Vec<_>>();
        self.and_many(&equal_bits)
    }

    pub fn lt(&mut self, left: &[String], right: &[String]) -> String {
        assert_eq!(left.len(), right.len());
        let mut less = self.constant(false);
        let mut equal_prefix = self.constant(true);
        for (a, b) in left.iter().zip(right).rev() {
            let not_a = self.not(a);
            let bit_less = self.and(&not_a, b);
            let first_difference_less = self.and(&equal_prefix, &bit_less);
            less = self.or(&less, &first_difference_less);
            let bit_equal = self.xnor(a, b);
            equal_prefix = self.and(&equal_prefix, &bit_equal);
        }
        less
    }

    pub fn select_word(
        &mut self,
        select: &str,
        when_true: &[String],
        when_false: &[String],
    ) -> Vec<String> {
        assert_eq!(when_true.len(), when_false.len());
        when_true
            .iter()
            .zip(when_false)
            .map(|(yes, no)| self.mux(select, yes, no))
            .collect()
    }

    /// Low `width` bits of an unsigned product. Partial rows are reduced with
    /// a balanced tree of ripple adders instead of a serial shift-add chain.
    pub fn mult(&mut self, left: &[String], right: &[String]) -> Vec<String> {
        assert_eq!(left.len(), right.len());
        let width = left.len();
        let zero = self.constant(false);
        let mut rows = (0..width)
            .map(|shift| {
                (0..width)
                    .map(|bit| {
                        if bit < shift {
                            zero.clone()
                        } else {
                            self.and(&left[bit - shift], &right[shift])
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        while rows.len() > 1 {
            let mut next = Vec::with_capacity((rows.len() + 1) / 2);
            for pair in rows.chunks(2) {
                next.push(if pair.len() == 2 {
                    self.add(&pair[0], &pair[1])
                } else {
                    pair[0].clone()
                });
            }
            rows = next;
        }
        rows.pop().unwrap()
    }

    pub fn connect(&mut self, output: &str, value: &str) {
        self.nodes
            .push(format!(".names {value} {output}\n1 1\n"));
    }

    pub fn connect_word(&mut self, outputs: &[String], values: &[String]) {
        assert_eq!(outputs.len(), values.len());
        for (output, value) in outputs.iter().zip(values) {
            self.connect(output, value);
        }
    }

    pub fn finish(self) -> String {
        format!(
            ".model {}\n.inputs {}\n.outputs {}\n{}.end\n",
            self.model,
            self.inputs.join(" "),
            self.outputs.join(" "),
            self.nodes.concat(),
        )
    }

    fn gate(&mut self, label: &str, inputs: &[&str], truth_rows: &str) -> String {
        let out = self.fresh(label);
        self.nodes.push(format!(
            ".names {} {out}\n{truth_rows}",
            inputs.join(" ")
        ));
        out
    }

    fn fresh(&mut self, label: &str) -> String {
        let name = format!("n{}_{}", self.next_node, label);
        self.next_node += 1;
        name
    }
}

fn bit_names(name: &str, width: usize) -> Vec<String> {
    (0..width).map(|bit| format!("{name}[{bit}]")).collect()
}

pub fn take_flag(args: &mut Vec<String>, flag: &str) -> bool {
    if let Some(index) = args.iter().position(|arg| arg == flag) {
        args.remove(index);
        true
    } else {
        false
    }
}

pub fn take_option(args: &mut Vec<String>, flag: &str) -> Result<Option<String>, Box<dyn Error>> {
    if let Some(index) = args.iter().position(|arg| arg == flag) {
        args.remove(index);
        if index >= args.len() {
            return Err(format!("{flag} requires a value").into());
        }
        Ok(Some(args.remove(index)))
    } else {
        Ok(None)
    }
}

pub fn parse_common_options(
    args: &mut Vec<String>,
    default_output_root: &str,
    default_mode: LutMode,
) -> Result<BenchmarkOptions, Box<dyn Error>> {
    let emit_only = take_flag(args, "--emit-only");
    let prepare_only = take_flag(args, "--prepare-only");
    if emit_only && prepare_only {
        return Err("--emit-only and --prepare-only cannot be used together".into());
    }
    let output_root = take_option(args, "--out-dir")?
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(default_output_root));
    let mode = take_option(args, "--mode")?
        .map(|value| LutMode::parse(&value))
        .transpose()?
        .unwrap_or(default_mode);
    let repeat = take_option(args, "--repeat")?
        .map(|value| value.parse::<usize>())
        .transpose()?
        .unwrap_or(1);
    if repeat == 0 {
        return Err("--repeat must be at least 1".into());
    }
    let seed = take_option(args, "--seed")?
        .map(|value| value.parse::<u64>())
        .transpose()?;
    Ok(BenchmarkOptions {
        output_root,
        mode,
        repeat,
        seed,
        emit_only,
        prepare_only,
    })
}

#[derive(Clone, Copy, Debug)]
pub struct BenchmarkTiming {
    pub input_selector_bootstrap_ms: f64,
    pub bdd_eval_ms: f64,
    pub inter_lut_cb_ms: f64,
    pub evaluation_only_ms: f64,
    pub pure_server_compute_ms: f64,
}

/// Write, preprocess, and execute one generated circuit through SALUS.
/// Returns `None` when execution was intentionally skipped, otherwise all
/// average server-side timing components reported by `execute_circuit`.
pub fn run_circuit(
    source: &CircuitSource,
    options: &BenchmarkOptions,
) -> Result<Option<BenchmarkTiming>, Box<dyn Error>> {
    let source_dir = options.output_root.join("blif");
    let circuit_dir = options.output_root.join("circuits").join(&source.name);
    fs::create_dir_all(&source_dir)?;
    let source_path = source_dir.join(format!("{}.blif", source.name));
    fs::write(&source_path, &source.blif)?;
    println!("source: {}", source_path.display());

    if options.emit_only {
        return Ok(None);
    }

    let generate = sibling_binary("generate_circuit")?;
    let mut generate_cmd = Command::new(&generate);
    generate_cmd.arg(&source_path).arg(&circuit_dir);
    if options.mode.needs_merge() {
        generate_cmd.arg("--with-merge");
    }
    run_checked(&mut generate_cmd, "circuit generation")?;

    if options.prepare_only {
        return Ok(None);
    }

    let execute = sibling_binary("execute_circuit")?;
    let mut execute_cmd = Command::new(&execute);
    execute_cmd
        .arg(&circuit_dir)
        .arg(options.mode.execute_flag())
        .arg("--repeat-random")
        .arg(options.repeat.to_string());
    if let Some(operation) = direct_comparison_operation(&source.name) {
        execute_cmd.arg("--direct-comparison").arg(operation);
    }
    if let Some(seed) = options.seed {
        execute_cmd.arg("--seed").arg(seed.to_string());
    }
    let output = execute_cmd.output().map_err(|error| {
        format!(
            "failed to launch {}: {error}",
            execute.as_path().display()
        )
    })?;
    print!("{}", String::from_utf8_lossy(&output.stdout));
    eprint!("{}", String::from_utf8_lossy(&output.stderr));
    if !output.status.success() {
        return Err(format!("FHE execution failed with status {}", output.status).into());
    }
    let stdout = String::from_utf8_lossy(&output.stdout);
    let timing = BenchmarkTiming {
        input_selector_bootstrap_ms: parse_reported_time(
            &stdout,
            "average_input_selector_bootstrap_ms:",
            "input_selector_bootstrap_ms:",
        )?,
        bdd_eval_ms: parse_reported_time(&stdout, "average_bdd_eval_ms:", "bdd_eval_ms:")?,
        inter_lut_cb_ms: parse_reported_time(
            &stdout,
            "average_inter_lut_cb_ms:",
            "inter_lut_cb_ms:",
        )?,
        evaluation_only_ms: parse_reported_time(
            &stdout,
            "average_evaluation_only_ms:",
            "evaluation_only_ms:",
        )?,
        pure_server_compute_ms: parse_reported_time(
            &stdout,
            "average_pure_server_compute_ms:",
            "pure_server_compute_ms:",
        )?,
    };
    Ok(Some(timing))
}

fn direct_comparison_operation(name: &str) -> Option<&str> {
    let (operation, width) = name.rsplit_once('_')?;
    if !matches!(operation, "gt" | "lt" | "eq" | "ge" | "le") {
        return None;
    }
    width.parse::<usize>().ok()?;
    Some(operation)
}

fn parse_reported_time(
    stdout: &str,
    average_prefix: &str,
    single_prefix: &str,
) -> Result<f64, Box<dyn Error>> {
    stdout
        .lines()
        .rev()
        .find_map(|line| {
            line.strip_prefix(average_prefix)
                .or_else(|| line.strip_prefix(single_prefix))
                .and_then(|value| value.trim().parse().ok())
        })
        .ok_or_else(|| format!("execute_circuit did not report {average_prefix}").into())
}

fn run_checked(command: &mut Command, label: &str) -> Result<(), Box<dyn Error>> {
    let status = command.status()?;
    if !status.success() {
        return Err(format!("{label} failed with status {status}").into());
    }
    Ok(())
}

fn sibling_binary(name: &str) -> Result<PathBuf, Box<dyn Error>> {
    let mut path = std::env::current_exe()?;
    path.set_file_name(format!("{name}{}", std::env::consts::EXE_SUFFIX));
    if path.is_file() {
        Ok(path)
    } else {
        Err(format!(
            "missing helper binary {}; first run `cargo build --release --bins`",
            path.display()
        )
        .into())
    }
}

pub fn ensure_no_unknown_args(args: &[String], usage: &str) -> Result<(), Box<dyn Error>> {
    if args.is_empty() {
        Ok(())
    } else {
        Err(format!("unknown arguments: {}\nusage: {usage}", args.join(" ")).into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blif::BlifCircuit;

    fn evaluate_binary(operation: &str, a: u64, b: u64) -> u64 {
        let width = 4;
        let mut builder = BlifBuilder::new("test");
        let left = builder.input_word("a", width);
        let right = builder.input_word("b", width);
        let outputs = builder.output_word("out", width);
        let values = match operation {
            "add" => builder.add(&left, &right),
            "sub" => builder.sub(&left, &right),
            "mult" => builder.mult(&left, &right),
            _ => unreachable!(),
        };
        builder.connect_word(&outputs, &values);
        let circuit = BlifCircuit::parse_str(&builder.finish()).unwrap();
        let inputs = (0..width)
            .map(|bit| (a >> bit) & 1)
            .chain((0..width).map(|bit| (b >> bit) & 1))
            .collect::<Vec<_>>();
        let signals = circuit.evaluate_primary_input_bits(&inputs).unwrap();
        circuit
            .output_values(&signals)
            .unwrap()
            .iter()
            .enumerate()
            .fold(0, |value, (bit, set)| value | (set << bit))
    }

    #[test]
    fn four_bit_arithmetic_is_exhaustively_correct() {
        for a in 0..16 {
            for b in 0..16 {
                assert_eq!(evaluate_binary("add", a, b), (a + b) & 15);
                assert_eq!(evaluate_binary("sub", a, b), a.wrapping_sub(b) & 15);
                assert_eq!(evaluate_binary("mult", a, b), (a * b) & 15);
            }
        }
    }

    #[test]
    fn four_bit_comparisons_are_exhaustively_correct() {
        for a in 0u64..16 {
            for b in 0u64..16 {
                let mut builder = BlifBuilder::new("compare");
                let left = builder.input_word("a", 4);
                let right = builder.input_word("b", 4);
                let lt_out = builder.output_bit("lt");
                let eq_out = builder.output_bit("eq");
                let lt = builder.lt(&left, &right);
                let eq = builder.eq(&left, &right);
                builder.connect(&lt_out, &lt);
                builder.connect(&eq_out, &eq);
                let circuit = BlifCircuit::parse_str(&builder.finish()).unwrap();
                let inputs = (0..4)
                    .map(|bit| (a >> bit) & 1)
                    .chain((0..4).map(|bit| (b >> bit) & 1))
                    .collect::<Vec<_>>();
                let signals = circuit.evaluate_primary_input_bits(&inputs).unwrap();
                let outputs = circuit.output_values(&signals).unwrap();
                assert_eq!(outputs, vec![(a < b) as u64, (a == b) as u64]);
            }
        }
    }

    #[test]
    fn recognizes_only_integer_comparison_benchmark_names() {
        assert_eq!(direct_comparison_operation("gt_256"), Some("gt"));
        assert_eq!(direct_comparison_operation("le_16"), Some("le"));
        assert_eq!(direct_comparison_operation("add_64"), None);
        assert_eq!(direct_comparison_operation("gt_custom"), None);
        assert_eq!(direct_comparison_operation("erc20_64"), None);
    }
}
