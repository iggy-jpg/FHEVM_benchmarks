use std::collections::HashSet;
use std::error::Error;
use std::fmt::Write;
use std::fs;

use salus::benchmark_runner::{
    ensure_no_unknown_args, parse_common_options, run_circuit, take_option, BlifBuilder,
    CircuitSource, LutMode,
};

const ALL_WIDTHS: &[usize] = &[4, 8, 16, 32, 64, 128, 256];
const ALL_OPS: &[&str] = &["add", "sub", "mult", "gt", "lt", "eq", "ge", "le"];

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1).collect::<Vec<_>>();
    let widths = take_option(&mut args, "--widths")?
        .map(|raw| parse_widths(&raw))
        .transpose()?
        .unwrap_or_else(|| ALL_WIDTHS.to_vec());
    let operations = take_option(&mut args, "--ops")?
        .map(|raw| parse_operations(&raw))
        .transpose()?
        .unwrap_or_else(|| ALL_OPS.iter().map(|op| (*op).to_owned()).collect());
    let options = parse_common_options(
        &mut args,
        "benchmark_artifacts/integer_ops",
        LutMode::Hybrid,
    )?;
    ensure_no_unknown_args(
        &args,
        "bench_integer_ops [--widths 4,8,...] [--ops add,sub,...] [--mode hybrid] [--repeat N] [--seed N] [--out-dir DIR] [--emit-only|--prepare-only]",
    )?;

    let mut results = Vec::new();
    for width in widths {
        for operation in &operations {
            let name = format!("{operation}_{width}");
            println!("\n=== {name} ===");
            let source = CircuitSource {
                name: name.clone(),
                blif: integer_operation_blif(operation, width),
            };
            if let Some(timing) = run_circuit(&source, &options)? {
                results.push((operation.clone(), width, timing));
            }
        }
    }

    if !results.is_empty() {
        let mut csv = concat!(
            "operation,width,average_input_selector_bootstrap_ms,",
            "average_bdd_eval_ms,average_inter_lut_cb_ms,",
            "average_evaluation_only_ms,average_pure_server_compute_ms\n"
        )
        .to_owned();
        for (operation, width, timing) in &results {
            writeln!(
                csv,
                "{operation},{width},{:.2},{:.2},{:.2},{:.2},{:.2}",
                timing.input_selector_bootstrap_ms,
                timing.bdd_eval_ms,
                timing.inter_lut_cb_ms,
                timing.evaluation_only_ms,
                timing.pure_server_compute_ms,
            )?;
        }
        println!("\n{csv}");
        fs::create_dir_all(&options.output_root)?;
        let results_path = options.output_root.join("results.csv");
        fs::write(&results_path, csv)?;
        println!("results: {}", results_path.display());
    }
    Ok(())
}

fn parse_widths(raw: &str) -> Result<Vec<usize>, Box<dyn Error>> {
    let allowed = ALL_WIDTHS.iter().copied().collect::<HashSet<_>>();
    let values = raw
        .split(',')
        .map(str::trim)
        .map(str::parse::<usize>)
        .collect::<Result<Vec<_>, _>>()?;
    if values.is_empty() || values.iter().any(|width| !allowed.contains(width)) {
        return Err("--widths accepts a comma-separated subset of 4,8,16,32,64,128,256".into());
    }
    Ok(values)
}

fn parse_operations(raw: &str) -> Result<Vec<String>, Box<dyn Error>> {
    let allowed = ALL_OPS.iter().copied().collect::<HashSet<_>>();
    let values = raw
        .split(',')
        .map(str::trim)
        .map(str::to_owned)
        .collect::<Vec<_>>();
    if values.is_empty() || values.iter().any(|op| !allowed.contains(op.as_str())) {
        return Err("--ops accepts a comma-separated subset of add,sub,mult,gt,lt,eq,ge,le".into());
    }
    Ok(values)
}

fn integer_operation_blif(operation: &str, width: usize) -> String {
    let mut circuit = BlifBuilder::new(format!("{operation}_{width}"));
    let a = circuit.input_word("a", width);
    let b = circuit.input_word("b", width);
    if matches!(operation, "add" | "sub" | "mult") {
        let outputs = circuit.output_word("out", width);
        let values = match operation {
            "add" => circuit.add(&a, &b),
            "sub" => circuit.sub(&a, &b),
            "mult" => circuit.mult(&a, &b),
            _ => unreachable!(),
        };
        circuit.connect_word(&outputs, &values);
    } else {
        let output = circuit.output_bit("out");
        let value = match operation {
            "gt" => circuit.lt(&b, &a),
            "lt" => circuit.lt(&a, &b),
            "eq" => circuit.eq(&a, &b),
            "ge" => {
                let lt = circuit.lt(&a, &b);
                circuit.not(&lt)
            }
            "le" => {
                let gt = circuit.lt(&b, &a);
                circuit.not(&gt)
            }
            _ => unreachable!("validated operation"),
        };
        circuit.connect(&output, &value);
    }
    circuit.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn comparison_has_scalar_output() {
        assert!(integer_operation_blif("ge", 16).contains(".outputs out"));
    }

    #[test]
    fn multiplication_is_modular_width() {
        assert!(integer_operation_blif("mult", 4).contains("out[3]"));
    }
}
