use std::error::Error;

use salus::benchmark_runner::{
    ensure_no_unknown_args, parse_common_options, run_circuit, take_option, BlifBuilder,
    CircuitSource, LutMode,
};

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1).collect::<Vec<_>>();
    let rows = take_option(&mut args, "--rows")?
        .map(|value| value.parse::<usize>())
        .transpose()?
        .unwrap_or(1);
    if rows == 0 {
        return Err("--rows must be at least 1".into());
    }
    let options = parse_common_options(
        &mut args,
        "benchmark_artifacts/tpch_q6_16",
        LutMode::Hybrid,
    )?;
    ensure_no_unknown_args(
        &args,
        "bench_tpch_q6 [--rows N] [--mode hybrid] [--repeat N] [--seed N] [--out-dir DIR] [--emit-only|--prepare-only]",
    )?;
    let source = CircuitSource {
        name: format!("tpch_q6_16_rows_{rows}"),
        blif: tpch_q6_blif(rows),
    };
    if let Some(timing) = run_circuit(&source, &options)? {
        println!(
            "tpch_q6_16_rows_{rows}_average_evaluation_only_ms: {:.2}",
            timing.evaluation_only_ms
        );
        println!(
            "tpch_q6_16_rows_{rows}_average_pure_server_compute_ms: {:.2}",
            timing.pure_server_compute_ms
        );
    }
    Ok(())
}

/// Generates a 16-bit modular/fixed-point form of TPC-H Q6. Dates are day
/// numbers, one year is 365, and one discount unit represents 0.01.
fn tpch_q6_blif(rows: usize) -> String {
    let mut circuit = BlifBuilder::new(format!("tpch_q6_16_rows_{rows}"));
    let query_date = circuit.input_word("query_date", 16);
    let query_discount = circuit.input_word("query_discount", 16);
    let query_quantity = circuit.input_word("query_quantity", 16);
    let year = circuit.constant_word(365, 16);
    let one = circuit.constant_word(1, 16);
    let zero = circuit.constant_word(0, 16);
    let end_date = circuit.add(&query_date, &year);
    let discount_low = circuit.sub(&query_discount, &one);
    let discount_high = circuit.add(&query_discount, &one);
    let mut contributions = Vec::with_capacity(rows);

    for row in 0..rows {
        let shipdate = circuit.input_word(&format!("l_shipdate_{row}"), 16);
        let discount = circuit.input_word(&format!("l_discount_{row}"), 16);
        let quantity = circuit.input_word(&format!("l_quantity_{row}"), 16);
        let extended_price = circuit.input_word(&format!("l_extendedprice_{row}"), 16);
        let date_lt_start = circuit.lt(&shipdate, &query_date);
        let date_ge_start = circuit.not(&date_lt_start);
        let date_lt_end = circuit.lt(&shipdate, &end_date);
        let discount_lt_low = circuit.lt(&discount, &discount_low);
        let discount_ge_low = circuit.not(&discount_lt_low);
        let discount_gt_high = circuit.lt(&discount_high, &discount);
        let discount_le_high = circuit.not(&discount_gt_high);
        let quantity_lt = circuit.lt(&quantity, &query_quantity);
        let matches = circuit.and_many(&[
            date_ge_start,
            date_lt_end,
            discount_ge_low,
            discount_le_high,
            quantity_lt,
        ]);
        let product = circuit.mult(&extended_price, &discount);
        contributions.push(circuit.select_word(&matches, &product, &zero));
    }

    while contributions.len() > 1 {
        let mut next = Vec::with_capacity((contributions.len() + 1) / 2);
        for pair in contributions.chunks(2) {
            next.push(if pair.len() == 2 {
                circuit.add(&pair[0], &pair[1])
            } else {
                pair[0].clone()
            });
        }
        contributions = next;
    }
    let revenue_value = contributions.pop().unwrap();
    let revenue = circuit.output_word("revenue", 16);
    circuit.connect_word(&revenue, &revenue_value);
    circuit.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generates_one_predicate_and_contribution_per_row() {
        let source = tpch_q6_blif(3);
        assert!(source.contains(".model tpch_q6_16_rows_3"));
        assert!(source.contains("l_extendedprice_2[15]"));
        assert!(source.contains(".outputs revenue[0]"));
    }
}
