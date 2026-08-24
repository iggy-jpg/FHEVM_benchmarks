use std::error::Error;

use salus::benchmark_runner::{
    ensure_no_unknown_args, parse_common_options, run_circuit, BlifBuilder, CircuitSource, LutMode,
};

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1).collect::<Vec<_>>();
    let options = parse_common_options(
        &mut args,
        "benchmark_artifacts/erc20_64",
        LutMode::Hybrid,
    )?;
    ensure_no_unknown_args(
        &args,
        "bench_erc20 [--mode hybrid] [--repeat N] [--seed N] [--out-dir DIR] [--emit-only|--prepare-only]",
    )?;
    let source = CircuitSource {
        name: "erc20_transaction_64".to_owned(),
        blif: erc20_blif(),
    };
    if let Some(timing) = run_circuit(&source, &options)? {
        println!(
            "erc20_64_average_evaluation_only_ms: {:.2}",
            timing.evaluation_only_ms
        );
        println!(
            "erc20_64_average_pure_server_compute_ms: {:.2}",
            timing.pure_server_compute_ms
        );
    }
    Ok(())
}

/// One private ERC-20 transfer:
/// lt + sub + 2*select + add.  `success` is one exactly when the sender has
/// enough balance. Failed transfers leave both balances unchanged.
fn erc20_blif() -> String {
    let mut circuit = BlifBuilder::new("erc20_transaction_64");
    let sender = circuit.input_word("sender_balance", 64);
    let receiver = circuit.input_word("receiver_balance", 64);
    let amount = circuit.input_word("amount", 64);
    let sender_out = circuit.output_word("sender_out", 64);
    let receiver_out = circuit.output_word("receiver_out", 64);
    let success_out = circuit.output_bit("success");

    let insufficient = circuit.lt(&sender, &amount); // FHE.lt
    let success = circuit.not(&insufficient);
    let debited = circuit.sub(&sender, &amount); // FHE.sub
    let selected_sender = circuit.select_word(&success, &debited, &sender); // select 1
    let zero = circuit.constant_word(0, 64);
    let accepted_amount = circuit.select_word(&success, &amount, &zero); // select 2
    let selected_receiver = circuit.add(&receiver, &accepted_amount); // FHE.add

    circuit.connect_word(&sender_out, &selected_sender);
    circuit.connect_word(&receiver_out, &selected_receiver);
    circuit.connect(&success_out, &success);
    circuit.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contains_required_operation_shape() {
        let source = erc20_blif();
        assert!(source.contains(".model erc20_transaction_64"));
        assert!(source.contains(".outputs sender_out[0]"));
        assert!(source.contains("receiver_out[63] success"));
    }
}
