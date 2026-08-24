//! 64-bit confidential ERC-20 transfer.
//! Circuit: LT + SUB + 2*SELECT + ADD.

use std::time::Instant;
use tfhe::prelude::*;
use tfhe::{generate_keys, set_server_key, ConfigBuilder, FheUint64};

fn main() {
    let (client_key, server_key) = generate_keys(ConfigBuilder::default());
    set_server_key(server_key);

    let clear_sender = 1_000_u64;
    let clear_receiver = 250_u64;
    let clear_amount = 300_u64;
    let sender = FheUint64::encrypt(clear_sender, &client_key);
    let receiver = FheUint64::encrypt(clear_receiver, &client_key);
    let amount = FheUint64::encrypt(clear_amount, &client_key);
    let zero = FheUint64::encrypt_trivial(0_u64);

    let begin = Instant::now();
    // insufficient=1 selects the original sender balance and a zero credit.
    let insufficient = sender.lt(&amount);
    let sender_candidate = &sender - &amount;
    let sender_out = insufficient.select(&sender, &sender_candidate);
    let credited_amount = insufficient.select(&zero, &amount);
    let receiver_out = &receiver + &credited_amount;
    let elapsed = begin.elapsed();

    let actual_sender: u64 = sender_out.decrypt(&client_key);
    let actual_receiver: u64 = receiver_out.decrypt(&client_key);
    let rejected = clear_sender < clear_amount;
    let expected_sender = if rejected {
        clear_sender
    } else {
        clear_sender - clear_amount
    };
    let expected_receiver = if rejected {
        clear_receiver
    } else {
        clear_receiver + clear_amount
    };
    let correct = actual_sender == expected_sender && actual_receiver == expected_receiver;

    println!(
        "ERC20_64 ms={:.3} sender_out={} receiver_out={} circuit=1xLT+1xSUB+2xSELECT+1xADD correct={}",
        elapsed.as_secs_f64() * 1_000.0,
        actual_sender,
        actual_receiver,
        correct
    );
    if !correct {
        std::process::exit(2);
    }
}
