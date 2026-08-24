//! TPC-H Query 6 with encrypted 16-bit input columns and a 32-bit revenue sum.

use std::time::Instant;
use tfhe::prelude::*;
use tfhe::{generate_keys, set_server_key, ConfigBuilder, FheUint16, FheUint32};

#[derive(Clone, Copy)]
struct Row {
    price: u16,
    discount: u16,
    ship_date: u16,
    quantity: u16,
}

fn main() {
    const DATE: u16 = 10_000;
    const NEXT_YEAR: u16 = DATE + 365;
    const DISCOUNT_LOW: u16 = 500;
    const DISCOUNT_HIGH: u16 = 700;
    const QUANTITY_LIMIT: u16 = 24;
    let rows = [
        Row { price: 1_000, discount: 600, ship_date: 10_100, quantity: 10 },
        Row { price: 2_000, discount: 900, ship_date: 10_100, quantity: 30 },
    ];

    let (client_key, server_key) = generate_keys(ConfigBuilder::default());
    set_server_key(server_key);
    let date = FheUint16::encrypt(DATE, &client_key);
    let next_year = FheUint16::encrypt(NEXT_YEAR, &client_key);
    let discount_low = FheUint16::encrypt(DISCOUNT_LOW, &client_key);
    let discount_high = FheUint16::encrypt(DISCOUNT_HIGH, &client_key);
    let quantity_limit = FheUint16::encrypt(QUANTITY_LIMIT, &client_key);
    let zero = FheUint32::encrypt_trivial(0_u32);
    let mut revenue = FheUint32::encrypt_trivial(0_u32);

    let encrypted_rows: Vec<_> = rows
        .iter()
        .map(|row| {
            (
                FheUint16::encrypt(row.price, &client_key),
                FheUint16::encrypt(row.discount, &client_key),
                FheUint16::encrypt(row.ship_date, &client_key),
                FheUint16::encrypt(row.quantity, &client_key),
            )
        })
        .collect();

    let begin = Instant::now();
    for (price, discount, ship_date, quantity) in encrypted_rows {
        let predicate = ship_date.ge(&date)
            & ship_date.lt(&next_year)
            & discount.ge(&discount_low)
            & discount.le(&discount_high)
            & quantity.lt(&quantity_limit);
        let price32 = FheUint32::cast_from(price);
        let discount32 = FheUint32::cast_from(discount);
        let term = price32 * discount32;
        revenue += predicate.select(&term, &zero);
    }
    let elapsed = begin.elapsed();

    let actual: u32 = revenue.decrypt(&client_key);
    let expected = rows
        .iter()
        .filter(|row| {
            row.ship_date >= DATE
                && row.ship_date < NEXT_YEAR
                && row.discount >= DISCOUNT_LOW
                && row.discount <= DISCOUNT_HIGH
                && row.quantity < QUANTITY_LIMIT
        })
        .fold(0_u32, |sum, row| {
            sum.wrapping_add(u32::from(row.price) * u32::from(row.discount))
        });
    let correct = actual == expected;
    println!(
        "TPCH_Q6_16 rows={} revenue={} expected={} ms={:.3} correct={}",
        rows.len(),
        actual,
        expected,
        elapsed.as_secs_f64() * 1_000.0,
        correct
    );
    if !correct {
        std::process::exit(2);
    }
}
