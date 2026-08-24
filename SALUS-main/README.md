# SALUS: Large-Scale Homomorphic Circuit Synthesis via Logic-Aware LUT Optimization

This code is the implementation of the paper "SALUS: Large-Scale Homomorphic Circuit Synthesis via Logic-Aware LUT Optimization".

## Requirements

To build and run this project, you will need:

- `git`
- `make`
- `gcc` or `clang`
- `rust` and `cargo`
- `python3`


## Third-Party Dependencies

This repository uses two third-party dependencies under [`third_party`](third_party):

- [`third_party/abc`](third_party/abc): <https://github.com/berkeley-abc/abc>
- [`third_party/refined-tfhe-lhe`](third_party/refined-tfhe-lhe): <https://github.com/KAIST-CryptLab/refined-tfhe-lhe>

After cloning the repository, initialize the submodules with:

```bash
git submodule update --init --recursive
```

Build the vendored `ABC` binary with:

```bash
make -C third_party/abc
```

## Building SALUS

Build the Rust binaries with:

```bash
cargo build --release
```

## Generating a Homomorphic Circuit

```bash
cargo run --release --bin generate_circuit -- <input.{v|blif}> <output_dir> [--with-merge]
```
Options:

- `--with-merge`
  - enable LUT merging for `hybrid` / `multi` execution
- `--max-abc-iters N`
  - limit the number of repeated mapping iterations

Example: generate a preprocessed circuit from [`testcircuit/c7552.v`](testcircuit/c7552.v)

```bash
make -C third_party/abc
cargo run --release --bin generate_circuit -- testcircuit/c7552.v HomCircuit/c7552 --with-merge
```

## Executing a Preprocessed Circuit

```bash
cargo run --release --bin execute_circuit -- <preprocessed_dir> [input_bits] [--single-lut-mode|--hybrid-lut-mode|--multi-lut-mode] [--repeat-random N]
```

Execution modes:

- `--single-lut-mode`
  - only single-output HomLUT
- `--multi-lut-mode`
  - only multi-output HomLUT
- `--hybrid-lut-mode`
  - both single-output and multi-output HomLUT

Examples:

```bash
cargo run --release --bin execute_circuit -- HomCircuit/c7552 --hybrid-lut-mode --repeat-random 1
```

```bash
cargo run --release --bin execute_circuit -- HomCircuit/divisor --hybrid-lut-mode --repeat-random 1
```

## Example Workflow

```bash

cargo run --release --bin generate_circuit -- testcircuit/c7552.v HomCircuit/c7552 --with-merge
cargo run --release --bin execute_circuit -- HomCircuit/c7552 --hybrid-lut-mode --repeat-random 1
```

## Integer and Application Benchmarks

Build all helper binaries first:

```bash
cargo build --release --bins
```

Run all eight unsigned integer operations at 4, 8, 16, 32, 64, 128, and 256 bits:

```bash
./target/release/bench_integer_ops --mode hybrid --repeat 1
```

Run the 64-bit private ERC-20 transfer and the 16-bit TPC-H Query 6 circuit:

```bash
./target/release/bench_erc20 --mode hybrid --repeat 1
./target/release/bench_tpch_q6 --rows 1 --mode hybrid --repeat 1
```

Useful benchmark options:

- `--emit-only`: only write generated gate-level BLIF; no ABC/SALUS execution.
- `--prepare-only`: generate and preprocess the circuit, but do not execute FHE.
- `--out-dir DIR`: choose the generated source and circuit directory.
- `--seed N`: make encryption/test inputs reproducible.
- Integer benchmark only: `--widths 16,64 --ops add,lt` selects a subset.
- Q6 only: `--rows N` controls the number of lineitem rows aggregated.

Integer arithmetic is modulo `2^N`; comparisons are unsigned. Q6 encodes dates
as day numbers, one year as 365, and `0.01` as one fixed-point discount unit.
The integer benchmark also writes its summary to `<out-dir>/results.csv`.

All three benchmark executables use the shared cryptographic configuration in
`src/benchmark_params.rs`: LWE dimension 920, GLWE dimension 1, polynomial size
2048, PBS decomposition 22/1, KS decomposition 3/5, message/carry modulus 4,
maximum noise level 5, native `2^64` ciphertext modulus, and Big-key input
encryption. The requested TUniform(45)/TUniform(17) distributions are represented
by equal-standard-deviation Gaussian noise because the pinned `tfhe = 0.5.4`
`IntLheParam` interface predates TUniform support. SALUS-specific automorphism,
scheme-switching, and circuit-bootstrap decomposition values remain route-specific.
