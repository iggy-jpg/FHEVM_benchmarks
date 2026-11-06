# FHE VM Benchmarks

This repository contains the benchmark implementations used to compare FHE
virtual-machine programming models. The benchmark scope follows the evaluation
in the accompanying paper: eight unsigned integer instructions and two
end-to-end programs.

## Benchmark scope

The instruction matrix contains:

- arithmetic: `ADD`, `SUB`, and `MUL`;
- comparisons: `GT`, `LT`, `GE`, `LE`, and `EQ`;
- precisions: 4, 8, 16, 32, 64, 128, and 256 bits where supported by the
  implementation.

The application benchmarks are:

- **P1 / ERC-20**: three 64-bit private token transfers;
- **P2 / TPC-H Q6**: a Q6-style filtered revenue aggregation with 16-bit input
  fields. The paper evaluates 32,768 rows; several executables accept a smaller
  row count for smoke tests.

The paper reports single-thread computation. Use `OMP_NUM_THREADS=1` and
`RAYON_NUM_THREADS=1` for comparable CPU latency runs. Key generation,
encryption, decryption, and setup are excluded only when the individual
benchmark says so; check each program's output before comparing numbers.

## Repository layout

| VM model | Implementation | Main benchmark paths |
| --- | --- | --- |
| SISD-Bool | AutoHog | `AutoHog-Coprocessor-dev/Coprocessor-operator/test_PERF_CORE.cpp`, `program_erc20_64.cpp`, `program_tpch_q6_16.cpp` |
| SISD-Carry | TFHE-rs | `tfhe-rs/tfhe-benchmark/benches/integer/bench.rs`, `tfhe-rs/tfhe/examples/erc20_64.rs`, `tpch_q6_16.rs` |
| SISD-Tree | SALUS | `SALUS-main/src/bin/bench_integer_ops.rs`, `bench_erc20.rs`, `bench_tpch_q6.rs` |
| SIMD-Poly | OpenFHE | `openfhe-development/benchmark/src/openfhe_requested_integer_benchmark.cpp` |
| SIMD-Bool | CKKS Boolean circuits | `ks256/src/boolean_operator_benchmark_cpu.cpp`, `transaction_64bit_cpu.cpp`, `tpch_q6_16bit_cpu.cpp` |
| SIMD-Carry | SEC-CKKS / Engorgio | `sec-ckks/src/pke/examples/carry_test.cpp`, `Engorgio/test/operation_performance_test.cpp` |
| Hybrid | ArcEDB | `ArcEDB/test/integer_perf_matrix_test.cpp`, `hybrid_transaction_test.cpp`, `relational_queries_test.cpp` |

All commands below assume Ubuntu Linux and start in the repository root.

## Common prerequisites

Install a C/C++ toolchain, CMake, Ninja or Make, OpenMP, GMP/NTL, and Rust:

```bash
sudo apt update
sudo apt install -y build-essential cmake ninja-build clang-14 lld-14 \
    libomp-14-dev libgmp-dev libntl-dev libssl-dev pkg-config python3
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
```

Some projects vendor large third-party source trees. A fresh Git clone should
therefore be obtained recursively when submodule metadata is present:

```bash
git clone --recursive <repository-url>
```

## 1. AutoHog (SISD-Bool)

Source directory: `AutoHog-Coprocessor-dev/`.

The default ZK backend requires a complete BLST tree at
`thirdparties/blst/`, including `bindings/`. Configure and build with:

```bash
cd AutoHog-Coprocessor-dev
CC=/usr/bin/clang-14 \
CXX=/usr/bin/clang++-14 \
AR=/usr/bin/llvm-ar-14 \
RANLIB=/usr/bin/llvm-ranlib-14 \
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_FLAGS_INIT="-fuse-ld=lld" \
    -DCMAKE_CXX_FLAGS_INIT="-fuse-ld=lld" \
    -DCMAKE_EXE_LINKER_FLAGS_INIT="-fuse-ld=lld" \
    -DCMAKE_SHARED_LINKER_FLAGS_INIT="-fuse-ld=lld"
cmake --build build --target test_PERF_CORE program_erc20_64 program_tpch_q6_16 -j"$(nproc)"
```

Run the complete instruction matrix and the two applications serially:

```bash
OMP_NUM_THREADS=1 ./build/Coprocessor-operator/test_PERF_CORE 3 4,8,16,32,64,128,256
OMP_NUM_THREADS=1 ./build/Coprocessor-operator/program_erc20_64 3
OMP_NUM_THREADS=1 ./build/Coprocessor-operator/program_tpch_q6_16 32768
```

`test_PERF_ALL <width> <epochs> <op>` is also available. Its operator names
include `add_u`, `sub_u`, `mul`, `gt_u`, `lt_u`, `eq`, `ge_u`, and `le_u`.

## 2. TFHE-rs (SISD-Carry)

Source directory: `tfhe-rs/`. The workspace requires the Rust version declared
in `Cargo.toml` (currently 1.91.1).

```bash
cd tfhe-rs
rustup toolchain install 1.91.1
rustup override set 1.91.1
cargo build --release -p tfhe --features integer --examples
cargo build --release -p tfhe-benchmark \
    --features "integer internal-keycache pbs-stats"
```

Run the CPU integer Criterion benchmark with one Rayon worker, or filter it by
a Criterion benchmark name after the final `--`:

```bash
RAYON_NUM_THREADS=1 cargo bench -p tfhe-benchmark --bench integer \
    --features "integer internal-keycache pbs-stats"
```

Run the two application examples:

```bash
RAYON_NUM_THREADS=1 cargo run --release -p tfhe --example erc20_64 --features integer
RAYON_NUM_THREADS=1 cargo run --release -p tfhe --example tpch_q6_16 --features integer
```

## 3. SALUS (SISD-Tree)

Source directory: `SALUS-main/`. SALUS vendors ABC and a Rust TFHE dependency.

```bash
cd SALUS-main
make -C third_party/abc -j"$(nproc)"
cargo build --release --bins
```

Run the eight instructions and application benchmarks:

```bash
RAYON_NUM_THREADS=1 ./target/release/bench_integer_ops --mode hybrid --repeat 3
RAYON_NUM_THREADS=1 ./target/release/bench_erc20 --mode hybrid --repeat 3
RAYON_NUM_THREADS=1 ./target/release/bench_tpch_q6 --rows 32768 --mode hybrid --repeat 3
```

For a quick check, use `--rows 1` for Q6 or select a subset such as
`--widths 16,64 --ops add,lt` for the instruction benchmark. See
`SALUS-main/README.md` for circuit generation and the single-, hybrid-, and
multi-LUT execution modes.

## 4. OpenFHE (SIMD-Poly)

Source directory: `openfhe-development/`.

```bash
cd openfhe-development
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_BENCHMARKS=ON -DBUILD_EXAMPLES=ON -DWITH_OPENMP=OFF
cmake --build build --target openfhe_requested_integer_benchmark -j"$(nproc)"
OMP_NUM_THREADS=1 ./build/bin/benchmark/openfhe_requested_integer_benchmark
```

This executable covers the precision ranges supported by its polynomial
configuration: ADD at 4--32 bits, MUL at 4--16 bits, and EQ/GT at 4--8 bits.

## 5. CKKS Boolean circuits (SIMD-Bool)

Source directory: `ks256/`. It links against an installed OpenFHE package.
Install the bundled OpenFHE build first:

```bash
cd openfhe-development
cmake -S . -B build-install -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED=ON -DWITH_OPENMP=OFF
cmake --build build-install -j"$(nproc)"
sudo cmake --install build-install
cd ../ks256
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
```

Run the instruction matrix and applications:

```bash
OMP_NUM_THREADS=1 ./build/boolean_operator_benchmark_cpu
OMP_NUM_THREADS=1 ./build/transaction_64bit_cpu
OMP_NUM_THREADS=1 ./build/tpch_q6_16bit_cpu 32768
```

The instruction program writes a copy of its output to
`boolean_operator_benchmark_cpu.log`.

## 6. SEC-CKKS and Engorgio (SIMD-Carry)

### SEC-CKKS

Source directory: `sec-ckks/`.

```bash
cd sec-ckks
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_EXAMPLES=ON -DBUILD_BENCHMARKS=ON -DWITH_OPENMP=OFF
cmake --build build --target carry_test -j"$(nproc)"
```

`carry_test` accepts `<total_bits> <segment_bits> <add|mul>`. For example:

```bash
OMP_NUM_THREADS=1 ./build/bin/examples/pke/carry_test 32 4 add
OMP_NUM_THREADS=1 ./build/bin/examples/pke/carry_test 32 4 mul
CARRY_JOBS=1 ./run_carry_4bit_benchmarks.sh 8 16 32 64 128 256
```

### Engorgio

Source directory: `Engorgio/`. Its OpenFHE fork is included under
`Engorgio/thirdparty/openfhe-development/`.

```bash
cd Engorgio
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DWITH_OPENMP=OFF
cmake --build build --target operation_performance_test -j"$(nproc)"
OMP_NUM_THREADS=1 ./build/bin/operation_performance_test
```

The current operation benchmark measures packed ADD, MUL, and GT at the bit
widths implemented in `test/operation_performance_test.cpp`.

## 7. ArcEDB (Hybrid)

Source directory: `ArcEDB/`.

```bash
cd ArcEDB
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    -DSEAL_THROW_ON_TRANSPARENT_CIPHERTEXT=OFF \
    -DSEAL_USE_INTEL_HEXL=ON
cmake --build build --target integer_perf_matrix_test \
    hybrid_transaction_test relational_queries_test -j"$(nproc)"
```

Run the instruction matrix and representative hybrid workloads:

```bash
OMP_NUM_THREADS=1 ./build/bin/integer_perf_matrix_test
OMP_NUM_THREADS=1 ./build/bin/hybrid_transaction_test
OMP_NUM_THREADS=1 ./build/bin/relational_queries_test
```

`integer_perf_matrix_test` evaluates ADD, SUB, MUL, GT, LT, EQ, GE, and LE at
4--256 bits. A single comparison pass produces all five comparison predicates,
so those five rows report a shared comparison latency.

## Reproducibility notes

- Use a Release build and record the compiler, CPU model, thread count, and
  exact Git commit.
- Run a warm-up before collecting measurements and use at least three measured
  repetitions when the executable supports it.
- Do not compare total process time with kernel-only timing. Read the timing
  boundary printed by each executable.
- Large packed-SIMD benchmarks can require substantial RAM and long key
  generation or precomputation time.
- Every benchmark that provides a correctness status should report success
  before its timing is retained.

