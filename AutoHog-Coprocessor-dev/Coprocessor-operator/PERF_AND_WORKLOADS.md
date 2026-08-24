# Core performance matrix and workload programs

The repository already has `test_PERF_ALL`, whose unsigned operator names are
`add_u`, `sub_u`, `mul`, `gt_u`, `lt_u`, `eq`, `ge_u`, and `le_u`. It accepts
any single positive bit width, but its default runner previously covered only
64 and 128 bits and mixed these operations with many unrelated operators.

`test_PERF_CORE` is the focused executable. By default it measures and checks
`add`, `sub`, `mult`, `gt`, `lt`, `eq`, `ge`, and `le` at all requested widths:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target test_PERF_CORE -j
./build/Coprocessor-operator/test_PERF_CORE
```

The optional arguments are the epoch count and a comma-separated width list:

```bash
./build/Coprocessor-operator/test_PERF_CORE 3 4,8,16,32,64,128,256
```

Build and run the two application workloads with:

```bash
cmake --build build --target program_erc20_64 program_tpch_q6_16 -j
./build/Coprocessor-operator/program_erc20_64 2
./build/Coprocessor-operator/program_tpch_q6_16 2
```

The ERC-20 argument is the epoch count. The Q6 argument is the number of input
rows. Timing excludes key generation, encryption, and decryption. Every
executable decrypts its result after timing and returns exit code 2 if a
correctness check fails.

The Q6 input columns and predicates are 16-bit unsigned fixed-point values.
The natural `extended_price * discount` result and aggregate are 32-bit values
modulo 2^32.
