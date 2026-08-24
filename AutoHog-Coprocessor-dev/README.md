# AutoHoG README

## Requirements

To set up the environment for this project, you will need the following:

```bash
git 
clang
cmake >= 3.16
python >= 3.10

```

## Building AutoHoG

You can build the AutoHoG for your machine by executing the following commands:

```bash
cd AutoHog-Coprocessor-dev
rm -rf build
mkdir build && cd build

CC=/usr/bin/clang-14 \
CXX=/usr/bin/clang++-14 \
AR=/usr/bin/llvm-ar-14 \
RANLIB=/usr/bin/llvm-ranlib-14 \
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_FLAGS_INIT="-fuse-ld=lld" \
      -DCMAKE_CXX_FLAGS_INIT="-fuse-ld=lld" \
      -DCMAKE_EXE_LINKER_FLAGS_INIT="-fuse-ld=lld" \
      -DCMAKE_SHARED_LINKER_FLAGS_INIT="-fuse-ld=lld" \
      ..

cmake --build . -- -j"$(nproc)"
cd ..

```

Supported backend/curve combinations:

- blst/bls12_381: `-DAUTOHOG_ZK_BACKEND=blst -DAUTOHOG_ZK_CURVE=bls12_381` (default)
- mcl/bls12_381: `-DAUTOHOG_ZK_BACKEND=mcl -DAUTOHOG_ZK_CURVE=bls12_381`
- mcl/bn_snark1: `-DAUTOHOG_ZK_BACKEND=mcl -DAUTOHOG_ZK_CURVE=bn_snark1`

When switching to a different backend or curve, please delete `build/CMakeCache.txt` before rebuilding.

## Operation-Examples Test

Operator tests are located in `./Coprocessor-operator`.

---


### **Automated Testing Framework**

The project includes a robust testing suite to verify the **correctness** and **performance** (benchmarking) of all FHE operators.

#### **1. CLI Test Tool (`main_testop.cpp`)**

The `main_testop` binary is the primary entry point for manual and automated testing. It handles key generation, data initialization, and invokes specific test modules for each operator.

* **Supported Operators**:
`FheAdd`, `FheSub`, `FheBitAnd`, `FheBitOr`, `FheBitXor`, `FheBitNot`, `FheEq`, `FheNe`, `FheGe`, `FheGt`, `FheLe`, `FheLt`, `FheIfThenElse`, `TrivialEncrypt`.
* **Key Features**:
* **Dual Mode**: Supports both `std` (standard LWE with noise) and `triv` (deterministic trivial encryption) modes.
* **Performance Benchmarking**: Compares **Serial** vs. **Parallel** (OpenMP) execution times.
* **Validation**: Automatically generates plaintext ground truth using `boost::multiprecision` to verify FHE results.



**Usage**:

```bash
./main_testop <op> <bit_width> <mode> [epochs] [detail]

```

* `op`: The operator name (e.g., `FheAdd`).
* `bit_width`: Operand size in bits (must be a multiple of 8, e.g., `32`, `64`).
* `mode`: `std` for standard FHE, `triv` for trivial encryption.
* `epochs`: (Optional) Number of iterations to average the performance results. Default is `1`.
* `detail`: (Optional) `1` to show per-epoch details, `0` for summary only. Default is `1`.

---

#### **2. Automation Script (`run_optests.sh`)**

This Bash script automates the execution of the entire operator suite. It is designed for regression testing and generating full performance reports.

* **Log Management**: Automatically creates a `./logs` directory and saves execution results into timestamped files (e.g., `fhe_test_20260315_1530.log`).
* **Batch Execution**: Iterates through all 14 supported operators in a single run.
* **Configurable**: Easily adjust `BIT_WIDTH`, `MODE`, and `EPOCHS` within the script header.

**Usage**:

```bash
# 1. Grant execution permission
chmod +x Coprocessor-operator/run_optests.sh

# 2. Run the full test suite
./Coprocessor-operator/run_optests.sh

```

---

#### **3. Example Output (Summary)**

The framework provides a structured summary at the end of each operator test:

```text
================= SUMMARY (FheAdd) =================
Data Size     : 8 Bytes (64 bits)
HomAdd (ser)  : 939.983 ms (avg)
HomAdd (pal)  : 45.624 ms (avg)
Speedup       : x20.58
----------------------------------------------------
Serial vs Plain 2^W : 100/100 (100%)
Parallel vs Serial  : 100/100 (100%)
====================================================

```

> [!NOTE]
> For comparison operators (e.g., `FheEq`, `FheGt`), the test suite automatically evaluates three distinct scenarios per epoch: **random values**, **guaranteed equal**, and **guaranteed unequal** to ensure full circuit coverage.

---


### **Encryption and Decryption**

* **Plaintext operand**:

```cpp
std::vector<uint8_t> Apl(W); 
// W-bit binary sequence in two’s complement (signed).
// Apl[0] is the MSB, Apl[W-1] is the LSB.

```

* **Ciphertext operand**:

```cpp
std::vector<TLWE<lvl2param>> A(W); 
// W-bit LWE vector.

```

* **Encryption**:

```cpp
for (int i = 0; i < W; ++i)
    A[i] = bootsSymEncrypt<lvl2param>(Apl[i], sk);

```

* **Decryption**:

```cpp
for (int i = 0; i < W; ++i)
    Apl[i] = bootsSymDecrypt<lvl2param>(A[i], sk);

```

* **Key generation**:

```cpp
SecretKey sk;
EvalKey ek;
ek.emplacebkfft<lvl02param>(sk);
ek.emplaceiksk<lvl20param>(sk);

```

### **Trivial Encryption**

Trivial encryption is used to encode **public values** (constants) into valid FHE ciphertext formats without introducing randomness or noise. Unlike standard encryption, trivial encryption is **deterministic**, meaning the same plaintext input will always produce the exact same ciphertext binary. This is essential for consensus scenarios where multiple nodes must initialize state (e.g., zero) identically.

* **Single Value**:

```cpp
// Encrypts a single byte/bit deterministically (Mask=0, Noise=0)
auto ct = key.trivialEncrypt(plain_val);

```

* **Batch Encryption**:

```cpp
// Encrypts a vector of bytes deterministically
auto vec = key.trivialEncryptBatch(plain_bytes);

```

* **Test Example**:
To verify correctness and determinism, run the provided test script:

```bash
# Usage: ./Coprocessor-operator/test_trivial <byte_count> [epochs]
./Coprocessor-operator/test_trivial 100 5

```

Here are the additions to your documentation, formatted to match your existing style.

### **Trivial Operations & Conversion**

Beyond basic encryption, the system supports **deterministic arithmetic** and **state conversion** on trivial ciphertexts. These operators ensure that computations on public values across different nodes yield bit-identical results (Essential for Consensus).

* **Deterministic Arithmetic**:
Performs logic-gate based arithmetic (Ripple Carry) on trivial ciphertexts without introducing noise.

```cpp
// Deterministic Addition (LHS + RHS)
auto res_add = key.trivialAdd(lhs_ct, rhs_ct);

// Deterministic Subtraction (LHS - RHS)
auto res_sub = key.trivialSub(lhs_ct, rhs_ct);

```

* **Privacy & Consensus Conversion**:
* **Trivial -> Standard**: Injects random noise into a trivial ciphertext to make it secure (indistinguishable from random) for transmission.
* **Standard -> Trivial**: "Cleans" a noisy standard ciphertext back into a deterministic trivial format ( Consensus Cleaning).



```cpp
// Inject Noise: Output is randomized (Secure)
auto std_ct = key.trivialToStandard(triv_ct);

// Clean Noise: Output is deterministic (Consensus)
auto clean_ct = key.standardToTrivial(std_ct);

```

* **Test Example**:
To verify the determinism of arithmetic operators and the full Privacy-Consensus lifecycle, run:

```bash
# Verifies Add/Sub determinism and Trivial<->Standard conversion
./Coprocessor-operator/test_trivial_ops.cpp

```

---

### **AHK Persistence & Packing**

The **AnnihilateKey (AHK)** allows a server to compress multiple LWE ciphertexts into a single TRLWE ciphertext **without** access to the Secret Key. This is critical for optimizing downstream bandwidth in Client-Server architectures.

* **Key Persistence**:
The AHK can be generated on the client side, saved to a file, and loaded on the server side independently.

```cpp
// Save AHK to a specific file path
key.save(sk_path, eval_path, "server_ahk.bin", pk_path, pk_v2_path);

// Load AHK on server (Secret Key not required if only packing is needed)
key.load(sk_path, eval_path, "server_ahk.bin", pk_path, pk_v2_path);

```

* **Packing & Unpacking**:

```cpp
// Packing (LWEs -> TRLWE): Compresses data using AHK
auto packed_trlwe = key.encryptData(raw_bits); 

// Unpacking (TRLWE -> LWEs): Extracts individual ciphertexts
auto unpacked_tlwes = decompress_ciphertexts(packed_trlwe);

```

* **Test Example**:
To verify that the AnnihilateKey can be saved/loaded identically and that packing/unpacking preserves data integrity:

```bash
# Verifies AHK persistence and Pack/Unpack correctness
./Coprocessor-operator/test_pack_unpack

```

### **HomADD & HomSUB**

```cpp
HomAddResult HomAdd(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    bool compute_overflow
);

HomAddResult HomSub(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    bool compute_overflow
);

struct HomAddResult {
    std::vector<TLWE<lvl2param>> SUM;              // W-bit result
    std::optional<TLWE<lvl2param>> ovf_ct;         // overflow flag (optional)
};

```

* **INPUT**:
* A_in, B_in: W-bit ciphertext vectors (LWE arrays)
* W: bit-width
* ek: evaluation key
* compute_overflow: whether to compute the overflow flag


* **OUTPUT**:
* SUM: W-bit result (LWE array)
* ovf_ct: overflow indicator (LWE), optional


### **HomMUL & HomDIV**

`HomMul` and `HomDiv` implement encrypted integer multiplication and division over
MSB-first `std::vector<TLWE<lvl2param>>` operands.

```cpp
std::vector<TLWE<lvl2param>> HomMul(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek
);

std::vector<TLWE<lvl2param>> HomMul_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    int threads_hint = -1
);

HomDivResult HomDiv(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek
);

HomDivResult HomDiv_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    int threads_hint = -1
);

struct HomDivResult {
    std::vector<TLWE<lvl2param>> quotient;   // W-bit quotient
    std::vector<TLWE<lvl2param>> remainder;  // W-bit remainder
};

```

* **INPUT**:
* A_in, B_in: W-bit ciphertext vectors (LWE arrays), MSB first
* W: input operand bit-width
* ek: evaluation key
* threads_hint: optional OpenMP thread hint for parallel variants


* **HomMul OUTPUT**:
* Returns a `2W`-bit ciphertext vector containing the full unsigned product.
* Internally, both operands are first cast from `W` bits to `2W` bits by zero-extending the MSB side with `CT_ZERO`, then multiplication is executed on the widened operands.
* Example: multiplying two 64-bit encrypted integers returns a 128-bit encrypted product.

```cpp
auto product128 = HomMul(A64, B64, 64, ek);
// product128.size() == 128

```

* **HomDiv OUTPUT**:
* Returns a `HomDivResult` containing a W-bit quotient and a W-bit remainder.
* For `B_in != 0`, the result satisfies `A = quotient * B + remainder` and `remainder < B`.
* For division by zero, the implementation returns `quotient = 0` and `remainder = A`.

```cpp
auto div_res = HomDiv(A, B, W, ek);
auto quotient = div_res.quotient;
auto remainder = div_res.remainder;

```

* **Test Example**:
To verify multiplication/division correctness with per-epoch decrypted detail:

```bash
# Multiplication only: full 2W-bit product semantics
./build/Coprocessor-operator/test_MULDIV mul 4 1 --show-epochs

# Division only, including a division-by-zero epoch when epochs >= 1
./build/Coprocessor-operator/test_MULDIV div 4 2 --show-epochs

# Both operators
./build/Coprocessor-operator/test_MULDIV both 4 1 --show-epochs

# Summary-only mode
./build/Coprocessor-operator/test_MULDIV both 4 10 --no-epoch-detail

```

In detail mode, each epoch prints the input bits, output bits, expected plaintext value, and decrypted plaintext value. Multiplication reports `input_width=<W>` and `output_width=<2W>` to make the full-product bit width explicit. Division reports W-bit quotient and W-bit remainder; for division by zero, the expected behavior is `quotient = 0` and `remainder = A`.



### **HomComp**

```cpp
TLWE<lvl2param> HomLE(A_in, B_in, W, ek);  // A <= B
TLWE<lvl2param> HomLT(A_in, B_in, W, ek);  // A <  B
TLWE<lvl2param> HomGE(A_in, B_in, W, ek);  // A >= B
TLWE<lvl2param> HomGT(A_in, B_in, W, ek);  // A >  B
TLWE<lvl2param> HomEQ(A_in, B_in, W, ek);  // A == B

```

* **INPUT**:
* A_in, B_in: W-bit ciphertext vectors (LWE arrays)
* W: bit-width
* ek: evaluation key


* **OUTPUT**:
* A single LWE ciphertext representing the comparison result (true / false).



## PBS-Error Correction

This project provides two lightweight noise-correction techniques for **lvl0 TLWE after IKS**
(i.e., **after `IdentityKeySwitch` and before `BlindRotate**`):

* **Drift**: mitigates systematic bias and heavy-tailed noise introduced by modulus switching (MS); **currently achieves ~128-bit error security**.
* **Centered**: recenters the phase around ±μ to reduce the risk of hitting the decision boundary; **currently achieves ~122-bit error security**.

```cpp
#include "drift_ms_tfhepp.hpp"
#include "drift_ms_params.hpp"

```

### **1. Drift key (if Drift is enabled)**

Drift requires an auxiliary **lvl0 key** consisting of several encryptions of zero:

```cpp
using DriftKeyLvl0 = driftms::DriftMsKey<lvl0param>;

DriftKeyLvl0 drift_key_lvl0;
{
    using TestParam = DriftMsTestParamLvl0;
    using P0        = lvl0param;

    drift_key_lvl0.r_sigma        = TestParam::r_sigma;
    drift_key_lvl0.bound          = TestParam::bound();
    drift_key_lvl0.input_variance = TestParam::input_variance();
    drift_key_lvl0.log_mod        = TestParam::log_mod;

    drift_key_lvl0.zeros.resize(TestParam::zeros_count);
    for (auto &z : drift_key_lvl0.zeros)
        z = tlweSymEncrypt<P0>(0, P0::α, sk.key.get<P0>());
}

```

### **2. Where to apply Drift / Centered**

**Fixed position: after IKS, before PBS**

```cpp
IdentityKeySwitch<lvl20param>(tlwe0, linear, ksk);

// Drift (optional)
driftms::improve_lwe_ciphertext_modulus_switch_noise_for_binary_key<lvl0param>(
    tlwe0, drift_key_lvl0
);

// Centered (optional)
centered::apply_centered_correction_to_body<lvl0param>(
    tlwe0, drift_key_lvl0.log_mod
);

// Then proceed to BlindRotate(...)

```

### **3. Example**

```bash
# code: Coprocessor-operator/and-chain-pack.cpp
./build/Coprocessor-operator/and-chain-pack

```

## HomDecompress (Packed TRLWE Decompression)

`HomDecompress` provides **three decompression pipelines** for converting a **packed TRLWE at lvl1** into a **packed TRLWE at lvl2**.

All three follow the same high-level pattern:

1. `TRLWE_Decompress` : packed TRLWE → `slot_count` TLWE slots
2. per-slot transform (parallelized with OpenMP)
3. `TLWE2TRLWEPacking` : TLWE slots → packed TRLWE (lvl2)

* **Packed: TRLWE<N=2048> (lvl2): 32.000 KB**
* **Compressed: TRLWE<N=1024> (lvl1): 8.000 KB**

> Fixed params: `Pin = lvl1param`, `Pout = lvl2param`.

### 1. Usage

**（1）PBS + IKS**

```cpp
double Decompress_PBSiks_lvl1_to_lvl2(
    TFHEpp::TRLWE<Pout>                                  &packed_out,
    const TFHEpp::TRLWE<Pin>                             &packed_in,
    const TFHEpp::KeySwitchingKey<TFHEpp::lvl10param>     &ksk10,     // 1 -> 0
    const TFHEpp::BootstrappingKeyFFT<TFHEpp::lvl02param> &bkfft,     // 0 -> 2
    const TFHEpp::AnnihilateKey<Pout>                     &ahk2,
    int                                                    slot_count
);

```

**Pipeline**

* Unpack: TRLWE(lvl1) → TLWE(lvl1)[S]
* Per slot (OMP): IKS(1->0) then PBS(0->2)
* Repack: TLWE(lvl2)[S] → TRLWE(lvl2)

**（2） PBS-only**

```cpp
double Decompress_PBS_lvl1_to_lvl2(
    TFHEpp::TRLWE<Pout>                                  &packed_out,
    const TFHEpp::TRLWE<Pin>                             &packed_in,
    const TFHEpp::BootstrappingKeyFFT<TFHEpp::lvl12param> &bkfft12,   // 1 -> 2
    const TFHEpp::AnnihilateKey<Pout>                     &ahk2,
    int                                                    slot_count
);

```

**Pipeline**

* Unpack: TRLWE(lvl1) → TLWE(lvl1)[S]
* Per slot (OMP): PBS(1->2) (requires bkfft12)
* Repack: TLWE(lvl2)[S] → TRLWE(lvl2)

**（3） Direct IKS**

```cpp
double Decompress_DirectIKS_lvl1_to_lvl2(
    TFHEpp::TRLWE<Pout>                               &packed_out,
    const TFHEpp::TRLWE<Pin>                          &packed_in,
    const TFHEpp::KeySwitchingKey<TFHEpp::lvl12param>  &ksk12,     // 1 -> 2
    const TFHEpp::AnnihilateKey<Pout>                  &ahk2,
    int                                                 slot_count,
    const DriftKeyLvl2                                 *drift_key_lvl2 = nullptr
);

```

**Pipeline**

* Unpack: TRLWE(lvl1) → TLWE(lvl1)[S]
* Per slot (OMP): IKS(1->2)
* (Optional) post-fix hooks (e.g., centered / drift on lvl2 TLWE)
* Repack: TLWE(lvl2)[S] → TRLWE(lvl2)

### 2. Examples

```bash
# code：Coprocessor-operator/Decompressed.cpp
./build/Coprocessor-operator/Decompressed

```

## FHE Consensus & LSB Masking

When concurrently invoking the underlying floating-point FFT engine in a multi-threaded/multi-processed environment, physical-level rounding errors of 1~2 bits are generated in the least significant bits (LSB) of the polynomial coefficients when transforming from the frequency domain back to the time domain. This is due to the dynamic scheduling of FMA floating-point instructions by the hardware, microarchitectural differences in heterogeneous cores (P-cores/E-cores), and thread context switching.

**Solution**: In both the application layer and the underlying circuits, forcefully truncate the pure noise region (smoothing out LSB errors) via bitwise operations. This trades a negligible constant noise overhead for physical-level consensus alignment.

------

### 1. Core Modification Guidelines

The modifications consist of two lines of defense: **Deep Circuit Truncation** and **Top-Level Exit Truncation**.

### （1） Surgical Masking

**Applicable Scenarios:** Deep logic files containing multi-layer nested PBS, such as `src/HomComp_unsigned.cpp`.

**Development Rules:** In any operator containing carry chains or tree reductions, **the sanitization function must be called immediately after every PBS execution**. This kills the error at the current layer, preventing it from flowing into the next layer and triggering an avalanche effect.

First, define the generic sanitization function at the top of the corresponding source file (within an anonymous namespace):

C++

```
namespace {
    using namespace TFHEpp;

    // 🌟 Generic LSB Sanitizer: Forcefully smooths out the last 4 bits of floating-point glitches
    template <typename TLWE_TYPE>
    inline void sanitize_tlwe(TLWE_TYPE& ct) {
        for (auto& val : ct) {
            val &= ~0xFU; 
        }
    }
}
```

Then, apply this sanitizer immediately after the PBS calls in each Gadget:

C++

```
// [Example 1: Comparison Gadget for parallel tree nodes]
static void GT2_unsigned_PBS(std::vector<TLWE<lvl2param>> &out, ...) {
    // ... Preceding linear combination logic ...
    PBS_Multi(out, linear, ek, lut);
    
    // 🌟 Core Interception: Kill the computation error of the current leaf node
    for (auto& c : out) sanitize_tlwe(c); 
}

// [Example 2: Serial carry chain Gadget]
static void HomFA3_carry_PBS(TLWE<lvl2param> &carry_out, ...) {
    // ... Preceding linear combination logic ...
    PBS_Single(carry_out, linear, ek, lut);
    
    // 🌟 Core Interception: Directly sanitize the single ciphertext
    sanitize_tlwe(carry_out);
}
```

### （2）Top-Level Masking

**Applicable Scenarios:** Top-Level dispatch interfaces that output results externally, such as `Coprocessor-operator/test_consensus.cpp`.

**Development Rules:** For simple operators with only a single-layer depth (e.g., ADD, SUB) or deep copy operators (e.g., OP=16), execute a non-discriminatory, catch-all sanitization right before `return` to ensure that the final ciphertext written to disk or transmitted over the network is 100% binary consistent.

C++

```
tlwe_vec_ptr compute(int op, const vector<tlwe_vec_ptr> &cts, FheKey &key) {
    tlwe_vec_ptr res;
    
    // ... [Execute various OP operator branch logic] ...
    
    // 🌟 [Consensus Lock]: Non-discriminatorily smooth out the last 4 bits of floating-point error before returning the final ciphertext
    if (res) { 
        // Note: Must use res->data to extract the actual array inside the wrapper
        for (auto& ct : res->data) {
            for (auto& val : ct) {
                val &= ~0xFU; 
            }
        }
    }
    return res;
}
```

------

### 2. Developer Checklist (Must Read)

When writing or extending new FHE operators, team members must verify the following two items:

- [ ] **Does it contain multi-layer nested PBS?** If developing a new multiplier or complex logic circuit, you must insert `sanitize_tlwe()` after every intermediate `PBS_Multi` or `PBS_Single` in the corresponding `HomOP.cpp` file.
- [ ] **Does it bypass the top-level sanitization?** Ensure that all ciphertext references exposed to the consensus layer have undergone the `val &= ~0xFU` mask filtering before leaving `compute` (or the related top-level interface).

### 4. Concurrency Stress Testing Guide

To verify the effectiveness of the LSB truncation and ensure that the binary consensus of node computation results is 100% aligned under high multi-threading concurrency, a one-click stress testing script is provided.

The script (e.g., `test_concurrency.sh`) will automatically:

1. Clean the execution environment (deleting old keys and temporary directories `tmp_*`).
2. Launch Coprocessor nodes to perform stress tests sequentially according to preset concurrency levels (default: `2, 4, 8, 16, 32`).
3. Capture all standard/error outputs and isolate the detailed execution logs of each round into the `test_logs/` directory.
4. Monitor the execution status, halting immediately if a crash or hash inconsistency occurs.



## ZK Workflow (RLWE Ciphertext Proof)

---

### 1. Components

#### 1.1 Trusted Setup

- Obtain a BLS12-381 CRS from an audited multi-party ceremony and convert it to the native CRS format
- Outputs:
  - `crs_lvl2.bin`
  - corresponding `SHA-256 digest`
- Independently verify the ceremony contribution chain and pin the final digest through a trusted channel
- `FheKey::keygen()` never generates a CRS; client and server must load the external artifact explicitly
- The bundled `gen_crs` is an insecure single-party test tool and must not be used for production

---

#### 1.2 Client / Prover

- Holds the public key
- Loads CRS and a trusted CRS SHA-256 digest
- Executes:

    encrypt + prove

- Sends:

    ciphertext + proof + metadata

---

#### 1.3 Server / Verifier

- Holds public key + CRS  
- Executes:

    verify (ZK proof)

- If verification passes:

    decrypt / circuit evaluation

---

### 2. Workflow

1. Multi-party Trusted Setup:

    obtain and independently verify a compatible BLS12-381 ceremony artifact

2. Save:

    crs_lvl2.bin

3. Record:

    SHA-256 digest

4. Client / Server load CRS and provide a trusted SHA-256 digest

5. Client:

    encrypt + prove

6. Server:

    verify

7. If valid:

    decrypt / circuit evaluation

---

### 3. Example

#### 3.1 Test-only single-party setup

    gen_crs crs_lvl2.bin --allow-insecure-single-party
      -> ZKProof::generate_public_params
      -> ZKProof::save_public_params
      -> ZKProof::compute_file_sha256_hex

This command is only for local regression tests. Without the explicit acknowledgement flag,
`gen_crs` refuses to create a single-party CRS. Production must use a compatible, audited
multi-party ceremony artifact; the original audit's BN254/Hermez `.ptau` file is not compatible
with the current BLST + BLS12-381 backend.

---

#### 3.2 Client

    FheKey::load(..., crs_file, crs_expected_sha256)
      -> ZKProof::load_public_params_by_expected_hash

    FheKey::encryptDataWithZk
      -> ZKPKE::encrypt_bits_with_zk
      -> ZKPKE::encrypt_with_zk_proof
      -> ZKProof::make_public_commit
      -> ZKProof::prove
      -> serialize_ciphertext_with_proof

---

#### 3.3 Server

    FheKey::load(..., crs_file, crs_expected_sha256)
      -> ZKProof::load_public_params_by_expected_hash

    FheKey::verifyZkProof
      -> deserialize_ciphertext_with_proof
      -> ZKProof::make_public_commit
      -> ZKProof::verify

or:

    loadSealedCiphertext
    -> deserialize ciphertext / proof / metadata

    then
    -> ZKPKE::verify_and_decompress
       -> ZKProof::verify
       -> TRLWE_Decompress

---

### 4. CRS Digest Verification

The CRS must be verified using a trusted digest. Do NOT rely on local `.sha256` files.

#### Recommended Approach

- Multi-party setup:
  - obtain and verify a compatible BLS12-381 `crs_lvl2.bin`
  - compute SHA-256 digest
  - distribute via trusted channel  

- Runtime:
  - provide expected digest
  - verify local CRS  

    actual_hash == expected_hash

For the C ABI, use the authenticated entry point:

```cpp
create_fhe_key_with_zk_sha256(
    pk, pk_len,
    crs, crs_len,
    trusted_sha256_hex, 64);
```

The legacy `create_fhe_key_with_zk` symbol is retained for binary compatibility but fails
closed. Likewise, `FheKey::load(..., crs_file)` and `FheKey::load_crs(crs, len)` reject a
non-empty CRS unless the caller supplies an independently trusted SHA-256 digest.

---

#### API

    ZKProof::load_public_params_by_expected_hash

---

#### Example
    AUTOHOG_CRS_SHA256=<trusted_digest> ./test_with_zk --crs crs_lvl2.bin

    ./test_with_zk --crs crs_lvl2.bin --crs-sha256 <trusted_digest>

---
