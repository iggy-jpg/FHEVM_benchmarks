#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <string>
#include <array>
#include <filesystem>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <cstdlib>

#include <tfhe++.hpp>
#include "Autohog_Cop.hpp"
#include "ZKPKE.h"
#include "zk_backend.h"
#include "SecureRandom.h"

namespace {

bool is_decimal_arg(const std::string& value)
{
    return !value.empty() &&
           std::all_of(value.begin(), value.end(), [](unsigned char ch) {
               return std::isdigit(ch) != 0;
           });
}

size_t ilog2_u64_local(uint64_t v)
{
    size_t lg = 0;
    while (v > 1) {
        v >>= 1;
        ++lg;
    }
    return lg;
}

void print_usage(const char* argv0)
{
    std::cout
        << "Usage: " << argv0 << " [num_bits] [--crs path] [--crs-sha256 hex]\n"
        << "  num_bits:      payload bits to encrypt/test (default: 128, range: 32..256)\n"
        << "  --crs path:    CRS file path (default: crs_lvl2.bin)\n"
        << "  --crs-sha256:  trusted SHA-256 digest for the CRS file\n"
        << "Environment:\n"
        << "  AUTOHOG_CRS_SHA256  trusted SHA-256 digest for the CRS file\n";
}

} // namespace

int main(int argc, char** argv) try {
    using namespace TFHEpp;
    using P = lvl2param;

    autohog::SecureRng engine;
    std::bernoulli_distribution coin(0.5);

    std::cout << "========================================\n";
    std::cout << "  TFHE with ZK Proof + Bootstrap Test\n";
    std::cout << "========================================\n\n";

    constexpr uint32_t MAX_NUM_BITS = 256;
    static_assert(MAX_NUM_BITS <= P::n, "MAX_NUM_BITS must be <= P::n");
    uint32_t num_bits = 128;
    bool saw_num_bits = false;
    std::string crs_file = "crs_lvl2.bin";
    std::string expected_crs_hash;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        }
        if (arg == "--crs") {
            if (i + 1 >= argc) {
                throw std::runtime_error("--crs requires a file path");
            }
            crs_file = argv[++i];
            continue;
        }
        if (arg == "--crs-sha256") {
            if (i + 1 >= argc) {
                throw std::runtime_error("--crs-sha256 requires a 64-hex digest");
            }
            expected_crs_hash = argv[++i];
            continue;
        }
        if (!saw_num_bits && is_decimal_arg(arg)) {
            unsigned long parsed = std::stoul(arg);
            if (parsed < 32 || parsed > MAX_NUM_BITS) {
                throw std::runtime_error("nbit must be in [32, 256]");
            }
            num_bits = static_cast<uint32_t>(parsed);
            saw_num_bits = true;
            continue;
        }
        throw std::runtime_error("Unrecognized argument: " + arg);
    }

    if (expected_crs_hash.empty()) {
        if (const char* env_hash = std::getenv("AUTOHOG_CRS_SHA256")) {
            expected_crs_hash = env_hash;
        }
    }

    std::cout << "Configuration:\n";
    std::cout << "  - Bits to encrypt: " << num_bits << "\n";
    std::cout << "  - Payload bits:     " << num_bits << "\n";
    std::cout << "  - Proved c2 slots:  pending CRS load\n";
    std::cout << "  - RLWE dimension:   " << P::n << "\n";
    std::cout << "  - Parameter set:    lvl2param\n";
    std::cout << "  - CRS file:         " << crs_file << "\n\n";

    std::chrono::duration<double, std::milli> t_keygen_total(0.0);
    std::chrono::duration<double, std::milli> t_crs_gen(0.0);
    std::chrono::duration<double, std::milli> t_encrypt_zk(0.0);
    std::chrono::duration<double, std::milli> t_proof_verify(0.0);
    std::chrono::duration<double, std::milli> t_pbs_total(0.0);
    std::chrono::duration<double, std::milli> t_pack(0.0);

    std::cout << "=== Step 1: Key Generation ===\n";
    SecretKey sk;
    EvalKey ek;

    {
        auto start = std::chrono::high_resolution_clock::now();
        ek.emplacebkfft<lvl02param>(sk);
        auto end = std::chrono::high_resolution_clock::now();
        t_keygen_total += std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start);
        std::cout << "BK Generation " << (end - start).count() / 1e6 << " ms\n";
    }
    {
        auto start = std::chrono::high_resolution_clock::now();
        ek.emplaceiksk<lvl20param>(sk);
        auto end = std::chrono::high_resolution_clock::now();
        t_keygen_total += std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start);
        std::cout << "IKS Generation " << (end - start).count() / 1e6 << " ms\n";
    }

    std::cout << "Keys generated successfully\n\n";

    std::cout << "=== Step 2: ZK CRS Loading ===\n";

    ZKProof::PublicParams<P> zk_params;
    constexpr uint64_t ZK_B = 1ULL << 16;
    constexpr uint64_t ZK_BR = static_cast<uint64_t>(P::n / 2 + 1);
    constexpr uint64_t ZK_T = 32;
    constexpr uint64_t ZK_MSBS = 0;
    constexpr size_t ZK_EFFECTIVE_T = ZK_T >> ZK_MSBS;
    const size_t ZK_STEP_T = ilog2_u64_local(ZK_EFFECTIVE_T);
    const size_t ZK_STEP_B = 1 + ilog2_u64_local(ZK_B);
    const size_t ZK_STEP_BR = 1 + ilog2_u64_local(ZK_BR);
    static_assert(ZK_EFFECTIVE_T == 32);

    {
        auto start = std::chrono::high_resolution_clock::now();
        if (!std::filesystem::exists(crs_file)) {
            std::cerr
                << "ERROR: Missing CRS file: " << crs_file << "\n"
                << "Generate it first with Coprocessor-operator/gen_crs.\n";
            return -1;
	    }

	    if (expected_crs_hash.empty()) {
            std::cerr
                << "ERROR: Missing trusted CRS SHA-256 digest.\n"
                << "Pass --crs-sha256 <64-hex> or set AUTOHOG_CRS_SHA256.\n"
                << "The runtime no longer trusts a sibling .sha256 file by itself.\n";
            return -1;
        }

        std::cout << "Found existing CRS file: " << crs_file << "\n";
        std::cout << "Using pinned CRS SHA-256: " << expected_crs_hash << "\n";
        try {
            auto loaded_params = ZKProof::load_public_params_by_expected_hash<P>(
                crs_file,
                expected_crs_hash
            );
            if (loaded_params.k < num_bits) {
                std::cerr
                    << "ERROR: CRS statement width " << loaded_params.k
                    << " is smaller than requested payload width " << num_bits << ".\n"
                    << "Re-generate the CRS with a larger k.\n";
                return -1;
            }
            if (loaded_params.d != P::n ||
                loaded_params.n != P::n + loaded_params.k * ZK_STEP_T + (P::n + loaded_params.k) * (ZK_STEP_B + ZK_STEP_BR) + 1 ||
                loaded_params.b != ZK_B ||
                loaded_params.b_r != ZK_BR ||
                loaded_params.q != 0 ||
                loaded_params.t != ZK_T ||
                loaded_params.msbs_zero_padding_bit_count != ZK_MSBS) {
                std::cerr
                    << "ERROR: CRS file is incompatible with the expected lvl2 parameters.\n"
                    << "Re-generate it with Coprocessor-operator/gen_crs.\n";
                return -1;
            }
            zk_params = std::move(loaded_params);
            std::cout << "CRS loaded from file after SHA-256 check\n";
        } catch (const std::exception& e) {
            std::cerr
                << "ERROR: Failed to load CRS file (" << e.what() << ").\n"
                << "Re-generate it or update the pinned SHA-256 digest.\n";
            return -1;
        }

        auto end = std::chrono::high_resolution_clock::now();
        t_crs_gen = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start);
        std::cout << "CRS ready in " << t_crs_gen.count() << " ms\n";
    }

    // A valid subgroup point is not sufficient: changing the order of two
    // powers must also be rejected by the CRS power-chain validation.
    {
        auto* g_list = static_cast<blst_p1_affine*>(zk_params.g_list_ptr);
        std::swap(g_list[0], g_list[1]);
        const bool accepted_tampered_chain = ZKProof::validate_public_params(zk_params);
        std::swap(g_list[0], g_list[1]);
        const bool restored_chain_is_valid = ZKProof::validate_public_params(zk_params);
        if (accepted_tampered_chain || !restored_chain_is_valid) {
            std::cerr << "ERROR: CRS power-chain tamper test failed.\n";
            return -1;
        }
        std::cout << "CRS power-chain tamper rejection test passed\n";
    }

    std::cout << "  - Witness dimension: " << zk_params.n << "\n";
    std::cout << "  - RLWE dimension:    " << zk_params.d << "\n";
    std::cout << "  - Proved c2 slots:   " << zk_params.k
              << " (payload uses " << num_bits << ", max app payload 256)\n";
    std::cout << "  - CRS max slots:     " << zk_params.k << "\n";
    printFileSize(crs_file);
    std::cout << "\n";

    std::cout << "=== Step 3: Plaintext Preparation ===\n";

    std::vector<uint8_t> plaintext(num_bits);
    for (uint32_t i = 0; i < num_bits; ++i) {
        plaintext[i] = static_cast<uint8_t>(coin(engine));
    }

    std::cout << "Plaintext (first 32 bits): ";
    for (uint32_t i = 0; i < std::min(32u, num_bits); ++i) {
        std::cout << int(plaintext[i]);
    }
    std::cout << (num_bits > 32 ? "...\n" : "\n");
    std::cout << "\n";

    std::cout << "=== Step 4: Client — Encrypt + Prove + Serialize ===\n";

    const std::string sealed_file = "sealed_ciphertext.bin";
    std::vector<uint8_t> metadata = {'Z', 'K', 'P', 'K', 'E'};

    auto rlwe_pk = PKE::generate_rlwe_pk_pair_seeded<P>(sk, engine);

    ZKPKE::EncryptionResult<P> enc_result;
    {
        auto start = std::chrono::high_resolution_clock::now();
        enc_result = ZKPKE::encrypt_bits_with_zk<P>(plaintext, rlwe_pk, engine, zk_params);
        auto end = std::chrono::high_resolution_clock::now();
        t_encrypt_zk = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start);
        std::cout << "Encrypt + Prove " << t_encrypt_zk.count() << " ms\n";
    }

    saveSealedCiphertext<P>(enc_result.ciphertext, enc_result.proof, metadata, sealed_file);
    printFileSize(sealed_file);
    std::cout << "\n";

    std::cout << "=== Step 5: Server — Load + Verify ZK ===\n";

    TRLWE<P> server_ct;
    ZKProof::Proof server_proof;
    std::vector<uint8_t> server_metadata;

    if (!loadSealedCiphertext<P>(server_ct, server_proof, server_metadata, sealed_file)) {
        std::cerr << "ERROR: Failed to load sealed ciphertext!\n";
        return -1;
    }
    std::cout << "Metadata: \"" << std::string(server_metadata.begin(), server_metadata.end()) << "\"\n";

    std::vector<TFHEpp::TLWE<P>> tlwe_vec;
    bool zk_valid = false;

    {
        auto start = std::chrono::high_resolution_clock::now();

        zk_valid = ZKPKE::verify_and_decompress<P>(
            server_proof, zk_params, rlwe_pk, server_ct, tlwe_vec, num_bits, server_metadata
        );

        auto end = std::chrono::high_resolution_clock::now();
        t_proof_verify = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start);
    }

    if (zk_valid) {
        std::cout << "✅ ZK Proof verification PASSED!\n";
    } else {
        std::cerr << "❌ ZK verification FAILED! Aborting.\n";
        return -1;
    }
    std::cout << "ZK Verification " << t_proof_verify.count() << " ms\n\n";

    std::cout << "=== Step 6: Security / Robustness Tests ===\n";

    auto expect_reject = [&](const std::string& name,
                             const ZKProof::Proof& proof,
                             const PKE::SeededRlwePkPair<P>& pk,
                             const TFHEpp::TRLWE<P>& ct,
                             const std::vector<uint8_t>& meta,
                             size_t num_bits) {
        std::vector<TFHEpp::TLWE<P>> dummy;
        bool ok = ZKPKE::verify_and_decompress<P>(
            proof, zk_params, pk, ct, dummy, num_bits, meta
        );

        if (!ok) {
            std::cout << "  ✅ PASS: " << name << "\n";
        } else {
            std::cerr << "  ❌ FAIL: " << name << " NOT rejected!\n";
            std::exit(-1);
        }
    };

    auto expect_accept = [&](const std::string& name,
                             const ZKProof::Proof& proof,
                             const PKE::SeededRlwePkPair<P>& pk,
                             const TFHEpp::TRLWE<P>& ct,
                             const std::vector<uint8_t>& meta,
                             size_t num_bits) {
        std::vector<TFHEpp::TLWE<P>> dummy;
        bool ok = ZKPKE::verify_and_decompress<P>(
            proof, zk_params, pk, ct, dummy, num_bits, meta
        );

        if (ok) {
            std::cout << "  ✅ PASS: " << name << "\n";
        } else {
            std::cerr << "  ❌ FAIL: " << name << " unexpectedly rejected!\n";
            std::exit(-1);
        }
    };

    {
        TRLWE<P> ct_tampered = server_ct;
        ct_tampered[0][0] ^= 1;
        std::cout << "Test 1: Flip c1[0] LSB\n";
        expect_reject("c1[0] tamper detected", server_proof, rlwe_pk, ct_tampered, server_metadata, num_bits);
    }

    {
        TRLWE<P> ct_tampered = server_ct;
        ct_tampered[P::k][0] ^= 1;
        std::cout << "Test 2: Flip c2[0] LSB\n";
        expect_reject("c2[0] tamper detected", server_proof, rlwe_pk, ct_tampered, server_metadata, num_bits);
    }

    {
        TRLWE<P> ct_tampered = server_ct;
        ct_tampered[P::k][1] ^= 1;
        std::cout << "Test 3: Flip c2[1] LSB (security critical)\n";
        expect_reject("c2[1] tamper detected", server_proof, rlwe_pk, ct_tampered, server_metadata, num_bits);
    }

    {
        TRLWE<P> ct_tampered = server_ct;
        ct_tampered[P::k][num_bits - 1] ^= 1;
        std::cout << "Test 4: Flip c2[" << num_bits - 1 << "] LSB\n";
        expect_reject("c2[last] tamper detected", server_proof, rlwe_pk, ct_tampered, server_metadata, num_bits);
    }

    {
        TRLWE<P> ct_tampered = server_ct;
        ct_tampered[P::k][P::n - 1] ^= 1;
        std::cout << "Test 4b: Flip c2[" << (P::n - 1) << "] LSB outside ZK statement\n";
        expect_accept("unused c2 tail does not affect proved output", server_proof, rlwe_pk, ct_tampered, server_metadata, num_bits);
    }

    {
        std::cout << "Test 4c: Request one output beyond the ZK statement\n";
        expect_reject("unproved output range rejected", server_proof, rlwe_pk, server_ct,
                      server_metadata, zk_params.k + 1);
    }

    {
        std::cout << "Test 5: Tampered metadata\n";
        std::vector<uint8_t> bad_meta = {'H', 'A', 'C', 'K'};
        expect_reject("metadata tamper detected", server_proof, rlwe_pk, server_ct, bad_meta, num_bits);
    }

    {
        std::cout << "Test 6: Replay attack (same pk, new ciphertext, old proof)\n";

        const auto& pk = rlwe_pk;

        const auto v2  = PKE::sample_gaussian_poly<P>(engine, P::α);
        const auto e0_2 = PKE::sample_gaussian_poly<P>(engine, P::α);
        std::array<TFHEpp::Polynomial<P>, P::k> e1_2{};
        for (size_t row = 0; row < P::k; ++row) {
            e1_2[row] = PKE::sample_gaussian_poly<P>(engine, P::α);
        }

        TFHEpp::Polynomial<P> m2{};
        for (size_t j = 0; j < P::n; ++j) {
            m2[j] = (j < plaintext.size() && plaintext[j]) ? P::μ : (typename P::T(0) - P::μ);
        }

        TFHEpp::TRLWE<P> ct2;
        TFHEpp::Polynomial<P> ai, tmp, ca, cb;

        for (size_t row = 0; row < P::k; ++row) {
            PKE::uniform_poly_from_seed<P>(ai, pk.seed, row);
            TFHEpp::PolyMul<P>(tmp, ai, v2);
            for (size_t j = 0; j < P::n; ++j) {
                ca[j] = tmp[j] + e1_2[row][j];
            }
            PKE::assign_poly_to_trlwe_comp<P>(ct2, row, ca);
        }

        TFHEpp::PolyMul<P>(cb, pk.b, v2);
        for (size_t j = 0; j < P::n; ++j) {
            cb[j] += e0_2[j] + m2[j];
        }
        PKE::assign_poly_to_trlwe_comp<P>(ct2, P::k, cb);

        expect_reject("replay attack correctly rejected", server_proof, pk, ct2, server_metadata, num_bits);
    }

    {
        std::cout << "Test 7: Proof serialize/deserialize round-trip\n";

        auto proof_bytes = server_proof.serialize();
        auto proof2 = ZKProof::Proof::deserialize(proof_bytes);

        expect_accept("round-trip proof accepted", proof2, rlwe_pk, server_ct, server_metadata, num_bits);

        auto proof_bytes_2 = proof2.serialize();
        if (proof_bytes != proof_bytes_2) {
            std::cerr << "  ❌ FAIL: proof round-trip bytes mismatch!\n";
            return -1;
        } else {
            std::cout << "  ✅ PASS: proof round-trip bytes identical\n";
        }
    }

    {
        std::cout << "Test 8: Tampered proof bytes (c_hat)\n";

        auto proof_bytes = server_proof.serialize();
        proof_bytes[autohog::zk::config::g2_serialized_size - 1] ^= 0x01;

        bool rejected = false;
        try {
            auto bad_proof = ZKProof::Proof::deserialize(proof_bytes);
            std::vector<TFHEpp::TLWE<P>> dummy;
            bool ok = ZKPKE::verify_and_decompress<P>(
                bad_proof, zk_params, rlwe_pk, server_ct, dummy, num_bits, server_metadata
            );
            rejected = !ok;
        } catch (const std::exception&) {
            rejected = true;
        }

        if (rejected) {
            std::cout << "  ✅ PASS: tampered c_hat proof rejected\n";
        } else {
            std::cerr << "  ❌ FAIL: tampered c_hat proof accepted!\n";
            return -1;
        }
    }

    {
        std::cout << "Test 9: Tampered proof bytes (c_y)\n";

        auto proof_bytes = server_proof.serialize();
        proof_bytes[autohog::zk::config::g2_serialized_size] ^= 0x01;

        bool rejected = false;
        try {
            auto bad_proof = ZKProof::Proof::deserialize(proof_bytes);
            std::vector<TFHEpp::TLWE<P>> dummy;
            bool ok = ZKPKE::verify_and_decompress<P>(
                bad_proof, zk_params, rlwe_pk, server_ct, dummy, num_bits, server_metadata
            );
            rejected = !ok;
        } catch (const std::exception&) {
            rejected = true;
        }

        if (rejected) {
            std::cout << "  ✅ PASS: tampered c_y proof rejected\n";
        } else {
            std::cerr << "  ❌ FAIL: tampered c_y proof accepted!\n";
            return -1;
        }
    }

    {
        std::cout << "Test 10: Truncated proof bytes\n";

        auto proof_bytes = server_proof.serialize();
        if (proof_bytes.size() < 10) {
            std::cerr << "  ❌ FAIL: proof too short unexpectedly\n";
            return -1;
        }
        proof_bytes.resize(proof_bytes.size() - 10);

        bool rejected = false;
        try {
            auto bad_proof = ZKProof::Proof::deserialize(proof_bytes);
            (void)bad_proof;
        } catch (const std::exception&) {
            rejected = true;
        }

        if (rejected) {
            std::cout << "  ✅ PASS: truncated proof rejected\n";
        } else {
            std::cerr << "  ❌ FAIL: truncated proof accepted!\n";
            return -1;
        }
    }

    {
        std::cout << "Test 11: Default all-infinity proof\n";
        ZKProof::Proof forged_proof;
        expect_reject("all-infinity proof rejected", forged_proof, rlwe_pk, server_ct, server_metadata, num_bits);
    }

    {
        std::cout << "Test 12: Generator-point forged proof\n";
        ZKProof::Proof forged_proof;
        *static_cast<blst_p2_affine*>(forged_proof.c_hat_ptr) = *blst_p2_affine_generator();
        *static_cast<blst_p1_affine*>(forged_proof.c_y_ptr) = *blst_p1_affine_generator();
        *static_cast<blst_p1_affine*>(forged_proof.pi_ptr) = *blst_p1_affine_generator();
        expect_reject("generator-point forged proof rejected", forged_proof, rlwe_pk, server_ct, server_metadata, num_bits);
    }

    {
        std::cout << "Test 13: Wrong public key binding\n";
        auto wrong_pk = PKE::generate_rlwe_pk_pair_seeded<P>(sk, engine);
        expect_reject("wrong public key rejected", server_proof, wrong_pk, server_ct, server_metadata, num_bits);
    }

    {
        std::cout << "Test 14: Non-binary message witness\n";

        auto malicious_ct = server_ct;
        auto malicious_witness = enc_result.witness;
        const size_t stmt_idx = num_bits - 1;
        const size_t witness_idx = enc_result.witness.m.size() - 1 - stmt_idx;
        const int64_t old_m = malicious_witness.m[witness_idx];
        malicious_witness.m[witness_idx] = 2;

        const __uint128_t decoded_q = (__uint128_t(1) << 64);
        const typename P::T delta_torus = static_cast<typename P::T>(decoded_q / ZK_T);
        const typename P::T adjustment = static_cast<typename P::T>(
            (__uint128_t)(malicious_witness.m[witness_idx] - old_m) * (__uint128_t)delta_torus
        );
        malicious_ct[P::k][stmt_idx] += adjustment;

        auto malicious_public_commit = ZKProof::make_public_commit<P>(rlwe_pk, malicious_ct, zk_params.k);
        std::array<uint8_t, ZKProof::SEED_LEN> malicious_seed{};
        for (auto& byte : malicious_seed) byte = static_cast<uint8_t>(engine() & 0xFF);
        try {
            auto malicious_proof = ZKProof::prove<P>(
                zk_params,
                malicious_public_commit,
                malicious_witness,
                server_metadata,
                ZKProof::ComputeLoad::Verify,
                malicious_seed
            );

            expect_reject("non-binary plaintext proof rejected", malicious_proof, rlwe_pk, malicious_ct, server_metadata, num_bits);
        } catch (const std::exception& e) {
            std::cout << "  ✅ PASS: non-binary plaintext witness rejected during proving (" << e.what() << ")\n";
        }
    }

    std::cout << "\nAll security / robustness tests PASSED!\n\n";

    std::cout << "=== Step 7: Decompress & Check Decryption ===\n";

    bool decrypt_ok = true;
    uint32_t errors = 0;

    for (uint32_t i = 0; i < num_bits; ++i) {
        uint8_t decrypted = bootsSymDecrypt<P>({tlwe_vec[i]}, sk).at(0);
        if (decrypted != plaintext[i]) {
            if (errors < 5) {
                std::cout << "  [Mismatch] bit " << i
                          << ": expected " << int(plaintext[i])
                          << ", got " << int(decrypted) << "\n";
            }
            decrypt_ok = false;
            errors++;
        }
    }

    if (decrypt_ok) {
        std::cout << "All " << num_bits << " bits decrypted correctly\n";
    } else {
        std::cout << "ERROR: " << errors << "/" << num_bits << " bits incorrect\n";
        return -1;
    }
    std::cout << "\n";

    std::cout << "=== Step 8: Bootstrap (Homomorphic Refresh) ===\n";

    std::vector<TLWE<P>> booted_vec(num_bits);

    {
        auto start = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for
        for (uint32_t i = 0; i < num_bits; ++i) {
            GateBootstrapping<lvl20param, lvl02param, P::μ>(booted_vec[i], tlwe_vec[i], ek);
        }

        auto end = std::chrono::high_resolution_clock::now();
        t_pbs_total = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start);

        std::cout << "Bootstrap " << num_bits << " bits: " << t_pbs_total.count() << " ms\n";
        std::cout << "Avg per bit: " << t_pbs_total.count() / num_bits << " ms\n";
    }

    bool pbs_ok = true;
    errors = 0;
    for (uint32_t i = 0; i < num_bits; ++i) {
        uint8_t after_pbs = bootsSymDecrypt<P>({booted_vec[i]}, sk).at(0);
        if (after_pbs != plaintext[i]) {
            if (errors < 5) {
                std::cout << "  [PBS Error] bit " << i
                          << ": expected " << int(plaintext[i])
                          << ", got " << int(after_pbs) << "\n";
            }
            pbs_ok = false;
            errors++;
        }
    }

    if (pbs_ok) {
        std::cout << "Bootstrap successful, all bits preserved\n";
    } else {
        std::cout << "ERROR: PBS corrupted " << errors << "/" << num_bits << " bits\n";
        return -1;
    }
    std::cout << "\n";

    std::cout << "========================================\n";
    std::cout << "           Performance Summary\n";
    std::cout << "========================================\n\n";

    std::cout << "Key Generation:       " << t_keygen_total.count() << " ms\n";
    std::cout << "CRS Loading / Check:  " << t_crs_gen.count() << " ms\n";
    std::cout << "Encryption + Prove:   " << t_encrypt_zk.count() << " ms\n";
    std::cout << "ZK Verification:      " << t_proof_verify.count() << " ms\n";
    std::cout << "Bootstrap (" << num_bits << " bits): " << t_pbs_total.count() << " ms\n";
    std::cout << "  Avg per bit:        " << t_pbs_total.count() / num_bits << " ms\n";
   
    std::cout << "File Sizes:\n";
    printFileSize(crs_file);
    printFileSize(sealed_file);
    std::cout << "\n";

    std::cout << "========================================\n";
    std::cout << "  All tests PASSED!\n";
    std::cout << "========================================\n";

    return 0;

} catch (const std::exception& e) {
    std::cerr << "\nFATAL ERROR: " << e.what() << "\n";
    return -1;
}
