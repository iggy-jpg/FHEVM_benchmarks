#include <iostream>
#include <vector>
#include <cassert>
#include <cstdint>
#include <random>
#include <algorithm>


#include <tfhe++.hpp>


#include "HomPKE.h"
#include "SecureRandom.h"


#define ASSERT_TRUE(condition, message) \
    do { \
        if (!(condition)) { \
            std::cerr << "[\033[31mFAIL\033[0m] " << message << std::endl; \
            std::exit(1); \
        } \
    } while (0)

#define LOG_PASS(test_name) \
    std::cout << "[\033[32mPASS\033[0m] " << test_name << std::endl;

// -----------------------------------------------------------------------------
// Regression test: 64-bit Torus scaling must not narrow 2^64 to uint64_t.
// -----------------------------------------------------------------------------
void test_torus_noise_scaling() {
    std::cout << "Testing 32/64-bit Torus noise scaling..." << std::endl;

    const std::uint64_t expected64 = std::uint64_t{1} << 24;
    const auto positive64 = PKE::torus_from_normalized_noise<std::uint64_t>(
        std::ldexp(1.0L, -40));
    const auto negative64 = PKE::torus_from_normalized_noise<std::uint64_t>(
        -std::ldexp(1.0L, -40));
    ASSERT_TRUE(positive64 == expected64,
                "64-bit Torus positive noise scaling is incorrect");
    ASSERT_TRUE(negative64 == std::uint64_t{0} - expected64,
                "64-bit Torus negative noise scaling is incorrect");

    const std::uint32_t expected32 = std::uint32_t{1} << 12;
    const auto positive32 = PKE::torus_from_normalized_noise<std::uint32_t>(
        std::ldexp(1.0L, -20));
    ASSERT_TRUE(positive32 == expected32,
                "32-bit Torus positive noise scaling is incorrect");

    LOG_PASS("32/64-bit Torus noise scaling");
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void test_lwe_pke() {
    std::cout << "Testing Standard LWE PKE..." << std::endl;
    
    
    TFHEpp::SecretKey sk;
    using P = TFHEpp::lvl0param; 
    
    
    auto pk = PKE::generate_public_key<P>(sk);
    
    autohog::SecureRng rng;

    
    {
        uint8_t msg = 0;
        auto ct = PKE::encrypt_public_key_bit<P>(pk, msg, rng);
        bool decrypted = TFHEpp::tlweSymDecrypt<P>(ct, sk.key.get<P>());
        ASSERT_TRUE(decrypted == (msg == 1), "LWE PKE failed for message 0");
    }

    
    {
        uint8_t msg = 1;
        auto ct = PKE::encrypt_public_key_bit<P>(pk, msg, rng);
        bool decrypted = TFHEpp::tlweSymDecrypt<P>(ct, sk.key.get<P>());
        ASSERT_TRUE(decrypted == (msg == 1), "LWE PKE failed for message 1");
    }

    LOG_PASS("Standard LWE PKE");
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void test_compressed_lwe_pke() {
    std::cout << "Testing Compressed (Seeded) LWE PKE..." << std::endl;

    TFHEpp::SecretKey sk;
    using P = TFHEpp::lvl0param;

    
    auto cpk = PKE::generate_compressed_public_key<P>(sk);

    autohog::SecureRng rng;

    int trials = 50;
    int correct = 0;

    for (int i = 0; i < trials; ++i) {
        uint8_t msg = i % 2; 
        auto ct = PKE::encrypt_with_compressed_public_key_bit<P>(cpk, msg, rng);
        bool decrypted = TFHEpp::tlweSymDecrypt<P>(ct, sk.key.get<P>());
        if (decrypted == (msg == 1)) correct++;
    }

    ASSERT_TRUE(correct == trials, "Compressed LWE PKE decryption failed in random trials");
    LOG_PASS("Compressed LWE PKE");
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void test_rlwe_pke() {
    std::cout << "Testing Standard RLWE PKE (Encrypting Poly)..." << std::endl;

    TFHEpp::SecretKey sk;
    using P = TFHEpp::lvl1param; 

    autohog::SecureRng rng;

    
    auto pk = PKE::generate_rlwe_pk_pair<P>(sk, rng);

    
    std::array<bool, P::n> msg_bits;
    for (auto& b : msg_bits) b = rng() % 2;

    
    auto ct = PKE::encrypt_bits_trlwe_with_rlwe_pair<P>(pk, msg_bits, rng);

    
    
    std::array<bool, P::n> decrypted = TFHEpp::trlweSymDecrypt<P>(ct, sk.key.get<P>());

    
    bool all_match = true;
    for (size_t i = 0; i < P::n; ++i) {
        if (decrypted[i] != msg_bits[i]) {
            all_match = false;
            break;
        }
    }

    ASSERT_TRUE(all_match, "RLWE PKE failed: decrypted bits do not match original");
    LOG_PASS("Standard RLWE PKE");
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void test_compressed_rlwe_pke() {
    std::cout << "Testing Compressed (Seeded) RLWE PKE..." << std::endl;

    TFHEpp::SecretKey sk;
    using P = TFHEpp::lvl1param;

    autohog::SecureRng rng;

    
    auto spk = PKE::generate_rlwe_pk_pair_seeded<P>(sk, rng);

    
    std::array<bool, P::n> msg_bits;
    for (auto& b : msg_bits) b = rng() % 2;

    
    auto ct = PKE::encrypt_bits_trlwe_with_seeded_rlwe_pair<P>(spk, msg_bits, rng);

    
    std::array<bool, P::n> decrypted = TFHEpp::trlweSymDecrypt<P>(ct, sk.key.get<P>());

    
    bool all_match = true;
    for (size_t i = 0; i < P::n; ++i) {
        if (decrypted[i] != msg_bits[i]) {
            all_match = false;
            break;
        }
    }

    ASSERT_TRUE(all_match, "Compressed RLWE PKE failed");
    LOG_PASS("Compressed RLWE PKE");
}

int main() {
    try {
        std::cout << "Starting PKE Tests..." << std::endl;
        std::cout << "=====================" << std::endl;

        test_torus_noise_scaling();
        test_lwe_pke();
        test_compressed_lwe_pke();
        test_rlwe_pke();
        test_compressed_rlwe_pke();

        std::cout << "=====================" << std::endl;
        std::cout << "All Tests Passed Successfully!" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
