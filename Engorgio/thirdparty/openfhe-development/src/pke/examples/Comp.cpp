//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*
  Example of evaluating arbitrary smooth functions with the Chebyshev approximation using CKKS.
 */

#include "openfhe.h"
#include <chrono>

using namespace lbcrypto;

void EvalLogisticExample();
Ciphertext<DCRTPoly> MultByInteger(const ConstCiphertext<DCRTPoly>& ciphertext, const int64_t constant) {
    auto result = ciphertext->Clone();
    for (auto& c : result->GetElements())
        c *= static_cast<typename DCRTPoly::Integer>(constant);
    return result;
}
void EvalFunctionExample();
template <typename T>
inline uint64_t CeilLog2(T x) {
    return static_cast<uint64_t>(std::ceil(std::log2(x)));
}
void print_moduli_chain(const DCRTPoly& poly) {
    int num_primes       = poly.GetNumOfElements();
    double total_bit_len = 0.0;
    for (int i = 0; i < num_primes; i++) {
        auto qi = poly.GetParams()->GetParams()[i]->GetModulus();
        std::cout << "q_" << i << ": " << qi.ConvertToInt() << ",  log q_" << i << ": "
                  << log(qi.ConvertToInt()) / log(2) << std::endl;
        total_bit_len += log(qi.ConvertToInt()) / log(2);
    }
    std::cout << "Total bit length: " << total_bit_len << std::endl;
}
std::vector<double> coeff1 = {1.5, -0.5};  // 1.5 x  - 0.5 x ^ 3
std::vector<double> coeff3 = {2.1875, -2.1875, 1.3125, -0.3125};
std::vector<double> coeff5 = {2.4609375 / 2, -3.28125 / 2, 2.953125 / 2, -1.40625 / 2, 0.2734375 / 2};
std::vector<double> g      = {double(4589 / 1024), double(-16577 / 1024), double(25614 / 1024), double(-12860 / 1024)};
// -0.2095 x ^ 15 + 1.692 x ^ 13 + -5.999 x ^ 11 + 12.22 x ^ 9 + -15.71 x ^ 7 + 13.2 x ^ 5 + -7.332 x ^ 3 + 3.142 x
void EvalPower(std::vector<double> coefficients, std::vector<Ciphertext<lbcrypto::DCRTPoly>>& power_basis,
               Ciphertext<lbcrypto::DCRTPoly>& result) {
    auto cc = power_basis[0]->GetCryptoContext();
    if (coefficients.size() == 1) {
        if (std::fabs(std::round(coefficients[0] * pow(2, 50))) > 1.) {
            // double qT_d = parms.coeff_modulus()[power_basis[0].coeff_modulus_size() - 1].value();
            // double cipher_scale = target_scale * qT_d;
            // std::cout << "coefficients.size() == 1,cipher_scale: " << std::log2(cipher_scale / power_basis[0].scale()) << std::endl;
            result = power_basis[0];
            // std::cout << "coefficients.size() == 1" << std::endl;
            result = cc->EvalMult(result, coefficients[0]);
            return;
        }
        else {
            return;
        }
    }
    std::vector<double> quotient, remainder;
    uint64_t degree = 2 * coefficients.size() - 1;
    uint64_t m      = CeilLog2(degree + 1);
    remainder.resize((1 << (m - 1)) / 2);
    quotient.resize(coefficients.size() - remainder.size());

    for (size_t i = 0; i < remainder.size(); i++) {
        remainder[i] = coefficients[i];
    }
    for (size_t i = 0; i < quotient.size(); i++) {
        quotient[i] = coefficients[i + remainder.size()];
    }
    Ciphertext<lbcrypto::DCRTPoly> cipher_quotient, cipher_remainder;
    EvalPower(quotient, power_basis, cipher_quotient);
    EvalPower(remainder, power_basis, cipher_remainder);
    // std::cout << "EvalPower Mult" << std::endl;
    result = cc->EvalMult(cipher_quotient, power_basis[m - 1]);
    // std::cout << "EvalAdd" << std::endl;
    result = cc->EvalAdd(result, cipher_remainder);
}

void poly_evaluate_power(Ciphertext<lbcrypto::DCRTPoly>& result, Ciphertext<lbcrypto::DCRTPoly>& x,
                         std::vector<double>& coefficients) {
    auto cc         = x->GetCryptoContext();
    uint64_t degree = coefficients.size() * 2 - 1;
    uint64_t m      = CeilLog2(degree + 1);

    std::vector<Ciphertext<lbcrypto::DCRTPoly>> power_basis(m);
    power_basis[0] = x;
    for (size_t i = 1; i < m; i++) {
        power_basis[i] = cc->EvalMult(power_basis[i - 1], power_basis[i - 1]);
    }

    EvalPower(coefficients, power_basis, result);
}

void Homround(Ciphertext<lbcrypto::DCRTPoly>& cipher) {
    auto cc = cipher->GetCryptoContext();
    cipher  = MultByInteger(cipher, 2.0);
    cipher  = cc->EvalAdd(cipher, -1.0);
    poly_evaluate_power(cipher, cipher, coeff1);
    // poly_evaluate_power(cipher, scale, cipher, coeff3, context, encoder, evaluator, relin_keys, decryptor);
    poly_evaluate_power(cipher, cipher, coeff5);
    cipher = cc->EvalAdd(cipher, 0.5);
}

// Ciphertext<lbcrypto::DCRTPoly> approxcomp(Ciphertext<lbcrypto::DCRTPoly>& cipher) {
//     auto cc = cipher->GetCryptoContext();
//     std::cout << "poly_evaluate_power coeff3" << std::endl;
//     for (int i = 0; i < 4; i++)
//         poly_evaluate_power(cipher, cipher, coeff3);
//     // poly_evaluate_power(cipher, scale, cipher, coeff3, context, encoder, evaluator, relin_keys, decryptor);
//     // std::cout << "poly_evaluate_power g" << std::endl;
//     for (int i = 0; i < 4; i++) {
//         std::cout << "poly_evaluate_power g" << std::endl;
//         poly_evaluate_power(cipher, cipher, g);
//     }

//     std::cout << "EvalAdd 1" << std::endl;
//     auto cipher_add = cc->EvalAdd(cipher, 1);
//     std::cout << "EvalMult cipher*0.5" << std::endl;
//     auto res = cc->EvalMult(cipher_add, 0.5);
//     return res;
// }

void EvalSignExample(uint32_t plain_bits) {
    std::cout << "--------------------------------- EVAL LOGISTIC FUNCTION ---------------------------------"
              << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    // We set a smaller ring dimension to improve performance for this example.
    // In production environments, the security level should be set to
    // HEStd_128_classic, HEStd_192_classic, or HEStd_256_classic for 128-bit, 192-bit,
    // or 256-bit security, respectively.
    parameters.SetSecurityLevel(HEStd_128_classic);
    // parameters.SetRingDim(1 << 16);
#if NATIVEINT == 128
    usint scalingModSize = 78;
    usint firstModSize   = 89;
#else
    usint scalingModSize = 59;
    usint firstModSize   = 60;
#endif
    parameters.SetScalingModSize(scalingModSize);
    parameters.SetFirstModSize(firstModSize);

    // Choosing a higher degree yields better precision, but a longer runtime.
    uint32_t polyDegree = 59;

    // The multiplicative depth depends on the polynomial degree.
    // See the FUNCTION_EVALUATION.md file for a table mapping polynomial degrees to multiplicative depths.
    uint32_t multDepth = 27;

    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.
    cc->Enable(ADVANCEDSHE);
    double precision  = (1 << (plain_bits - 1) - 1);
    double lowerBound = -precision;
    double upperBound = precision;
    double bound      = 3;
    auto keyPair      = cc->KeyGen();
    // We need to generate mult keys to run Chebyshev approximations.
    cc->EvalMultKeyGen(keyPair.secretKey);
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message_1bit(-1, 1);
    std::uniform_int_distribution<int> message_part_1(-2, -1);
    std::uniform_int_distribution<int> message_part_2(1, 2);
    // std::uniform_int_distribution<int> message_part_16_1(-lowerBound, -1);
    std::uniform_int_distribution<int> message_part_16_2(1, upperBound);

    usint ringDim = cc->GetRingDimension();
    int length    = 128;
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl << std::endl;
    std::vector<std::complex<double>> input0_1(length), input_4bit(length), input_16bit(length);
    for (int i = 0; i < length; i++) {
        input0_1[i] = double(message_1bit(engine));
    }
    for (int i = 0; i < length / 2; i++) {
        input_4bit[2 * i]     = double(message_part_1(engine));
        input_4bit[2 * i + 1] = double(message_part_2(engine));
    }
    for (int i = 0; i < length / 2; i++) {
        input_16bit[2 * i]     = double(message_part_16_2(engine));
        input_16bit[2 * i + 1] = double(message_part_16_2(engine));
    }
    size_t encodedLength = input0_1.size();
    std::cout << "MakeCKKSPackedPlaintext " << std::endl << std::endl;
    Plaintext plaintext    = cc->MakeCKKSPackedPlaintext(input0_1);
    Plaintext plaintext_4  = cc->MakeCKKSPackedPlaintext(input_4bit);
    Plaintext plaintext_16 = cc->MakeCKKSPackedPlaintext(input_16bit);
    std::cout << "Encrypt " << std::endl << std::endl;
    auto ciphertext    = cc->Encrypt(keyPair.publicKey, plaintext);
    auto ciphertext_4  = cc->Encrypt(keyPair.publicKey, plaintext_4);
    auto ciphertext_16 = cc->Encrypt(keyPair.publicKey, plaintext_16);

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    std::cout << "EvalSign " << std::endl << std::endl;
    auto result = cc->EvalSign(ciphertext, bound, lowerBound, upperBound, polyDegree);

    // auto result_4  = cc->EvalSign(ciphertext_4, bound, lowerBound, upperBound, polyDegree);
    // auto result_16 = cc->EvalSign(ciphertext_16, bound, lowerBound, upperBound, polyDegree);
    // std::cout << "Homround " << std::endl << std::endl;
    // Homround(result_16);
    // Homround(result_4);
    // Homround(result);
    // end = std::chrono::system_clock::now();
    // //round
    // // result_4  = cc->EvalRound(result_4, bound, 0, 1, polyDegree);
    // // result_16 = cc->EvalRound(result_16, bound, -5, 5, polyDegree);
    // Plaintext plaintextDec, plaintextDec_4, plaintextDec_16;

    // // std::cout << "Decrypt plaintextDec" << std::endl;
    // cc->Decrypt(keyPair.secretKey, result, &plaintextDec);
    // plaintextDec->SetLength(encodedLength);

    // // std::cout << "Decrypt plaintextDec_4" << std::endl;
    // cc->Decrypt(keyPair.secretKey, result_4, &plaintextDec_4);
    // plaintextDec_4->SetLength(encodedLength);

    // // std::cout << "Decrypt plaintextDec_16" << std::endl;
    // cc->Decrypt(keyPair.secretKey, result_16, &plaintextDec_16);
    // plaintextDec_16->SetLength(encodedLength);

    // std::vector<std::complex<double>> expectedOutput_1;
    // std::vector<std::complex<double>> expectedOutput_4;
    // std::vector<std::complex<double>> expectedOutput_16;
    // std::vector<std::complex<double>> finalResult = plaintextDec->GetCKKSPackedValue();
    // double err                                    = 0;

    // for (int i = 0; i < (int)encodedLength; i++) {
    //     if (fabs(input0_1[i].real()) > 0)
    //         expectedOutput_1.push_back(static_cast<double>(input0_1[i].real() > 0));
    //     else
    //         expectedOutput_1.push_back(0.5);
    //     err += std::abs(expectedOutput_1[i] - finalResult[i]);
    // }
    // std::cout << "---------------------------------0bit ---------------------------------" << std::endl;
    // std::cout << "2^" << std::log2(err / encodedLength) << std::endl;
    // std::cout << "input\n\t" << input0_1 << std::endl;
    // std::cout << "Expected output\n\t" << expectedOutput_1 << std::endl;
    // std::cout << "Actual output\n\t" << finalResult << std::endl << std::endl;
    // finalResult = plaintextDec_4->GetCKKSPackedValue();
    // err         = 0.;
    // for (int i = 0; i < (int)encodedLength; i++) {
    //     expectedOutput_4.push_back(static_cast<double>(input_4bit[i].real() > 0));
    //     err += std::fabs(expectedOutput_4[i] - finalResult[i]);
    // }
    // std::cout << "---------------------------------3bit ---------------------------------" << std::endl;

    // std::cout << "2^" << std::log2(err / encodedLength) << std::endl;
    // // std::cout << "input\n\t" << input_16bit << std::endl;
    // // std::cout << "Expected output\n\t" << expectedOutput_1 << std::endl;
    // // std::cout << "Actual output\n\t" << finalResult << std::endl << std::endl;
    // finalResult = plaintextDec_16->GetCKKSPackedValue();
    // err         = 0.;
    // for (int i = 0; i < (int)encodedLength; i++) {
    //     expectedOutput_16.push_back(static_cast<double>(input_16bit[i].real() > 0));
    //     err += std::fabs(expectedOutput_16[i] - finalResult[i]);
    // }
    // std::cout << "---------------------------------11bit ---------------------------------" << std::endl;
    // std::cout << "2^" << std::log2(err / encodedLength) << std::endl;
    // std::cout << "time:" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 3 << "ms"
    //           << std::endl;
}

void EvalgreaterExample(uint32_t plain_bits) {
    std::cout << "\nplain_bits greater: " << plain_bits << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    // We set a smaller ring dimension to improve performance for this example.
    // In production environments, the security level should be set to
    // HEStd_128_classic, HEStd_192_classic, or HEStd_256_classic for 128-bit, 192-bit,
    // or 256-bit security, respectively.
    parameters.SetSecurityLevel(HEStd_128_classic);
    // parameters.SetRingDim(1 << 16);
#if NATIVEINT == 128
    usint scalingModSize = 78;
    usint firstModSize   = 89;
#else
    usint scalingModSize = 40;
    usint firstModSize   = 60;
#endif
    parameters.SetScalingModSize(scalingModSize);
    parameters.SetFirstModSize(firstModSize);

    // Choosing a higher degree yields better precision, but a longer runtime.
    uint32_t polyDegree = 35;

    // The multiplicative depth depends on the polynomial degree.
    // See the FUNCTION_EVALUATION.md file for a table mapping polynomial degrees to multiplicative depths.
    uint32_t multDepth = 20;

    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.
    cc->Enable(ADVANCEDSHE);
    double precision  = (1 << (plain_bits - 1) - 1);
    double lowerBound = -precision;
    double upperBound = precision;
    double bound      = 3;
    auto keyPair      = cc->KeyGen();
    // We need to generate mult keys to run Chebyshev approximations.
    cc->EvalMultKeyGen(keyPair.secretKey);
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    // std::uniform_int_distribution<int> message_part_16_1(-2048, -1);
    std::uniform_int_distribution<int> message(0, precision);
    usint ringDim = cc->GetRingDimension();
    int length    = ringDim / 2;

    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;
    std::vector<std::complex<double>> input_a(length), input_b(length);

    for (int i = 0; i < length / 2; i++) {
        input_a[i] = double(message(engine));
        input_b[i] = double(message(engine));
    }
    size_t encodedLength  = input_a.size();
    Plaintext plaintext_a = cc->MakeCKKSPackedPlaintext(input_a);
    Plaintext plaintext_b = cc->MakeCKKSPackedPlaintext(input_b);
    auto ciphertext_a     = cc->Encrypt(keyPair.publicKey, plaintext_a);
    auto ciphertext_b     = cc->Encrypt(keyPair.publicKey, plaintext_b);
    //comp(a>b),  a-b
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    std::cout << "sub " << std::endl;
    auto ciphertext_sub = cc->EvalSub(ciphertext_a, ciphertext_b);
    std::cout << "scaledown " << std::endl;
    auto ciphertext_sub_scaledown = cc->EvalMult(ciphertext_sub, double(1 / precision));
    std::cout << "approxcomp " << std::endl;
    // auto result_sub = approxcomp(ciphertext_sub_scaledown);
    auto result_sub = cc->EvalSign(ciphertext_sub_scaledown, bound, -1, 1, polyDegree);

    // auto result_sub               = cc->EvalSign(ciphertext_sub, bound, lowerBound, upperBound, polyDegree);
    // Homround(result_sub);

    end = std::chrono::system_clock::now();

    Plaintext plaintextDec, plaintextDec_comp_res, plaintextDec_comp_res_equal;

    cc->Decrypt(keyPair.secretKey, result_sub, &plaintextDec_comp_res);

    plaintextDec_comp_res->SetLength(encodedLength);
    std::vector<std::complex<double>> expectedOutput_compres;
    std::vector<std::complex<double>> finalResult = plaintextDec_comp_res->GetCKKSPackedValue();
    double err                                    = 0;
    int wrong                                     = 0;
    for (int i = 0; i < (int)encodedLength; i++) {
        if (input_a[i].real() != input_b[i].real())
            expectedOutput_compres.push_back(static_cast<double>(input_a[i].real() > input_b[i].real()));
        else
            expectedOutput_compres.push_back(0.5);
        if (fabs(expectedOutput_compres[i] - finalResult[i].real()) > 0.5) {
            wrong++;
            std::cout << "input_a:" << input_a[i] << ", input_b:" << input_b[i] << ", res :" << finalResult[i]
                      << std::endl;
        }

        err += std::fabs(expectedOutput_compres[i] - finalResult[i].real());
    }

    std::cout << "number of levels remaining before comp: " << multDepth - ciphertext_a->GetLevel() << std::endl;
    std::cout << "number of levels remaining after comp: " << multDepth - result_sub->GetLevel() << std::endl;
    std::cout << "---------------------------------compres---------------------------------" << std::endl;
    std::cout << "2^" << std::log2(err / encodedLength) << std::endl;
    std::cout << "wrong number:" << wrong << std::endl;
    std::cout << encodedLength
              << "slots total time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << "ms" << std::endl;
    std::cout << encodedLength << "slots amortize time: "
              << double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) /
                     double(encodedLength)
              << "ms" << std::endl
              << std::endl;
}

// 1, 0.5, 0
void Eval_modular_greater_test(uint32_t plain_bits, uint32_t block) {
    std::cout << plain_bits << " plain_bits " << block << " modular greater test " << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    // We set a smaller ring dimension to improve performance for this example.
    // In production environments, the security level should be set to
    // HEStd_128_classic, HEStd_192_classic, or HEStd_256_classic for 128-bit, 192-bit,
    // or 256-bit security, respectively.
    parameters.SetSecurityLevel(HEStd_128_classic);
    // parameters.SetRingDim(1 << 16);
#if NATIVEINT == 128
    usint scalingModSize = 78;
    usint firstModSize   = 89;
#else
    usint scalingModSize = 59;
    usint firstModSize   = 60;
#endif
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetScalingModSize(scalingModSize);
    parameters.SetFirstModSize(firstModSize);
    uint32_t polyDegree = 0;
    uint32_t multDepth  = 50;
    // Choosing a higher degree yields better precision, but a longer runtime.
    // if (plain_bits >= 8)
    //     polyDegree = 119;
    // else {
    //     polyDegree = 27;
    //     multDepth  = 20;
    // }
    // The multiplicative depth depends on the polynomial degree.
    // See the FUNCTION_EVALUATION.md file for a table mapping polynomial degrees to multiplicative depths.

    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.
    cc->Enable(ADVANCEDSHE);
    double precision  = (1 << (plain_bits - 1) - 1);
    double lowerBound = -precision;
    double upperBound = precision;
    double bound      = 3;
    auto keyPair      = cc->KeyGen();
    // We need to generate mult keys to run Chebyshev approximations.
    cc->EvalMultKeyGen(keyPair.secretKey);
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    // std::uniform_int_distribution<int> message_part_16_1(-2048, -1);
    std::uniform_int_distribution<int> message(0, precision);
    usint ringDim = cc->GetRingDimension();
    int length    = ringDim / 2;

    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;
    std::cout << "using polyDegree " << polyDegree << std::endl;
    std::vector<std::vector<std::complex<double>>> input_a_mod, input_b_mod;
    std::vector<double> input_a_plain(length, 0), input_b_plain(length, 0);
    //a0,a1...  b0,b1...

    for (int i = 0; i < block; i++) {
        std::vector<std::complex<double>> a(length), b(length);
        for (int j = 0; j < length; j++) {
            a[j] = double(message(engine));
            b[j] = double(message(engine));
        }
        input_a_mod.push_back(a);
        input_b_mod.push_back(b);
    }

    for (int i = 0; i < length; i++) {
        for (int j = 0; j < block; j++) {
            input_a_plain[i] += input_a_mod[j][i].real() * std::pow(2, j * plain_bits);
            input_b_plain[i] += input_b_mod[j][i].real() * std::pow(2, j * plain_bits);
        }
    }
    // for (int j = 0; j < 8; j++) {
    //     std::cout << "input_a:" << input_a_plain[j] << ", input_b:" << input_b_plain[j] << std::endl;
    // }
    double totalerr      = 0;
    double totaltime     = 0;
    size_t encodedLength = input_a_plain.size();

    for (int num_test = 0; num_test < 10; num_test++) {
        std::vector<Plaintext> plaintext_a, plaintext_b;
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_a, ciphertext_b;

        std::chrono::system_clock::time_point start, end;
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_sub(block), result_sub(block), result_equal_comp(block),
            ciphertext_comp_greater(block);
        for (int i = 0; i < block; i++) {
            Plaintext plain_a = cc->MakeCKKSPackedPlaintext(input_a_mod[i]);
            Plaintext plain_b = cc->MakeCKKSPackedPlaintext(input_b_mod[i]);
            ciphertext_a.push_back(cc->Encrypt(keyPair.publicKey, plain_a));
            ciphertext_b.push_back(cc->Encrypt(keyPair.publicKey, plain_b));
        }
        start = std::chrono::system_clock::now();
        for (int i = block - 1; i >= 0; i--) {
            //comp(a>b) 1, 0.5, 0

            ciphertext_sub[i] = cc->EvalSub(ciphertext_a[i], ciphertext_b[i]);
            result_sub[i]     = cc->EvalSign(ciphertext_sub[i], bound, lowerBound, upperBound, polyDegree);
            Homround(result_sub[i]);

            //comp(a==b) 1, 0

            auto ciphertext_sub_neg = cc->EvalNegate(result_sub[i]);
            auto result_sub_neg     = cc->EvalAdd(ciphertext_sub_neg, 1.0);

            auto result_equal    = cc->EvalMult(result_sub_neg, result_sub[i]);
            result_equal_comp[i] = MultByInteger(result_equal, 4.0);
            // ciphertext_comp_greater[i] = cc->EvalSub(result_sub[i], result_equal_comp[i]);
        }

        Ciphertext<lbcrypto::DCRTPoly> comp_res;
        for (int i = 0; i < block; i++) {
            // std::cout << "mux:" << i << std::endl;
            if (i == 0)
                comp_res = result_sub[i];
            else {
                // std::cout << "EvalSub" << std::endl;
                auto ai_sub_bi = cc->EvalSub(comp_res, result_sub[i]);
                // std::cout << "EvalMult" << std::endl;
                auto equal_mul_sub = cc->EvalMult(result_equal_comp[i], ai_sub_bi);
                // std::cout << "EvalAdd" << std::endl;
                comp_res = cc->EvalAdd(equal_mul_sub, result_sub[i]);
            }
        }

        end = std::chrono::system_clock::now();
        Plaintext plaintextDec, plaintextDec_comp_res, plaintextDec_comp_res_equal;

        cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec_comp_res);

        plaintextDec_comp_res->SetLength(encodedLength);
        std::vector<double> expectedOutput_compres;
        std::vector<std::complex<double>> finalResult = plaintextDec_comp_res->GetCKKSPackedValue();
        double err                                    = 0;
        int wrong                                     = 0;

        for (int i = 0; i < int(encodedLength); i++) {
            if (input_a_plain[i] != input_b_plain[i])
                expectedOutput_compres.push_back(static_cast<double>(input_a_plain[i] > input_b_plain[i]));
            else
                expectedOutput_compres.push_back(0.5);
            if (fabs(expectedOutput_compres[i] - finalResult[i].real()) > 0.5) {
                wrong++;
                std::cout << "input_a:" << input_a_plain[i] << ", input_b:" << input_b_plain[i]
                          << ", res :" << finalResult[i] << std::endl;
            }

            err += std::fabs(expectedOutput_compres[i] - finalResult[i].real());
            totalerr += err / encodedLength;
        }
        totaltime +=
            double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / double(encodedLength));
        // std::cout << "number of levels remaining before comp: " << multDepth - ciphertext_a[0]->GetLevel() << std::endl;
        // std::cout << "number of levels remaining after comp: " << multDepth - comp_res->GetLevel() << std::endl;
        // std::cout << "---------------------------------compres---------------------------------" << std::endl;
        if (wrong > 0) {
            std::cout << "avg error: 2^" << std::log2(err / encodedLength) << std::endl;
            std::cout << "wrong number:" << wrong << std::endl;
        }

        // std::cout << encodedLength << "slots total time: "
        //           << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / double(encodedLength)
        //           << "ms" << std::endl;
        // std::cout << encodedLength << "slots amortize time: "
        //           << double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) /
        //                  double(encodedLength)
        //           << "ms" << std::endl
        //           << std::endl;
    }
    std::cout << encodedLength << "slots amortize time: " << double((totaltime / 10.0)) << "ms" << std::endl
              << std::endl;
    std::cout << "avg error: 2^" << std::log2(totalerr / 10) << std::endl;

    // for (int j = 0; j < 8; j++) {
    //     std::cout << "input_a:" << input_a_plain[j] << ", input_b:" << input_b_plain[j] << ", res :" << finalResult[j]
    //               << std::endl;
    // }
}

void Eval_modular_Strictly_greater_test(uint32_t plain_bits, uint32_t block) {
    std::cout << plain_bits << " plain_bits " << block << " Strictly modular greater test " << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    // We set a smaller ring dimension to improve performance for this example.
    // In production environments, the security level should be set to
    // HEStd_128_classic, HEStd_192_classic, or HEStd_256_classic for 128-bit, 192-bit,
    // or 256-bit security, respectively.
    parameters.SetSecurityLevel(HEStd_128_classic);
    // parameters.SetRingDim(1 << 16);
#if NATIVEINT == 128
    usint scalingModSize = 78;
    usint firstModSize   = 89;
#else
    usint scalingModSize = 50;
    usint firstModSize   = 60;
#endif
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetScalingModSize(scalingModSize);
    parameters.SetFirstModSize(firstModSize);
    parameters.SetRingDim(65536 * 2);
    uint32_t polyDegree = 59;
    uint32_t multDepth  = 20 + block;
    // uint32_t multDepth = 50;
    // Choosing a higher degree yields better precision, but a longer runtime.
    // if (plain_bits >= 8 && block > 1)
    //     polyDegree = 119;
    // else if (plain_bits >= 8 && block == 1)
    //     polyDegree = 119;
    // else {
    //     polyDegree = 27;
    //     multDepth  = 20;
    // }

    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.
    cc->Enable(ADVANCEDSHE);
    double precision  = (1 << (plain_bits - 1) - 1);
    double lowerBound = -precision;
    double upperBound = precision;
    double bound      = 3;
    auto keyPair      = cc->KeyGen();
    // We need to generate mult keys to run Chebyshev approximations.
    cc->EvalMultKeyGen(keyPair.secretKey);
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    // std::uniform_int_distribution<int> message_part_16_1(-2048, -1);
    std::uniform_int_distribution<int> message(0, precision);
    usint ringDim = cc->GetRingDimension();
    int length    = ringDim / 2;

    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;
    std::cout << "using polyDegree " << polyDegree << std::endl;
    std::vector<std::vector<std::complex<double>>> input_a_mod, input_b_mod;
    std::vector<double> input_a_plain(length, 0), input_b_plain(length, 0);
    //a0,a1...  b0,b1...

    for (int i = 0; i < block; i++) {
        std::vector<std::complex<double>> a(length), b(length);
        for (int j = 0; j < length; j++) {
            a[j] = double(message(engine));
            b[j] = double(message(engine));
        }
        input_a_mod.push_back(a);
        input_b_mod.push_back(b);
    }

    for (int i = 0; i < length; i++) {
        for (int j = 0; j < block; j++) {
            input_a_plain[i] += input_a_mod[j][i].real() * std::pow(2, j * plain_bits);
            input_b_plain[i] += input_b_mod[j][i].real() * std::pow(2, j * plain_bits);
        }
    }
    // for (int j = 0; j < 8; j++) {
    //     std::cout << "input_a:" << input_a_plain[j] << ", input_b:" << input_b_plain[j] << std::endl;
    // }
    double totalerr      = 0;
    double totaltime     = 0;
    size_t encodedLength = input_a_plain.size();
    double num_test      = 1;
    for (int num = 0; num < int(num_test); num++) {
        std::vector<Plaintext> plaintext_a, plaintext_b;
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_a, ciphertext_b;

        std::chrono::system_clock::time_point start, end;
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_sub(block), result_sub(block), result_equal_comp(block),
            ciphertext_comp_greater(block);
        for (int i = 0; i < block; i++) {
            Plaintext plain_a = cc->MakeCKKSPackedPlaintext(input_a_mod[i]);
            Plaintext plain_b = cc->MakeCKKSPackedPlaintext(input_b_mod[i]);
            ciphertext_a.push_back(cc->Encrypt(keyPair.publicKey, plain_a));
            ciphertext_b.push_back(cc->Encrypt(keyPair.publicKey, plain_b));
        }
        start = std::chrono::system_clock::now();
        for (int i = block - 1; i >= 0; i--) {
            //comp(a>b) 1, 0.5, 0

            ciphertext_sub[i] = cc->EvalSub(ciphertext_a[i], ciphertext_b[i]);
            result_sub[i]     = cc->EvalSign(ciphertext_sub[i], bound, lowerBound, upperBound, polyDegree);
            Homround(result_sub[i]);

            //comp(a==b) 1, 0

            auto ciphertext_sub_neg = cc->EvalNegate(result_sub[i]);
            auto result_sub_neg     = cc->EvalAdd(ciphertext_sub_neg, 1.0);

            auto result_equal          = cc->EvalMult(result_sub_neg, result_sub[i]);
            result_equal_comp[i]       = MultByInteger(result_equal, 2.0);
            ciphertext_comp_greater[i] = cc->EvalSub(result_sub[i], result_equal_comp[i]);
            result_equal_comp[i]       = MultByInteger(result_equal_comp[i], 2.0);
            // Homround(result_equal_comp[i]);
            // Homround(ciphertext_comp_greater[i]);
        }

        Ciphertext<lbcrypto::DCRTPoly> comp_res;
        for (int i = 0; i < block; i++) {
            // std::cout << "mux:" << i << std::endl;
            if (i == 0)
                comp_res = ciphertext_comp_greater[i];
            else {
                // std::cout << "EvalSub" << std::endl;
                auto ai_sub_bi = cc->EvalSub(comp_res, ciphertext_comp_greater[i]);
                // std::cout << "EvalMult" << std::endl;
                auto equal_mul_sub = cc->EvalMult(result_equal_comp[i], ai_sub_bi);
                // std::cout << "EvalAdd" << std::endl;
                comp_res = cc->EvalAdd(equal_mul_sub, ciphertext_comp_greater[i]);
            }
        }
        // Homround(comp_res);
        std::cout << "number of levels remaining after comp: " << multDepth - comp_res->GetLevel() << std::endl;
        end = std::chrono::system_clock::now();
        Plaintext plaintextDec, plaintextDec_comp_res, plaintextDec_comp_res_equal;

        cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec_comp_res);

        plaintextDec_comp_res->SetLength(encodedLength);
        std::vector<double> expectedOutput_compres;
        std::vector<std::complex<double>> finalResult = plaintextDec_comp_res->GetCKKSPackedValue();
        double err                                    = 0;
        int wrong                                     = 0;

        for (int i = 0; i < int(encodedLength); i++) {
            expectedOutput_compres.push_back(static_cast<double>(input_a_plain[i] > input_b_plain[i]));

            if (fabs(expectedOutput_compres[i] - finalResult[i].real()) > 0.5) {
                wrong++;
                std::cout << "input_a:" << input_a_plain[i] << ", input_b:" << input_b_plain[i]
                          << ", res :" << finalResult[i] << std::endl;
            }
            err += std::fabs(expectedOutput_compres[i] - finalResult[i].real());
        }
        totalerr += err / encodedLength;
        totaltime +=
            double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / double(encodedLength));

        if (wrong > 0) {
            std::cout << "avg error: 2^" << std::log2(err / encodedLength) << std::endl;
            std::cout << "wrong number:" << wrong << std::endl;
        }
    }

    std::cout << encodedLength << "slots amortize time: " << double((totaltime / num_test)) << "ms" << std::endl
              << std::endl;
    std::cout << "avg error: 2^" << std::log2(totalerr / num_test) << std::endl;

    // for (int j = 0; j < 8; j++) {
    //     std::cout << "input_a:" << input_a_plain[j] << ", input_b:" << input_b_plain[j] << ", res :" << finalResult[j]
    //               << std::endl;
    // }
}

void Eval_modular_Strictly_equal_test(uint32_t plain_bits, uint32_t block) {
    std::cout << plain_bits << " plain_bits " << block << " Strictly modular equal test " << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    // We set a smaller ring dimension to improve performance for this example.
    // In production environments, the security level should be set to
    // HEStd_128_classic, HEStd_192_classic, or HEStd_256_classic for 128-bit, 192-bit,
    // or 256-bit security, respectively.
    parameters.SetSecurityLevel(HEStd_128_classic);
    // parameters.SetRingDim(1 << 16);
#if NATIVEINT == 128
    usint scalingModSize = 78;
    usint firstModSize   = 89;
#else
    usint scalingModSize = 40;
    usint firstModSize   = 60;
#endif
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetScalingModSize(scalingModSize);
    parameters.SetFirstModSize(firstModSize);

    uint32_t polyDegree = 59;
    // uint32_t multDepth  = 23;
    uint32_t multDepth = 20 + block;
    // Choosing a higher degree yields better precision, but a longer runtime.
    // if (plain_bits >= 8 && block > 1)
    //     polyDegree = 247;
    // else if (plain_bits >= 8 && block == 1)
    //     polyDegree = 119;
    // else {
    //     polyDegree = 27;
    //     multDepth  = 20;
    // }

    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.
    cc->Enable(ADVANCEDSHE);
    double precision  = (1 << (plain_bits - 1) - 1);
    double lowerBound = -precision;
    double upperBound = precision;
    double bound      = 3;
    auto keyPair      = cc->KeyGen();
    // We need to generate mult keys to run Chebyshev approximations.
    cc->EvalMultKeyGen(keyPair.secretKey);
    const std::vector<DCRTPoly>& ckks_pk = keyPair.publicKey->GetPublicElements();
    std::cout << "Moduli chain of pk: " << std::endl;
    print_moduli_chain(ckks_pk[0]);
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    // std::uniform_int_distribution<int> message_part_16_1(-2048, -1);
    std::uniform_int_distribution<int> message(0, precision);
    usint ringDim = cc->GetRingDimension();
    int length    = ringDim / 2;

    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;
    std::cout << "using polyDegree " << polyDegree << std::endl;
    std::vector<std::vector<std::complex<double>>> input_a_mod, input_b_mod;
    std::vector<double> input_a_plain(length, 0), input_b_plain(length, 0);
    //a0,a1...  b0,b1...

    for (int i = 0; i < block; i++) {
        std::vector<std::complex<double>> a(length), b(length);
        for (int j = 0; j < length; j++) {
            a[j] = double(message(engine));
            b[j] = double(message(engine));
        }
        input_a_mod.push_back(a);
        input_b_mod.push_back(b);
    }

    for (int i = 0; i < length; i++) {
        for (int j = 0; j < block; j++) {
            input_a_plain[i] += input_a_mod[j][i].real() * std::pow(2, j * plain_bits);
            input_b_plain[i] += input_b_mod[j][i].real() * std::pow(2, j * plain_bits);
        }
    }
    // for (int j = 0; j < 8; j++) {
    //     std::cout << "input_a:" << input_a_plain[j] << ", input_b:" << input_b_plain[j] << std::endl;
    // }
    double totalerr      = 0;
    double totaltime     = 0;
    size_t encodedLength = input_a_plain.size();
    double num_test      = 1;
    for (int num = 0; num < int(num_test); num++) {
        std::vector<Plaintext> plaintext_a, plaintext_b;
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_a, ciphertext_b;

        std::chrono::system_clock::time_point start, end;
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_sub(block), result_sub(block), result_equal_comp(block),
            ciphertext_comp_greater(block);
        for (int i = 0; i < block; i++) {
            Plaintext plain_a = cc->MakeCKKSPackedPlaintext(input_a_mod[i]);
            Plaintext plain_b = cc->MakeCKKSPackedPlaintext(input_b_mod[i]);
            ciphertext_a.push_back(cc->Encrypt(keyPair.publicKey, plain_a));
            ciphertext_b.push_back(cc->Encrypt(keyPair.publicKey, plain_b));
        }
        start = std::chrono::system_clock::now();
        for (int i = block - 1; i >= 0; i--) {
            //comp(a>b) 1, 0.5, 0

            ciphertext_sub[i] = cc->EvalSub(ciphertext_a[i], ciphertext_b[i]);
            result_sub[i]     = cc->EvalSign(ciphertext_sub[i], bound, lowerBound, upperBound, polyDegree);
            Homround(result_sub[i]);

            //comp(a==b) 1, 0

            auto ciphertext_sub_neg = cc->EvalNegate(result_sub[i]);
            auto result_sub_neg     = cc->EvalAdd(ciphertext_sub_neg, 1.0);

            auto result_equal          = cc->EvalMult(result_sub_neg, result_sub[i]);
            result_equal_comp[i]       = MultByInteger(result_equal, 4.0);
            ciphertext_comp_greater[i] = cc->EvalSub(result_sub[i], result_equal_comp[i]);
        }

        Ciphertext<lbcrypto::DCRTPoly> comp_res;
        for (int i = 0; i < block; i++) {
            // std::cout << "mux:" << i << std::endl;
            if (i == 0)
                comp_res = result_equal_comp[i];
            else {
                comp_res = cc->EvalMult(comp_res, result_equal_comp[i]);
            }
        }
        Homround(comp_res);
        end = std::chrono::system_clock::now();
        std::cout << "number of levels remaining after comp: " << multDepth - comp_res->GetLevel() << std::endl;
        Plaintext plaintextDec, plaintextDec_comp_res, plaintextDec_comp_res_equal;

        cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec_comp_res);

        plaintextDec_comp_res->SetLength(encodedLength);
        std::vector<double> expectedOutput_compres;
        std::vector<std::complex<double>> finalResult = plaintextDec_comp_res->GetCKKSPackedValue();
        double err                                    = 0;
        int wrong                                     = 0;

        for (int i = 0; i < int(encodedLength); i++) {
            expectedOutput_compres.push_back(static_cast<double>(input_a_plain[i] == input_b_plain[i]));

            if (fabs(expectedOutput_compres[i] - finalResult[i].real()) > 0.1) {
                wrong++;
                std::cout << "input_a:" << input_a_plain[i] << ", input_b:" << input_b_plain[i]
                          << ", res :" << finalResult[i] << std::endl;
            }
            // std::cout << "input_a:" << input_a_plain[i] << ", input_b:" << input_b_plain[i]
            //           << ", res :" << finalResult[i] << std::endl;
            err += std::fabs(expectedOutput_compres[i] - finalResult[i].real());
        }
        totalerr += err / encodedLength;
        totaltime +=
            double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / double(encodedLength));

        if (wrong > 0) {
            std::cout << "avg error: 2^" << std::log2(err / encodedLength) << std::endl;
            std::cout << "wrong number:" << wrong << std::endl;
        }
    }
    std::cout << encodedLength << "slots amortize time: " << double((totaltime / num_test)) << "ms" << std::endl
              << std::endl;
    std::cout << "avg error: 2^" << std::log2(totalerr / num_test) << std::endl;

    // for (int j = 0; j < 8; j++) {
    //     std::cout << "input_a:" << input_a_plain[j] << ", input_b:" << input_b_plain[j] << ", res :" << finalResult[j]
    //               << std::endl;
    // }
}

void EvalequalExample(uint32_t plain_bits) {
    std::cout << "\nplain_bits equal comp: " << plain_bits << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    // We set a smaller ring dimension to improve performance for this example.
    // In production environments, the security level should be set to
    // HEStd_128_classic, HEStd_192_classic, or HEStd_256_classic for 128-bit, 192-bit,
    // or 256-bit security, respectively.
    parameters.SetSecurityLevel(HEStd_128_classic);
    // parameters.SetRingDim(1 << 16);
#if NATIVEINT == 128
    usint scalingModSize = 78;
    usint firstModSize   = 89;
#else
    usint scalingModSize = 59;
    usint firstModSize   = 60;
#endif
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetScalingModSize(scalingModSize);
    parameters.SetFirstModSize(firstModSize);
    std::vector<uint32_t> levelBudget = {4, 4};
    // Choosing a higher degree yields better precision, but a longer runtime.
    uint32_t polyDegree = 1007;

    // The multiplicative depth depends on the polynomial degree.
    // See the FUNCTION_EVALUATION.md file for a table mapping polynomial degrees to multiplicative depths.
    uint32_t multDepth = 25;
    // uint32_t levelsAvailableAfterBootstrap = 30;
    // uint32_t numIterations                 = 1;
    // usint multDepth =
    //     levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist) + (numIterations - 1);
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(FHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.
    cc->Enable(ADVANCEDSHE);

    double precision  = (1 << (plain_bits - 1) - 1);
    double lowerBound = -precision;
    double upperBound = precision;
    double bound      = 3;
    auto keyPair      = cc->KeyGen();
    usint ringDim     = cc->GetRingDimension();
    int length        = ringDim / 2;
    // We need to generate mult keys to run Chebyshev approximations.
    cc->EvalMultKeyGen(keyPair.secretKey);
    // cc->EvalBootstrapSetup(levelBudget);
    // cc->EvalBootstrapKeyGen(keyPair.secretKey, length);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    // std::uniform_int_distribution<int> message_part_16_1(-2048, -1);
    std::uniform_int_distribution<int> message(0, precision);

    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;
    std::vector<std::complex<double>> input_a(length), input_b(length);

    for (int i = 0; i < length / 2; i++) {
        input_a[i] = double(message(engine));
        input_b[i] = double(message(engine));
    }
    size_t encodedLength  = input_a.size();
    Plaintext plaintext_a = cc->MakeCKKSPackedPlaintext(input_a);
    Plaintext plaintext_b = cc->MakeCKKSPackedPlaintext(input_b);
    auto ciphertext_a     = cc->Encrypt(keyPair.publicKey, plaintext_a);
    auto ciphertext_b     = cc->Encrypt(keyPair.publicKey, plaintext_b);
    //comp(a>b),  a-b

    std::cout << "number of levels remaining before comp: " << multDepth - ciphertext_a->GetLevel() << std::endl;
    std::chrono::system_clock::time_point start, end;
    start               = std::chrono::system_clock::now();
    auto ciphertext_sub = cc->EvalSub(ciphertext_a, ciphertext_b);
    // auto ciphertext_sub       = cc->EvalMult(ciphertext_sub_scale, double(1/precision);
    auto result_sub = cc->EvalSign(ciphertext_sub, bound, lowerBound, upperBound, polyDegree);
    std::cout << "number of levels remaining before Homround: " << multDepth - result_sub->GetLevel() << std::endl;
    Homround(result_sub);
    std::cout << "number of levels remaining after Homround: " << multDepth - result_sub->GetLevel() << std::endl;
    // auto ciphertext_greater_equal_2 = cc->EvalMult(result_sub, 2.0);

    //comp(a==b)

    auto ciphertext_sub_neg = cc->EvalNegate(result_sub);
    auto result_sub_neg     = cc->EvalAdd(ciphertext_sub_neg, 1.0);

    auto result_equal   = cc->EvalMult(result_sub_neg, result_sub);
    auto result_equal_4 = MultByInteger(result_equal, 2.0);
    // auto result_equal_4 = cc->EvalMult(result_equal, 4.0);
    // Homround(result_equal_4);

    // result_equal_4 = cc->EvalBootstrap(result_equal_4);

    auto ciphertext_comp_greater = cc->EvalSub(result_sub, result_equal_4);
    end                          = std::chrono::system_clock::now();
    std::cout << "number of levels remaining after equal comp: " << multDepth - result_equal_4->GetLevel() << std::endl;
    // auto ciphertext_comp_greater = cc->EvalSub(ciphertext_greater_equal_2, result_equal_4);
    // ciphertext_comp_greater      = cc->EvalMult(ciphertext_comp_greater, 0.5);
    // Homround(ciphertext_comp_greater);
    std::cout << "number of levels remaining after greater comp: " << multDepth - ciphertext_comp_greater->GetLevel()
              << std::endl;
    auto ciphertext_res_greater = cc->EvalMult(ciphertext_a, ciphertext_comp_greater);
    auto ciphertext_res_equal   = cc->EvalMult(ciphertext_a, result_equal_4);

    // auto ciphertext_neg_sub        = cc->EvalNegate(result_total);
    // auto ciphertext_neg_sub_puls_1 = cc->EvalAdd(ciphertext_neg_sub, 1);
    // auto ciphertext_b_comp         = cc->EvalMult(ciphertext_b, ciphertext_neg_sub_puls_1);
    // auto ciphertext_res            = cc->EvalAdd(ciphertext_b_comp, ciphertext_a_comp);

    Plaintext plaintextDec, plaintextDec_result_equal, plaintextDec_result_greater, plaintextDec_comp_res,
        plaintextDec_greater_res;

    cc->Decrypt(keyPair.secretKey, ciphertext_res_greater, &plaintextDec_greater_res);
    cc->Decrypt(keyPair.secretKey, ciphertext_res_equal, &plaintextDec_comp_res);
    cc->Decrypt(keyPair.secretKey, result_equal_4, &plaintextDec_result_equal);
    cc->Decrypt(keyPair.secretKey, ciphertext_comp_greater, &plaintextDec_result_greater);
    plaintextDec_comp_res->SetLength(encodedLength);
    plaintextDec_result_equal->SetLength(encodedLength);
    plaintextDec_result_greater->SetLength(encodedLength);
    plaintextDec_greater_res->SetLength(encodedLength);
    std::vector<std::complex<double>> Result_comp_equal   = plaintextDec_result_equal->GetCKKSPackedValue();
    std::vector<std::complex<double>> Result_comp_greater = plaintextDec_result_greater->GetCKKSPackedValue();
    std::vector<std::complex<double>> Result_equal        = plaintextDec_comp_res->GetCKKSPackedValue();
    std::vector<std::complex<double>> Result_greater      = plaintextDec_greater_res->GetCKKSPackedValue();
    double err                                            = 0;
    double err_greater                                    = 0;
    double err_comp_equal                                 = 0;
    double err_comp_greater                               = 0;
    int wrong                                             = 0;
    std::vector<double> expectedOutput_Result_equal, expectedOutput_Result_greater, expectedOutput_Comp_equal,
        expectedOutput_Comp_greater;
    for (int i = 0; i < (int)encodedLength; i++) {
        // expectedOutput_compres.push_back(static_cast<double>(input_a[i].real() == input_b[i].real()));
        expectedOutput_Result_equal.push_back(
            static_cast<double>(input_a[i].real() == input_b[i].real() ? input_a[i].real() : 0));
        expectedOutput_Result_greater.push_back(input_a[i].real() > input_b[i].real() ? input_a[i].real() : 0);
        expectedOutput_Comp_equal.push_back(input_a[i].real() == input_b[i].real() ? 1 : 0);
        expectedOutput_Comp_greater.push_back(input_a[i].real() > input_b[i].real() ? 1 : 0);
        // std::cout << "input_a:" << input_a[i].real() << ", input_b:" << input_b[i].real()
        //           << ", res :" << Result_comp_greater[i].real() << std::endl;
        // if (fabs(expectedOutput_Result_equal[i] - Result_equal[i].real()) > 0.1 ||
        //     fabs(expectedOutput_Result_greater[i] - Result_greater[i].real()) > 0.1) {
        //     wrong++;
        //     // std::cout << "input_a:" << input_a[i].real() << ", input_b:" << input_b[i].real()
        //     //           << ", Result_comp_greater:" << Result_comp_greater[i].real()
        //     //           << ", Result_comp_equal:" << Result_comp_equal[i].real() << ", res :" << Result_equal[i].real()
        //     //           << ", res_greater :" << Result_greater[i].real() << std::endl;
        // }
        if (std::fabs(expectedOutput_Comp_equal[i] - Result_comp_equal[i].real()) > 0.5 ||
            std::fabs(expectedOutput_Comp_greater[i] - Result_comp_greater[i].real()) > 0.5) {
            wrong++;
            // std::cout << "input_a:" << input_a[i].real() << ", input_b:" << input_b[i].real()
            //           << ", Result_comp_greater:" << Result_comp_greater[i].real()
            //           << ", Result_comp_equal:" << Result_comp_equal[i].real() << ", res :" << Result_equal[i].real()
            //           << ", res_greater :" << Result_greater[i].real() << std::endl;
        }
        err += std::fabs(expectedOutput_Result_equal[i] - Result_equal[i].real());
        err_greater += std::fabs(expectedOutput_Result_greater[i] - Result_greater[i].real());
        err_comp_equal += std::fabs(expectedOutput_Comp_equal[i] - Result_comp_equal[i].real());
        err_comp_greater += std::fabs(expectedOutput_Comp_greater[i] - Result_comp_greater[i].real());
    }

    std::cout << "---------------------------------compres---------------------------------" << std::endl;
    std::cout << "error equal res: 2^" << std::log2(err / encodedLength) << std::endl;
    std::cout << "error equal comp: 2^" << std::log2(err_comp_equal / encodedLength) << std::endl;
    std::cout << "error greater res: 2^" << std::log2(err_greater / encodedLength) << std::endl;
    std::cout << "error greater comp: 2^" << std::log2(err_comp_greater / encodedLength) << std::endl;
    std::cout << "wrong number:" << wrong << std::endl;
    std::cout << encodedLength << "slots equal total time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
    std::cout << encodedLength << "slots amortize time: "
              << double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) /
                     double(encodedLength)
              << "ms" << std::endl
              << std::endl;
}

int main(int argc, char* argv[]) {
    // EvalLogisticExample();
    // EvalFunctionExample();
    // EvalSatExample();
    // EvalgreaterExample(8);
    // EvalSignExample(9);
    // Eval_modular_Strictly_greater_test(8, 2);
    Eval_modular_Strictly_equal_test(8, 2);
    // for (int i = 1; i < 9; i++) {
    //     Eval_modular_greater_test(8, i);
    // }
    // for (int i = 1; i < 9; i++) {
    //     Eval_modular_Strictly_greater_test(8, i);
    // }
    // for (int i = 1; i < 9; i++) {
    //     Eval_modular_Strictly_equal_test(8, i);
    // }

    // EvalgreaterExample(16);
    // for (int i = 8; i < 36; i++)
    // EvalcompExample(i);
    // EvalmaxExample(7);
    // EvalequalExample(8);
    // EvalDeadZoneExample();
    // EvalRelayExample();
    // EvalDeadRelayExample();
    // EvalDeadSatExample();
    return 0;
}
