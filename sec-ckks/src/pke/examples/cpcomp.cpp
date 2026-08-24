//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2025, CEA-List
//
// All rights reserved.
//
// Author TPOC: jules.dumezy@cea.fr
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

#define PROFILE
#include "math/hermite.h"
#include "openfhe.h"
#include <chrono>
#include <functional>
using namespace lbcrypto;
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
void FuncBootstrapExample(std::function<double(double)> f, int bits, int p);
void CarrymulExample(int bits, int p, int k);

void MVFuncBootstrap_kx_parallel_Example(int bits, int p, int opt_depth) {
    std::cout << "---------------------EvalMVFBT Parallel Test--------------------" << std::endl << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    // For secure computations, use parameters.SetSecurityLevel(HEStd_128_classic);
    // and remove the two lines below.
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 11);

    parameters.SetNumLargeDigits(3);
    parameters.SetKeySwitchTechnique(HYBRID);

#if NATIVEINT == 128
    OPENFHE_THROW("NATIVEINT 128 is not supported, use 64.");
#endif
    // Currently, only FIXEDMANUAL mode is supported for functional bootstrapping.
    ScalingTechnique rescaleTech = FIXEDMANUAL;
    usint dcrtBits               = 49;
    usint firstMod               = 50;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    std::vector<uint32_t> levelBudget = {4, 2};

    std::vector<uint32_t> bsgsDim = {0, 0};

    // usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 4 + bits + 2;
    usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 2 + opt_depth + 3;
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    cc->Enable(FBTS);

    usint ringDim = cc->GetRingDimension();
    std::cout << "CKKS scheme is using ring dimension " << ringDim << " depth: " << depth << std::endl << std::endl;

    // Precomputations for bootstrapping
    int numSlots = ringDim / 2;
    cc->EvalFuncBootstrapSetup(levelBudget, bsgsDim, numSlots, bits);

    // Key Generation
    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);
    std::vector<int32_t> indexList;
    for (int i = 0; i < std::log(p); i++) {
        indexList.push_back(pow(2, i));
    }
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, p - 1);
    cc->EvalRotateKeyGen(keyPair.secretKey, indexList);
    // Input vector [0, 1, ..., p-2, p-1]
    // std::complex<double> input = 1;
    std::vector<std::complex<double>> x;
    for (int i = 0; i < numSlots; ++i) {
        x.push_back(std::complex<double>((i) % p, (i + 1) % p));
    }

    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(x, 1, depth - levelBudget[1] - 1, nullptr, numSlots);
    ptxt->SetLength(numSlots);
    // std::cout << "Input = " << x << std::endl;
    // Plaintext ptxt1            = cc->MakeCKKSPackedPlaintext(1.0, 1, depth - levelBudget[1] - 1, nullptr, numSlots);
    Ciphertext<DCRTPoly> ctxt = cc->Encrypt(keyPair.publicKey, ptxt);
    DCRTPoly& poly            = ctxt->GetElements()[0];
    std::cout << "ctxt Modulus:" << std::endl;
    print_moduli_chain(poly);
    // Order of the hermite interpolation, should be between 1 and 3
    int hermite_order = 1;
    std::vector<std::complex<double>> coe_1, coe_2, coe_3, coe_4, coe_rs_4bit;
    std::chrono::system_clock::time_point start, end;
    start    = std::chrono::system_clock::now();
    auto f_1 = [](double x) -> double {
        return (double)((int)x % 16);
    };
    auto f_2 = [](double x) -> double {
        return (double)((((int)x % 256) - ((int)x % 16)) / 16);
    };
    auto f_3 = [](double x) -> double {
        return (double)((((int)x % 4096) - ((int)x % 256)) / 256);
    };
    auto f_4 = [](double x) -> double {
        return (double)(((int)x - ((int)x % 4096)) / 4096);
    };
    auto f_rs_4bit = [](double x) -> double {
        return (double)(((int)x - ((int)x % 16)) / 16);
    };

    coe_1       = hermiteInterpOrder1Coeff([f_1](double x) -> double { return 0.5 * f_1(x); }, p);
    coe_2       = hermiteInterpOrder1Coeff([f_2](double x) -> double { return 0.5 * f_2(x); }, p);
    coe_3       = hermiteInterpOrder1Coeff([f_3](double x) -> double { return 0.5 * f_3(x); }, p);
    coe_4       = hermiteInterpOrder1Coeff([f_4](double x) -> double { return 0.5 * f_4(x); }, p);
    coe_rs_4bit = hermiteInterpOrder1Coeff([f_rs_4bit](double x) -> double { return 0.5 * f_rs_4bit(x); }, p);
    std::vector<std::vector<std::complex<double>>> coe_1234 = {coe_1, coe_2, coe_3, coe_4};
    std::vector<std::vector<std::complex<double>>> coe_1_rs = {coe_1, coe_rs_4bit};
    end                                                     = std::chrono::system_clock::now();
    double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

    // We evaluate f(ctxt) while bootstrapping
    auto use_img = false;
    std::cout << "use_img: " << use_img << std::endl;
    std::vector<Ciphertext<DCRTPoly>> ctxtResult;
    start = std::chrono::system_clock::now();
    cc->EvalMultInPlace(ctxt, 1);
    cc->ModReduceInPlace(ctxt);
    auto num_test = 1;
    for (int i = 0; i < num_test; i++)
        ctxtResult = cc->EvalMVFuncBootstrap_kx_parallel_opt(ctxt, p, hermite_order, coe_1234, opt_depth, use_img);

    end       = std::chrono::system_clock::now();
    totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "level after fbt:" << depth - ctxtResult[0]->GetLevel() << std::endl;

    std::cout << "EvalFuncBootstrap: " << totaltime / num_test << " ms" << std::endl;
    std::cout << "EvalFuncBootstrap amortize: " << totaltime / (num_test * ringDim) << " ms" << std::endl;

    std::cout << "noise deg after fbt:" << ctxtResult[0]->GetNoiseScaleDeg() << std::endl;

    Plaintext result;
    cc->Decrypt(keyPair.secretKey, ctxtResult[0], &result);
    result->SetLength(numSlots);
    std::vector<std::complex<double>> finalResult = result->GetCKKSPackedValue();
    // if (int(x[0].real()) % 16 != int(finalResult[0].real())) {
    //     std::cout << "x: " << x[0].real() << ", int(x[i]) % 16: " << int(x[0].real()) % 16
    //               << ", finalResult[i].real(): " << finalResult[0].real()
    //               << ", int(finalResult[i].real()) % 16: " << int(finalResult[0].real()) % 16 << std::endl;
    // }
    // std::cout << "x: " << x[1].real() << ", int(x[i]) % 256: " << int(x[1].real()) % 256
    //           << ", finalResult[i].real(): " << finalResult[1].real() << std::endl;
    int error_num = 0;
    for (int i = 0; i < numSlots; i++) {
        // std::cout << "x.real(): " << x[i].real() << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
        if (f_1(x[i].real()) != int(finalResult[i].real())) {
            error_num++;
            std::cout << "x.real(): " << x[i].real() << ",f_1(x.real): " << f_1(x[i].real())
                      << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
        }
        // if (int(x[i].imag()) % 16 != int(finalResult[i].imag())) {
        //     error_num++;
        //     std::cout << "x.imag(): " << x[i].imag() << ", int(x[i]).imag() % 16: " << int(x[i].imag()) % 16
        //               << ", finalResult[i].imag(): " << finalResult[i].imag()
        //               << ", int(finalResult[i].imag()) % 16: " << int(finalResult[i].imag()) % 16 << std::endl;
        // }
        // error_num++;
        // std::cout << "i: " << i << ", x: " << x[i].real() << ", int(x[i]) % 16: " << int(x[i].real()) % 16
        //           << ", finalResult[i].real(): " << finalResult[i].real()
        //           << ", int(finalResult[i].real()) % 16: " << int(finalResult[i].real()) % 16 << std::endl;

        // }
    }
    std::cout << "error_num:" << error_num << std::endl;
}

void CP_comparison_Example(int bits, int p, int opt_depth, int comprange) {
    std::cout << "---------------------EvalMV CP FBT Parallel Test--------------------" << std::endl << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    // For secure computations, use parameters.SetSecurityLevel(HEStd_128_classic);
    // and remove the two lines below.
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 16);

    parameters.SetNumLargeDigits(3);
    parameters.SetKeySwitchTechnique(HYBRID);

#if NATIVEINT == 128
    OPENFHE_THROW("NATIVEINT 128 is not supported, use 64.");
#endif
    // Currently, only FIXEDMANUAL mode is supported for functional bootstrapping.
    ScalingTechnique rescaleTech = FIXEDMANUAL;
    usint dcrtBits               = 52;
    usint firstMod               = 53;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    std::vector<uint32_t> levelBudget = {4, 2};

    std::vector<uint32_t> bsgsDim = {0, 0};

    // usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 4 + bits + 2;
    usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 2 + opt_depth + 3;
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    cc->Enable(FBTS);

    usint ringDim = cc->GetRingDimension();
    std::cout << "CKKS scheme is using ring dimension " << ringDim << " depth: " << depth << std::endl << std::endl;

    // Precomputations for bootstrapping
    int numSlots = ringDim / 2;
    cc->EvalFuncBootstrapSetup(levelBudget, bsgsDim, numSlots, bits);

    // Key Generation
    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);
    std::vector<int32_t> indexList;
    for (int i = 0; i < std::log(p); i++) {
        indexList.push_back(pow(2, i));
    }
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, p - 1);
    cc->EvalRotateKeyGen(keyPair.secretKey, indexList);
    // Input vector [0, 1, ..., p-2, p-1]
    // std::complex<double> input = 1;
    std::vector<std::complex<double>> x;
    for (int i = 0; i < numSlots; ++i) {
        x.push_back(std::complex<double>((i) % p, (i + 1) % p));
    }

    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(x, 1, depth - levelBudget[1], nullptr, numSlots);
    ptxt->SetLength(numSlots);
    // std::cout << "Input = " << x << std::endl;
    // Plaintext ptxt1            = cc->MakeCKKSPackedPlaintext(1.0, 1, depth - levelBudget[1] - 1, nullptr, numSlots);
    Ciphertext<DCRTPoly> ctxt = cc->Encrypt(keyPair.publicKey, ptxt);
    DCRTPoly& poly            = ctxt->GetElements()[0];
    std::cout << "ctxt Modulus:" << std::endl;
    print_moduli_chain(poly);
    // Order of the hermite interpolation, should be between 1 and 3
    int hermite_order = 1;
    std::vector<std::vector<std::complex<double>>> coe_comp;
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    for (int i = 0; i < comprange; i++) {
        auto f = [i](double x) -> double {
            return (int)x - i == 0;
        };

        coe_comp.push_back(hermiteInterpOrder1Coeff([f](double x) -> double { return 0.5 * f(x); }, p));
    }
    // auto f = [](double x) -> double {
    //     return (int)x - 10 == 0;
    // };
    // coe_comp.push_back(hermiteInterpOrder1Coeff([f](double x) -> double { return 0.5 * f(x); }, p));
    end              = std::chrono::system_clock::now();
    double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

    // We evaluate f(ctxt) while bootstrapping
    auto use_img = true;
    std::cout << "use_img: " << use_img << std::endl;
    std::vector<Ciphertext<DCRTPoly>> ctxtResult;
    start = std::chrono::system_clock::now();
    // cc->EvalMultInPlace(ctxt, 1);
    // cc->ModReduceInPlace(ctxt);
    auto num_test = 5;
    for (int i = 0; i < num_test; i++) {
        std::cout << "test: " << i << std::endl;
        ctxtResult = cc->EvalMVFuncBootstrap_kx_parallel_opt(ctxt, p, hermite_order, coe_comp, opt_depth, use_img);
    }
    end       = std::chrono::system_clock::now();
    totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "level after fbt:" << depth - ctxtResult[0]->GetLevel() << std::endl;

    std::cout << "Eval " << comprange << " CP comparisons time: " << totaltime / num_test << " ms" << std::endl;
    std::cout << "Eval " << comprange << " CP comparisons amortize: " << totaltime / (num_test * ringDim * comprange)
              << " ms " << std::endl;

    std::cout << "noise deg after fbt:" << ctxtResult[0]->GetNoiseScaleDeg() << std::endl;
    for (int func = 0; func < int(coe_comp.size()); func++) {
        Plaintext result;
        cc->Decrypt(keyPair.secretKey, ctxtResult[func], &result);
        result->SetLength(numSlots);
        std::vector<std::complex<double>> finalResult = result->GetCKKSPackedValue();
        int error_num                                 = 0;
        for (int i = 0; i < numSlots; i++) {
            if ((int(x[i].real()) == func) != int(finalResult[i].real())) {
                error_num++;
                std::cout << "x.real(): " << x[i].real() << ",func(): " << func
                          << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
                // std::cout << "x.real(): " << x[i].real() << ",func: " << func
                //           << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
            }
            if ((int(x[i].imag()) == func) != int(finalResult[i].imag())) {
                error_num++;
                std::cout << "x.imag(): " << x[i].imag() << ",func: " << func
                          << ", finalResult[i].imag(): " << finalResult[i].imag() << std::endl;
                // std::cout << "x.real(): " << x[i].imag() << ",func: " << func
                //   << ", finalResult[i].real(): " << finalResult[i].imag() << std::endl;
            }
            // error_num++;
            // std::cout << "i: " << i << ", x: " << x[i].real() << ", int(x[i]) % 16: " << int(x[i].real()) % 16
            //           << ", finalResult[i].real(): " << finalResult[i].real()
            //           << ", int(finalResult[i].real()) % 16: " << int(finalResult[i].real()) % 16 << std::endl;

            // }
        }
        std::cout << "comp with " << func << " error_num: " << error_num << std::endl;
    }
}

int main(int argc, char* argv[]) {
    int bits      = 8;
    int opt_depth = 3;
    int p         = pow(2, bits);
    int comprange = 8;

    CP_comparison_Example(bits, p, opt_depth, comprange);
}