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
#include <fstream>
#include <iostream>
#include <vector>
#include <complex>
// #include "cryptocontext.h"
#include <chrono>
using namespace lbcrypto;
using namespace std::complex_literals;
// Coefficients for the function std::exp(1i * Pi/2.0 * x) in [-16, 16] of degree 46
// Need two double-angle iterations to get std::exp(1i * 2Pi * x)
static const inline std::vector<std::complex<double>> coeff_exp_16_double_46{
    0.22393566906777329084,
    4.72e-18 - 0.22176384914036376128i,
    0.24158307546266130639 - 7.09e-18i,
    -0.18331470851313935722i,
    0.28534623846463541552 + 5.91e-18i,
    -2.362e-18 - 0.092486179824488604084i,
    0.32214532018151836867 - 4.72e-18i,
    0.061326880477941263237i,
    0.28798365357787297780 - 4.72e-18i,
    1.417e-17 + 0.24466296846427090794i,
    0.112756709876059055264 - 7.087e-18i,
    0.33439190718203853914i,
    -0.17995397739265364678 + 2.36e-18i,
    9.45e-18 + 0.16254851699551037258i,
    -0.34811157721125479680 - 2.36e-18i,
    2.95e-18 - 0.22527723082929962395i,
    -0.079206690817228031509 - 3.543e-18i,
    9.45e-18 - 0.32612632178540534866i,
    0.36198254675123692214 - 4.13e-18i,
    -4.72e-18 + 0.19237548287066727482i,
    0.071116210979946081761 - 9.449e-18i,
    -7.09e-18 + 0.30556044798491310832i,
    -0.43951407397686892420 + 2.36e-18i,
    4.72e-18 - 0.46389876376571992367i,
    0.40955141151976887093 - 3.54e-18i,
    0.31828681535788988510i,
    -0.22366008829505132360 - 4.72e-18i,
    4.72e-18 - 0.14446909676096356123i,
    0.086745018497585937856 - 7.087e-18i,
    2.362e-18 + 0.048813481993878263254i,
    -0.025904132260782003483 - 2.362e-18i,
    2.3622e-18 - 0.0130280784432669702322i,
    0.0062348555293589751417 - 7.0865e-18i,
    -1.41731e-17 + 0.0028488507881147717878i,
    -0.00124638777412581667689 + 6.49599e-18i,
    -1.77163e-18 - 0.00052341839132928819674i,
    0.00021144315086689963591 - 2.36218e-18i,
    -4.724353e-18 + 0.000082321616249500458072i,
    -0.000030941853907365078754,
    -9.4487066e-18 - 0.0000112448146643207249387i,
    3.9566691190946415917e-6 + 4.7243533e-18i,
    -4.72435330e-18 + 1.34965353351963424859e-6i,
    -4.4681665430474632002e-7 - 7.08652994e-18i,
    3.54326497e-18 - 1.4370869470312489006e-7i,
    4.4978580236218894379e-8 + 2.362176648e-18i,
    -4.7243532963e-18 + 1.35960019696524527036e-8i,
    -4.3910916203385461663e-9 - 4.7243532963e-18i};

void SaveCoefficients(const std::string& filename, const std::vector<std::complex<double>>& coeffs) {
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile.is_open()) {
        throw std::runtime_error("无法打开文件进行写入: " + filename);
    }

    
    size_t size = coeffs.size();
    outfile.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

    
    
    
    if (size > 0) {
        outfile.write(reinterpret_cast<const char*>(coeffs.data()), size * sizeof(std::complex<double>));
    }

    outfile.close();
    std::cout << "[IO] 已保存 " << size << " 个系数到 " << filename << std::endl;
}

std::vector<std::complex<double>> LoadCoefficients(const std::string& filename) {
    std::ifstream infile(filename, std::ios::binary);
    if (!infile.is_open()) {
        throw std::runtime_error("无法打开文件进行读取: " + filename);
    }

    
    size_t size = 0;
    infile.read(reinterpret_cast<char*>(&size), sizeof(size_t));

    if (size == 0) {
        return {};
    }

    
    std::vector<std::complex<double>> coeffs(size);

    
    infile.read(reinterpret_cast<char*>(coeffs.data()), size * sizeof(std::complex<double>));

    if (!infile) {
        throw std::runtime_error("读取文件时发生错误或文件格式不正确: " + filename);
    }

    infile.close();
    std::cout << "[IO] 已从 " << filename << " 读取 " << size << " 个系数" << std::endl;
    return coeffs;
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
void FuncBootstrapExample(std::function<double(double)> f, int bits, int p);
void CarrymulExample(int bits, int p, int k);

Ciphertext<DCRTPoly> Conjugate(ConstCiphertext<DCRTPoly> ciphertext,
                               const std::map<usint, EvalKey<DCRTPoly>>& evalKeyMap) {
    const std::vector<DCRTPoly>& cv = ciphertext->GetElements();
    usint N                         = cv[0].GetRingDimension();

    std::vector<usint> vec(N);
    PrecomputeAutoMap(N, 2 * N - 1, &vec);

    auto algo = ciphertext->GetCryptoContext()->GetScheme();

    Ciphertext<DCRTPoly> result = ciphertext->Clone();

    algo->KeySwitchInPlace(result, evalKeyMap.at(2 * N - 1));

    std::vector<DCRTPoly>& rcv = result->GetElements();

    rcv[0] = rcv[0].AutomorphismTransform(2 * N - 1, vec);
    rcv[1] = rcv[1].AutomorphismTransform(2 * N - 1, vec);

    return result;
}

void FuncBootstrapExample(std::function<double(double)> f, int bits, int p) {
    std::cout << "---------------------EvalFBT Test--------------------" << std::endl << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    // For secure computations, use parameters.SetSecurityLevel(HEStd_128_classic);
    // and remove the two lines below.
    parameters.SetSecurityLevel(HEStd_128_classic);
    // parameters.SetRingDim(1 << 11);

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

    usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 4 + bits;
    // usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 4;
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
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, pow(2, 15) - 1);
    for (int i = 0; i < std::log(p); i++) {
        indexList.push_back(message(engine));
    }

    cc->EvalRotateKeyGen(keyPair.secretKey, indexList);
    std::vector<std::complex<double>> x;
    for (int i = 0; i < numSlots; ++i) {
        x.push_back(std::complex<double>((i) % 32768, (i + 1) % 32768));
    }
    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(x, 1, depth - levelBudget[1], nullptr, numSlots);
    ptxt->SetLength(numSlots);
    // std::cout << "Input = " << x << std::endl;

    Ciphertext<DCRTPoly> ctxt = cc->Encrypt(keyPair.publicKey, ptxt);

    // Order of the hermite interpolation, should be between 1 and 3
    int hermite_order = 1;

    // We evaluate f(ctxt) while bootstrapping
    std::chrono::system_clock::time_point start, end;
    Ciphertext<DCRTPoly> ctxtResult;
    auto num_test = 1;
    start         = std::chrono::system_clock::now();
    for (int i = 0; i < num_test; i++)
        ctxtResult = cc->EvalFuncBootstrap(ctxt, f, p, hermite_order);
    end = std::chrono::system_clock::now();
    std::cout << "level after fbt:" << depth - ctxtResult->GetLevel() << std::endl;
    std::cout << "noise deg ctxtResult:" << ctxtResult->GetNoiseScaleDeg() << std::endl;
    auto poly = ctxtResult->GetElements()[0];
    std::cout << "ctxtResult Modulus:" << std::endl;
    print_moduli_chain(poly);

    double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "EvalFuncBootstrap: " << totaltime / num_test << " ms" << std::endl;
    std::cout << "EvalFuncBootstrap amortize: " << totaltime / (num_test * ringDim) << " ms" << std::endl;
    Plaintext result;
    cc->Decrypt(keyPair.secretKey, ctxtResult, &result);
    result->SetLength(numSlots);
    std::vector<std::complex<double>> finalResult = result->GetCKKSPackedValue();
    int error_num                                 = 0;
    for (int i = 0; i < numSlots; i++) {
        // std::cout << "x.real(): " << x[i].real() << ", int(x[i]).real() % 16: " << int(x[i].real()) % 16
        //           << ", finalResult[i].real(): " << finalResult[i].real()
        //           << ", int(finalResult[i].real()) % 16: " << int(finalResult[i].real()) % 16 << std::endl;
        if (int(x[i].real()) % 16 != int(finalResult[i].real())) {
            error_num++;
            std::cout << "x.real(): " << x[i].real() << ", int(x[i]).real() % 16: " << int(x[i].real()) % 16
                      << ", finalResult[i].real(): " << finalResult[i].real()
                      << ", int(finalResult[i].real()) % 16: " << int(finalResult[i].real()) % 16 << std::endl;
        }
        // if (int(x[i].imag()) % 16 != int(finalResult[i].imag())) {
        //     error_num++;
        //     std::cout << "x.imag(): " << x[i].imag() << ", int(x[i]).imag() % 16: " << int(x[i].imag()) % 16
        //               << ", finalResult[i].imag(): " << finalResult[i].imag()
        //               << ", int(finalResult[i].imag()) % 16: " << int(finalResult[i].imag()) % 16 << std::endl;
        // }
        // std::cout << "x: " << x[i].real() << ", int(x[i]) % 16: " << int(x[i].real()) % 16
        //           << ", finalResult[i].real(): " << finalResult[i].real()
        //           << ", int(finalResult[i].real()) % 16: " << int(finalResult[i].real()) % 16 << std::endl;
    }
    std::cout << "error_num:" << error_num << std::endl;
    std::cout << "level after fbt:" << depth - ctxtResult->GetLevel() << std::endl << std::endl << std::endl;
}

void FuncBootstrap_kx_packed_Example(std::function<double(double)> f, int bits, int p) {
    std::cout << "---------------------EvalFBT Parallel Test (packed)--------------------" << std::endl << std::endl;
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
    usint dcrtBits               = 59;
    usint firstMod               = 60;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    std::vector<uint32_t> levelBudget = {4, 2};

    std::vector<uint32_t> bsgsDim = {0, 0};

    // usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 4 + bits + 2;
    usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 4 + 2;
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    // cc->Enable(FBTS);

    usint ringDim = cc->GetRingDimension();
    std::cout << "CKKS scheme is using ring dimension " << ringDim << " depth: " << depth << std::endl << std::endl;

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
    int hermite_order = 1;
    std::vector<std::complex<double>> hermite_coefficients;
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    switch (hermite_order) {
        case 1:
            hermite_coefficients = hermiteInterpOrder1Coeff([f](double x) -> double { return 0.5 * f(x); }, p);
            break;
        case 2:
            hermite_coefficients = hermiteInterpOrder2Coeff([f](double x) -> double { return 0.5 * f(x); }, p);
            break;
        case 3:
            hermite_coefficients = hermiteInterpOrder3Coeff([f](double x) -> double { return 0.5 * f(x); }, p);
            break;
        default:
            OPENFHE_THROW("Order of Hermite interpolation must be between 1 and 3");
            break;
    }
    end                = std::chrono::system_clock::now();
    int hermite_degree = hermite_coefficients.size();

    auto num_input   = numSlots / hermite_degree;
    double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "p: " << p << std::endl;
    std::cout << "num_input: " << num_input << std::endl;
    std::cout << "hermite_degree: " << hermite_degree << std::endl;
    std::vector<std::complex<double>> input;
    for (int i = 0; i < num_input; i++) {
        input.push_back(std::complex<double>(i % p, 0));
    }

    std::vector<std::complex<double>> x;
    for (int i = 0; i < num_input; i++) {
        for (int j = 0; j < hermite_degree; j++) {
            x.push_back((std::complex<double>(input[i].real() * j)));
        }
    }
    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(x, 1, depth - levelBudget[1], nullptr, numSlots);
    ptxt->SetLength(numSlots);
    Ciphertext<DCRTPoly> ctxt = cc->Encrypt(keyPair.publicKey, ptxt);

    start           = std::chrono::system_clock::now();
    auto ctxtResult = cc->EvalFuncBootstrap_kx(ctxt, p, hermite_order, hermite_coefficients);
    // int K           = 16;
    // int powR        = pow(2, 4);
    // auto f          = [K, powR](double x) -> std::complex<double> {
    //     return std::exp(std::complex<double>(0, 2 * M_PI * K * x / powR));
    // };
    end       = std::chrono::system_clock::now();
    totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "EvalFuncBootstrap: " << totaltime << " ms" << std::endl;
    Plaintext result;
    cc->Decrypt(keyPair.secretKey, ctxtResult, &result);
    result->SetLength(numSlots);
    std::vector<std::complex<double>> finalResult = result->GetCKKSPackedValue();
    // int error_num                                 = 0;
    for (int i = hermite_degree; i < 3 * hermite_degree; i++) {
        std::cout << "input[i]: " << x[i].real() << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
    }
}

void FuncBootstrap_kx_parallel_Example(std::function<double(double)> f, int bits, int p, int opt_depth) {
    std::cout << "---------------------EvalFBT Parallel Test--------------------" << std::endl << std::endl;
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
    Ciphertext<DCRTPoly> ctxt  = cc->Encrypt(keyPair.publicKey, ptxt);
    Ciphertext<DCRTPoly> ctxt1 = cc->Encrypt(keyPair.publicKey, ptxt);
    DCRTPoly& poly             = ctxt->GetElements()[0];
    std::cout << "ctxt Modulus:" << std::endl;
    print_moduli_chain(poly);
    // Order of the hermite interpolation, should be between 1 and 3
    int hermite_order = 1;
    std::vector<std::complex<double>> hermite_coefficients;
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    switch (hermite_order) {
        case 1:
            hermite_coefficients = hermiteInterpOrder1Coeff([f](double x) -> double { return 0.5 * f(x); }, p);
            break;
        case 2:
            hermite_coefficients = hermiteInterpOrder2Coeff([f](double x) -> double { return 0.5 * f(x); }, p);
            break;
        case 3:
            hermite_coefficients = hermiteInterpOrder3Coeff([f](double x) -> double { return 0.5 * f(x); }, p);
            break;
        default:
            OPENFHE_THROW("Order of Hermite interpolation must be between 1 and 3");
            break;
    }
    end              = std::chrono::system_clock::now();
    double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "hermite_coefficients time: " << totaltime << " ms" << std::endl;

    std::cout << "hermite_coefficients size " << hermite_coefficients.size() << std::endl;

    int zero_coe = 0;
    for (int i = 0; i < int(hermite_coefficients.size()); i++) {
        if (std::abs(hermite_coefficients[i].real()) < 1e-6 && std::abs(hermite_coefficients[i].imag()) < 1e-6) {
            zero_coe++;
        }
        if (std::abs(hermite_coefficients[i].real()) > 1e-6 || std::abs(hermite_coefficients[i].imag()) > 1e-6) {
            // std::cout << "coe not zero in " << i << std::endl;
        }
        // std::cout << hermite_coefficients[i] << ", " << std::endl;
    }
    std::cout << "zero_coe: " << zero_coe << std::endl;
    std::cout << "opt_depth " << opt_depth << std::endl << std::endl;
    // We evaluate f(ctxt) while bootstrapping
    auto use_img = true;
    std::cout << "use_img: " << use_img << std::endl;
    Ciphertext<DCRTPoly> ctxtResult;
    start = std::chrono::system_clock::now();
    cc->EvalMultInPlace(ctxt, 1);
    cc->ModReduceInPlace(ctxt);
    auto num_test = 1;
    for (int i = 0; i < num_test; i++)
        ctxtResult =
            cc->EvalFuncBootstrap_kx_parallel_opt(ctxt, p, hermite_order, hermite_coefficients, opt_depth, use_img);

    end       = std::chrono::system_clock::now();
    totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "level after fbt:" << depth - ctxtResult->GetLevel() << std::endl;
    // poly = ctxtResult->GetElements()[0];
    // print_moduli_chain(poly);
    std::cout << "EvalFuncBootstrap: " << totaltime / num_test << " ms" << std::endl;
    std::cout << "EvalFuncBootstrap amortize: " << totaltime / (num_test * ringDim) << " ms" << std::endl;
    // cc->ModReduceInPlace(ctxtResult);

    // ctxtResult      = cc->ModReduce(ctxtResult);
    // auto res        = cc->EvalMult(ctxtResult, 1);
    // DCRTPoly& poly1 = res->GetElements()[0];

    // std::cout << "res Modulus:" << std::endl;
    // print_moduli_chain(poly1);
    // cc->ModReduceInPlace(res);
    std::cout << "noise deg after fbt:" << ctxtResult->GetNoiseScaleDeg() << std::endl;
    // cc->ModReduceInPlace(res);
    // cc->RescaleInPlace(res);
    Plaintext result;
    cc->Decrypt(keyPair.secretKey, ctxtResult, &result);
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
        if (f(x[i].real()) != int(finalResult[i].real())) {
            error_num++;
            std::cout << "x.real(): " << x[i].real() << ",f(x.real): " << f(x[i].real())
                      << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
        }
        if (f(x[i].imag()) != int(finalResult[i].imag())) {
            error_num++;
            std::cout << "x.imag(): " << x[i].imag() << ",f(x.imag): " << f(x[i].imag())
                      << ", finalResult[i].imag(): " << finalResult[i].imag() << std::endl;
        }
        // error_num++;
        // std::cout << "i: " << i << ", x: " << x[i].real() << ", int(x[i]) % 16: " << int(x[i].real()) % 16
        //           << ", finalResult[i].real(): " << finalResult[i].real()
        //           << ", int(finalResult[i].real()) % 16: " << int(finalResult[i].real()) % 16 << std::endl;

        // }
    }
    std::cout << "error_num:" << error_num << std::endl;
}

void MVFuncBootstrap_kx_parallel_Example(int bits, int p, int opt_depth) {
    std::cout << "---------------------EvalMVFBT Parallel Test--------------------" << std::endl << std::endl;
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
    for (int i = -numSlots / 2; i < numSlots / 2; ++i) {
        x.push_back(std::complex<double>((i) % (p / 2), (i + 1) % (p / 2)));
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
        return ((int)x % 16);
    };
    auto f_2 = [](double x) -> double {
        return ((((int)x % 256) - ((int)x % 16)) / 16);
    };
    auto f_3 = [](double x) -> double {
        return ((((int)x % 4096) - ((int)x % 256)) / 256);
    };
    auto f_4 = [](double x) -> double {
        return (((int)x - ((int)x % 4096)) / 4096);
    };
    auto f_rs_4bit = [](double x) -> double {
        return (((int)x - ((int)x % 16)) / 16);
    };

    // auto f_1 = [](double x) -> double {
    //     return (double)((int)x == 16);
    // };
    // auto f_2 = [](double x) -> double {
    //     return (double)((int)x == 32);
    // };
    // auto f_3 = [](double x) -> double {
    //     return (double)((int)x == 64);
    // };
    // auto f_4 = [](double x) -> double {
    //     return (double)((int)x == 128);
    // };
    // auto f_rs_4bit = [](double x) -> double {
    //     return (double)(((int)x - ((int)x % 16)) / 16);
    // };
    // coe_1       = hermiteInterpOrder1Coeffshift([f_1](double x) -> double { return 0.5 * f_1(x); }, p);
    // coe_2       = hermiteInterpOrder1Coeffshift([f_2](double x) -> double { return 0.5 * f_2(x); }, p);
    // coe_3       = hermiteInterpOrder1Coeffshift([f_3](double x) -> double { return 0.5 * f_3(x); }, p);
    // coe_4       = hermiteInterpOrder1Coeffshift([f_4](double x) -> double { return 0.5 * f_4(x); }, p);
    // coe_rs_4bit = hermiteInterpOrder1Coeffshift([f_rs_4bit](double x) -> double { return 0.5 * f_rs_4bit(x); }, p);
    coe_1       = hermiteInterpOrder1Coeff([f_1](double x) -> double { return 0.5 * f_1(x); }, p);
    coe_2       = hermiteInterpOrder1Coeff([f_2](double x) -> double { return 0.5 * f_2(x); }, p);
    coe_3       = hermiteInterpOrder1Coeff([f_3](double x) -> double { return 0.5 * f_3(x); }, p);
    coe_4       = hermiteInterpOrder1Coeff([f_4](double x) -> double { return 0.5 * f_4(x); }, p);
    coe_rs_4bit = hermiteInterpOrder1Coeff([f_rs_4bit](double x) -> double { return 0.5 * f_rs_4bit(x); }, p);

    int zero_coe = 0;
    for (int i = 0; i < int(coe_1.size()); i++) {
        if (std::abs(coe_1[i].real()) < 1e-6 && std::abs(coe_1[i].imag()) < 1e-6) {
            zero_coe++;
        }
        if (std::abs(coe_1[i].real()) > 1e-6 || std::abs(coe_1[i].imag()) > 1e-6) {
            std::cout << "coe_1 not zero in " << i << std::endl;
        }
    }
    std::cout << "coe1 zero_coe: " << zero_coe << std::endl;

    zero_coe = 0;
    for (int i = 0; i < int(coe_2.size()); i++) {
        if (std::abs(coe_2[i].real()) < 1e-6 && std::abs(coe_2[i].imag()) < 1e-6) {
            zero_coe++;
        }
        if (std::abs(coe_2[i].real()) > 1e-6 || std::abs(coe_2[i].imag()) > 1e-6) {
            std::cout << "coe_2 not zero in " << i << std::endl;
        }
    }
    std::cout << "coe_2 zero_coe: " << zero_coe << std::endl;

    zero_coe = 0;
    for (int i = 0; i < int(coe_3.size()); i++) {
        if (std::abs(coe_3[i].real()) < 1e-6 && std::abs(coe_3[i].imag()) < 1e-6) {
            zero_coe++;
        }
        // if (std::abs(coe_3[i].real()) > 1e-6 || std::abs(coe_3[i].imag()) > 1e-6) {
        //     std::cout << "coe_3 not zero in " << i << std::endl;
        // }
    }
    std::cout << "coe_3 zero_coe: " << zero_coe << std::endl;

    zero_coe = 0;
    for (int i = 0; i < int(coe_4.size()); i++) {
        if (std::abs(coe_4[i].real()) < 1e-6 && std::abs(coe_4[i].imag()) < 1e-6) {
            zero_coe++;
        }
        // if (std::abs(coe_4[i].real()) > 1e-6 || std::abs(coe_4[i].imag()) > 1e-6) {
        //     std::cout << "coe_4 not zero in " << i << std::endl;
        // }
    }
    std::cout << "coe_4 zero_coe: " << zero_coe << std::endl;

    zero_coe = 0;
    for (int i = 0; i < int(coe_rs_4bit.size()); i++) {
        if (std::abs(coe_rs_4bit[i].real()) < 1e-6 && std::abs(coe_rs_4bit[i].imag()) < 1e-6) {
            zero_coe++;
        }
    }
    std::cout << "coe_rs_4bit zero_coe: " << zero_coe << std::endl;

    std::vector<std::vector<std::complex<double>>> coe_1234 = {coe_1, coe_2, coe_3, coe_4};
    std::vector<std::vector<std::complex<double>>> coe_1_rs = {coe_1, coe_rs_4bit};
    end                                                     = std::chrono::system_clock::now();
    double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

    // We evaluate f(ctxt) while bootstrapping
    auto use_img = true;
    std::cout << "use_img: " << use_img << std::endl;
    std::vector<Ciphertext<DCRTPoly>> ctxtResult;
    start = std::chrono::system_clock::now();
    cc->EvalMultInPlace(ctxt, 1);
    cc->ModReduceInPlace(ctxt);
    auto num_test = 5;
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
    int error_num                                 = 0;
    for (int i = 0; i < numSlots; i++) {
        // std::cout << "x.real(): " << x[i].real() << ",f_1(x.real): " << f_1(x[i].real())
        //           << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
        if (f_1(x[i].real()) != int(finalResult[i].real())) {
            error_num++;
            std::cout << "x.real(): " << x[i].real() << ",f_1(x.real): " << f_1(x[i].real())
                      << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
        }
        if (f_1(x[i].imag()) != int(finalResult[i].imag())) {
            error_num++;
            std::cout << "x.imag(): " << x[i].imag() << ",f_1(x.imag): " << f_1(x[i].imag())
                      << ", finalResult[i].imag(): " << finalResult[i].imag() << std::endl;
        }
    }
    std::cout << "func1 error_num:" << error_num << std::endl;

    cc->Decrypt(keyPair.secretKey, ctxtResult[1], &result);
    result->SetLength(numSlots);
    finalResult = result->GetCKKSPackedValue();
    error_num   = 0;
    for (int i = 0; i < numSlots; i++) {
        // std::cout << "x.real(): " << x[i].real() << ",f_2(x.real): " << f_2(x[i].real())
        //           << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
        if (f_2(x[i].real()) != int(finalResult[i].real())) {
            error_num++;
            std::cout << "x.real(): " << x[i].real() << ",f_2(x.real): " << f_2(x[i].real())
                      << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
        }
        if (f_2(x[i].imag()) != int(finalResult[i].imag())) {
            error_num++;
            std::cout << "x.imag(): " << x[i].imag() << ",f_2(x.imag): " << f_2(x[i].imag())
                      << ", finalResult[i].imag(): " << finalResult[i].imag() << std::endl;
        }
    }
    std::cout << "func2 error_num:" << error_num << std::endl;

    cc->Decrypt(keyPair.secretKey, ctxtResult[2], &result);
    result->SetLength(numSlots);
    finalResult = result->GetCKKSPackedValue();
    error_num   = 0;
    for (int i = 0; i < numSlots; i++) {
        // std::cout << "x.real(): " << x[i].real() << ",f_2(x.real): " << f_2(x[i].real())
        //           << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
        if (f_3(x[i].real()) != int(finalResult[i].real())) {
            error_num++;
            std::cout << "x.real(): " << x[i].real() << ",f_3(x.real): " << f_3(x[i].real())
                      << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
        }
        if (f_3(x[i].imag()) != int(finalResult[i].imag())) {
            error_num++;
            std::cout << "x.imag(): " << x[i].imag() << ",f_3(x.imag): " << f_3(x[i].imag())
                      << ", finalResult[i].imag(): " << finalResult[i].imag() << std::endl;
        }
    }
    std::cout << "func3 error_num:" << error_num << std::endl;

    cc->Decrypt(keyPair.secretKey, ctxtResult[3], &result);
    result->SetLength(numSlots);
    finalResult = result->GetCKKSPackedValue();
    error_num   = 0;
    for (int i = 0; i < numSlots; i++) {
        // std::cout << "x.real(): " << x[i].real() << ",f_2(x.real): " << f_2(x[i].real())
        //           << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
        if (f_4(x[i].real()) != int(finalResult[i].real())) {
            error_num++;
            std::cout << "x.real(): " << x[i].real() << ",f_4(x.real): " << f_4(x[i].real())
                      << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
        }
        if (f_4(x[i].imag()) != int(finalResult[i].imag())) {
            error_num++;
            std::cout << "x.imag(): " << x[i].imag() << ",f_4(x.imag): " << f_4(x[i].imag())
                      << ", finalResult[i].imag(): " << finalResult[i].imag() << std::endl;
        }
    }
    std::cout << "func3 error_num:" << error_num << std::endl;
}

void FuncBootstrap_kx_sparse_Example(std::function<double(double)> f, int bits, int p, int opt_depth, int mod) {
    std::cout << "---------------------EvalFBT Parallel Sparse Test--------------------" << std::endl << std::endl;
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

    std::vector<uint32_t> levelBudget = {3, 3};

    std::vector<uint32_t> bsgsDim = {0, 0};

    // usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 4 + bits + 2;
    usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 2 + opt_depth + 3 + 4;
    // usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 2 + opt_depth + 3 + 4 - 11;
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    // cc->Enable(FBTS);

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
    int error_precision = 14;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, p - 1);
    std::uniform_real_distribution<double> message_e(-1.0 / pow(2, error_precision),
                                                     -1.0 / (pow(2, error_precision) + 1));
    cc->EvalRotateKeyGen(keyPair.secretKey, indexList);
    // Input vector [0, 1, ..., p-2, p-1]
    // std::complex<double> input = 1;
    std::vector<std::complex<double>> x;
    for (int i = 0; i < numSlots; ++i) {
        x.push_back(std::complex<double>(message(engine), message(engine)));

        // x.push_back(std::complex<double>(double(message(engine) + message_e(engine)),
        //  double(message(engine) + message_e(engine))));
    }

    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(x, 1, depth - levelBudget[1] - 1, nullptr, numSlots);
    ptxt->SetLength(numSlots);
    // std::cout << "Input = " << x << std::endl;
    // Plaintext ptxt1            = cc->MakeCKKSPackedPlaintext(1.0, 1, depth - levelBudget[1] - 1, nullptr, numSlots);
    Ciphertext<DCRTPoly> ctxt  = cc->Encrypt(keyPair.publicKey, ptxt);
    Ciphertext<DCRTPoly> ctxt1 = cc->Encrypt(keyPair.publicKey, ptxt);
    DCRTPoly& poly             = ctxt->GetElements()[0];
    std::cout << "ctxt Modulus:" << std::endl;
    print_moduli_chain(poly);
    // Order of the hermite interpolation, should be between 1 and 3
    int hermite_order = 1;
    std::chrono::system_clock::time_point start, end;
    std::vector<std::complex<double>> hermite_coefficients;
    std::string filename = "testdata/coeffcient/o1/" + std::to_string(bits) + "bit_mod_" + std::to_string(mod) + ".bin";
    std::ifstream check_file(filename);
    if (check_file.good()) {
        check_file.close();
        std::cout << "检测到缓存文件，正在加载..." << std::endl;
        try {
            hermite_coefficients = LoadCoefficients(filename);
        }
        catch (const std::exception& e) {
            std::cerr << "加载失败: " << e.what() << "，将重新计算。" << std::endl;
            
            hermite_coefficients = hermiteInterpOrder1Coeff([f](double x) -> double { return 0.5 * f(x); }, p);
        }
    }
    else {
        
        std::cout << "未检测到缓存文件，正在计算..." << std::endl;
        check_file.close();
        hermite_coefficients = hermiteInterpOrder1Coeff([f](double x) -> double { return 0.5 * f(x); }, p);
        SaveCoefficients(filename, hermite_coefficients);
    }
    // hermite_coefficients = hermiteInterpOrder3Coeff([f](double x) -> double { return 0.5 * f(x); }, p);
    std::cout << "hermite_coefficients size " << hermite_coefficients.size() << std::endl;

    int zero_coe = 0;
    for (int i = 0; i < int(hermite_coefficients.size()); i++) {
        if (std::abs(hermite_coefficients[i].real()) < 1e-6 && std::abs(hermite_coefficients[i].imag()) < 1e-6) {
            zero_coe++;
        }
        if (std::abs(hermite_coefficients[i].real()) > 1e-6 || std::abs(hermite_coefficients[i].imag()) > 1e-6) {
            // std::cout << "coe not zero in " << i << std::endl;
        }
        // std::cout << hermite_coefficients[i] << ", " << std::endl;
    }
    std::cout << "zero_coe: " << zero_coe << std::endl;
    // We evaluate f(ctxt) while bootstrapping
    auto use_img = false;
    std::cout << "use_img: " << use_img << std::endl;
    Ciphertext<DCRTPoly> ctxtResult;
    start = std::chrono::system_clock::now();
    cc->EvalMultInPlace(ctxt, 1);
    cc->ModReduceInPlace(ctxt);
    auto num_test = 1;
    for (int i = 0; i < num_test; i++)
        ctxtResult = cc->EvalFuncBootstrap_kx_parallel_opt_sparse(ctxt, p, hermite_order, hermite_coefficients,
                                                                  opt_depth, use_img);

    end              = std::chrono::system_clock::now();
    double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "level after fbt:" << depth - ctxtResult->GetLevel() << std::endl;

    std::cout << "EvalFuncBootstrap: " << totaltime / num_test << " ms" << std::endl;
    std::cout << "EvalFuncBootstrap amortize: " << totaltime / (num_test * ringDim) << " ms" << std::endl;

    std::cout << "result level: " << ctxtResult->GetLevel() << std::endl;
    Plaintext result;
    cc->Decrypt(keyPair.secretKey, ctxtResult, &result);
    result->SetLength(numSlots);
    std::vector<std::complex<double>> finalResult = result->GetCKKSPackedValue();
    int error_num                                 = 0;
    double max_error                              = 0;
    double mean_error                             = 0;
    for (int i = 0; i < 16; i++) {
        // std::cout << "x.real(): " << x[i].real() << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
        if (fabs(f(x[i].real()) - finalResult[i].real()) > 0.1) {
            error_num++;
            std::cout << "x.real(): " << x[i].real() << ",f(x.real): " << f(x[i].real())
                      << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
        }
        // if (fabs(f(x[i].imag()) - finalResult[i].imag()) > 0.1) {
        //     error_num++;
        //     std::cout << "x.imag(): " << x[i].imag() << ",f(x.imag): " << f(x[i].imag())
        //               << ", finalResult[i].imag(): " << finalResult[i].imag() << std::endl;
        // }
        max_error = std::fabs(f(x[i].real()) - finalResult[i].real()) > max_error ?
                        std::fabs(f(x[i].real()) - finalResult[i].real()) :
                        max_error;
        // max_error = std::fabs(f(x[i].imag()) - finalResult[i].imag()) > max_error ?
        //                 std::fabs(f(x[i].imag()) - finalResult[i].imag()) :
        // max_error;
        mean_error += std::fabs(f(x[i].real()) - finalResult[i].real());
        // mean_error += std::fabs(f(x[i].imag()) - finalResult[i].imag());
    }
    std::cout << "max_error: 2^" << std::log2(max_error) << std::endl;
    std::cout << "mean_error: 2^" << std::log2(mean_error / numSlots) << std::endl;

    std::cout << "error_num:" << error_num << std::endl;
}

void Eval_mod_test(std::function<double(double)> f, int bits, int p, int opt_depth) {
    std::cout << "---------------------Eval mod Test--------------------" << std::endl << std::endl;
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

    std::vector<uint32_t> levelBudget = {3, 3};

    std::vector<uint32_t> bsgsDim = {0, 0};

    // usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 4 + bits + 2;
    usint depth = levelBudget[0] + levelBudget[1] + 2 + bits + 10;
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    // cc->Enable(FBTS);

    int ringDim = cc->GetRingDimension();
    std::cout << "CKKS scheme is using ring dimension " << ringDim << " depth: " << depth << std::endl << std::endl;

    // Precomputations for bootstrapping
    int numSlots = ringDim / 2;
    // cc->EvalFuncBootstrapSetup(levelBudget, bsgsDim, numSlots, bits);

    // Key Generation
    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    // cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);
    std::vector<int32_t> indexList;
    for (int i = 0; i <= std::log2(ringDim) + 1; i++) {
        indexList.push_back(pow(2, i));
        // std::cout << "rot key for " << pow(2, i) << std::endl;
    }
    indexList.push_back(2 * ringDim - 1);
    cc->EvalRotateKeyGen(keyPair.secretKey, indexList);
    std::vector<std::complex<double>> x;
    // random input form (0,p)
    int error_precision = 20;
    std::cout << "error_precision: " << error_precision << std::endl;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, p - 1);
    std::uniform_real_distribution<double> message_e(-1.0 / pow(2, error_precision),
                                                     -1.0 / (pow(2, error_precision) + 1));

    std::vector<size_t> nonzero_indices;
    const double threshold = 1e-9;
    // auto hermite_coefficients = hermiteInterpOrder1Coeff([f](double x) -> double { return 0.5 * f(x); }, p);
    // std::string filename = "testdata/coeffcient/" + std::to_string(bits) + "bit_mod_32.bin";
    // std::string filename = "testdata/coeffcient/o1/" + std::to_string(bits) + "bit_mod_2.bin";

    std::string filename = "testdata/coeffcient/o1/" + std::to_string(bits) + "bit_mod_16.bin";

    std::vector<std::complex<double>> hermite_coefficients;

    std::ifstream check_file(filename);

    if (check_file.good()) {
        check_file.close();
        std::cout << "检测到缓存文件，正在加载..." << std::endl;
        try {
            hermite_coefficients = LoadCoefficients(filename);
        }
        catch (const std::exception& e) {
            std::cerr << "加载失败: " << e.what() << "，将重新计算。" << std::endl;
            
            hermite_coefficients = hermiteInterpOrder1Coeff([f](double x) -> double { return 0.5 * f(x); }, p);
        }
    }
    else {
        
        std::cout << "未检测到缓存文件，正在计算..." << std::endl;
        check_file.close();
        hermite_coefficients = hermiteInterpOrder1Coeff([f](double x) -> double { return 0.5 * f(x); }, p);
        SaveCoefficients(filename, hermite_coefficients);
    }
    // hermite_coefficients = hermiteInterpOrder3Coeff([f](double x) -> double { return 0.5 * f(x); }, p);
    std::cout << "hermite_coefficients size " << hermite_coefficients.size() << std::endl;
    int zero_coe = 0;
    for (int i = 0; i < int(hermite_coefficients.size()); i++) {
        if (std::abs(hermite_coefficients[i].real()) < 1e-6 && std::abs(hermite_coefficients[i].imag()) < 1e-6) {
            zero_coe++;
        }
        if (std::abs(hermite_coefficients[i].real()) > 1e-6 || std::abs(hermite_coefficients[i].imag()) > 1e-6) {
            // std::cout << "coe not zero in " << i << std::endl;
        }
        // std::cout << hermite_coefficients[i] << ", " << std::endl;
    }
    std::cout << "zero_coe: " << zero_coe << std::endl;
    for (size_t i = 0; i < hermite_coefficients.size(); ++i) {
        if (std::abs(hermite_coefficients[i]) > threshold) {
            nonzero_indices.push_back(i);
        }
    }
    int step = 0;
    for (size_t idx : nonzero_indices) {
        if (idx > 0) {
            step = idx;
            break;
        }
    }
    std::cout << "step: " << step << std::endl;
    for (int i = 0; i < numSlots; ++i) {
        // x.push_back(std::complex<double>(double(message(engine) + message_e(engine)), 0));

        x.push_back(std::complex<double>(double(message(engine)), 0));
    }
    // for (int i = 0; i < 16; ++i) {
    //     std::cout << "x[" << i << "]: " << std::setprecision(10) << x[i] << std::endl;
    // }
    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(x, 1, 0, nullptr, numSlots);
    ptxt->SetLength(numSlots);
    Ciphertext<DCRTPoly> ctxt = cc->Encrypt(keyPair.publicKey, ptxt);

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    ctxt = cc->EvalMult(ctxt, 1.0 / (p));
    cc->ModReduceInPlace(ctxt);
    // ctxt         = cc->GetScheme()->MultByInteger(ctxt, step);
    int num_test = 10;
    int K        = 1;
    int R_FUNC   = std::log2(step / K);
    int powR_val = pow(2, R_FUNC);
    auto f_exp   = [K, powR_val, step](double x) -> std::complex<double> {
        return std::exp(std::complex<double>(0, 2 * M_PI * x * K));
    };
    start = std::chrono::system_clock::now();
    Ciphertext<DCRTPoly> ctxtResult;
    for (int i = 0; i < num_test; i++) {
        uint32_t poly_degree = 16;
        auto ctxtExp         = cc->EvalChebyshevFunction(f_exp, ctxt, -1, 1, poly_degree);
        // auto ctxtExp = cc->EvalChebyshevSeries(ctxt, coeff_exp_16_double_46, -1, 1);

        // auto ctxtExp = cc->EvalChebyshevSeries(ctxt, coeff_exp_16_double_46, -1, 1);
        // std::cout << "EvalChebyshevFunction finish " << std::endl;

        // Plaintext expresult;
        // cc->Decrypt(keyPair.secretKey, ctxtExp, &expresult);
        // expresult->SetLength(numSlots);
        // std::vector<std::complex<double>> exp_Result = expresult->GetCKKSPackedValue();
        // for (int i = 0; i < 65; i++) {
        //     std::cout << "x.real(): " << x[i].real() << ",f_exp(x.real): " << f_exp(x[i].real() / (p)).real()
        //               << ", exp_Result[i].real(): " << exp_Result[i].real() << std::endl;
        // }
        for (int k = 0; k < R_FUNC; ++k) {
            cc->EvalSquareInPlace(ctxtExp);
            cc->ModReduceInPlace(ctxtExp);
        }
        // end              = std::chrono::system_clock::now();
        // double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
        // std::cout << " exp time: " << totaltime << " ms" << std::endl;
        // start = std::chrono::system_clock::now();
        // auto ctxtInterp = cc->EvalHermiteFunction([f](double x) -> double { return 0.5 * f(x); }, ctxtExp, p, 1);
        auto ctxtInterp = cc->EvalHermiteFunction_kx_sparse_PS(ctxtExp, hermite_coefficients, step);
        auto conj       = cc->GetScheme()->EvalConj(ctxtInterp);
        cc->EvalAddInPlace(ctxtInterp, conj);
        ctxtResult = ctxtInterp;
    }
    end              = std::chrono::system_clock::now();
    double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

    end = std::chrono::system_clock::now();

    Plaintext result;
    cc->Decrypt(keyPair.secretKey, ctxtResult, &result);
    result->SetLength(numSlots);
    std::vector<std::complex<double>> finalResult = result->GetCKKSPackedValue();
    int error_num                                 = 0;
    double max_error                              = 0;
    double mean_error                             = 0;
    for (int i = 0; i < 16; i++) {
        // std::cout << "x.real(): " << x[i].real() << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
        if (fabs(f(x[i].real()) - finalResult[i].real()) > 0.1) {
            error_num++;
            std::cout << "x.real(): " << x[i].real() << ",f(x.real): " << f(x[i].real())
                      << ", finalResult[i].real(): " << finalResult[i].real() << std::endl;
        }
        max_error = std::fabs(f(x[i].real()) - finalResult[i].real()) > max_error ?
                        std::fabs(f(x[i].real()) - finalResult[i].real()) :
                        max_error;
        mean_error += std::fabs(f(x[i].real()) - finalResult[i].real());
    }
    std::cout << "max_error: 2^" << std::log2(max_error) << std::endl;
    std::cout << "mean_error: 2^" << std::log2(mean_error / numSlots) << std::endl;
    std::cout << "error_num:" << error_num << std::endl;
    std::cout << " hermite time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms" << std::endl;
    std::cout << " totaltime time: " << totaltime / num_test << " ms" << std::endl;
    std::cout << " amortize time: " << totaltime / (numSlots * num_test) << " ms" << std::endl;
    std::cout << "level after fbt:" << depth - ctxtResult->GetLevel() << std::endl;
}

std::vector<std::complex<double>> hermiteInterpOrder1Coefftest(std::function<double(double)> f, int p) {
    std::vector<double> y;
    for (int i = 0; i < p; ++i)
        y.push_back(f(i));

    std::vector<std::complex<double>> alpha(p, 0.0);

    double sum_y = 0;
    for (int i = 0; i < p; ++i)
        sum_y += y[i];

    alpha[0] = sum_y / (double)p;

    for (int k = 1; k < p; ++k) {
        for (int l = 0; l < p; ++l) {
            std::complex<double> exp_term = exp(std::complex<double>(0, -2 * M_PI * k * l / p));
            alpha[k] += y[l] * exp_term;
        }
        alpha[k] *= 2.0 * (p - k) / (double)(p * p);
    }
    return alpha;
}

double EvalPoly_plain(double x, const std::vector<std::complex<double>>& alpha, int p) {
    std::complex<double> sum = 0;

    for (int k = 0; k < p; ++k) {
        std::complex<double> basis = exp(std::complex<double>(0, 2 * M_PI * k * x / p));
        sum += alpha[k] * basis;
    }

    return sum.real();  
}

int main(int argc, char* argv[]) {
    int bits = 15;
    int mod  = 16;
    int p    = pow(2, bits);
    auto f   = [mod](double x) -> double {
        return ((int)x % mod);
    };
    int mod_depth = 4;

    Eval_mod_test(f, bits, p, mod_depth);
    // f is the function that is evaluated during bootstrapping
    // auto f = [](double x) -> double {
    //     // return (double)((int)x % 16);
    //     return (double)(((int)x - ((int)x % 4096)) / 4096);
    // };
    // auto f = [](double x) -> double {
    //     return (int)x == 256;
    // };

    /*auto f = [](double x) -> double { return x * x; };*/
    // FuncBootstrapExample(f, bits, p);
    // FuncBootstrap_kx_packed_Example(f, bits, p);

    // int opt_depth = 4;
    // FuncBootstrap_kx_parallel_Example(f, bits, p, opt_depth);
    // MVFuncBootstrap_kx_parallel_Example(bits, p, opt_depth);

    // FuncBootstrap_kx_sparse_Example(f, bits, p, mod_depth, mod);

    // int k = 4;
    // CarryExample(bits, p, k);
    // CarrymulExample16bits1(bits, p, k);
    // CarrymulExample32bits(bits, p, k);
    // for (int bits = 6; bits < 18; bits++) {
    //     auto f = [](double x) -> double {
    //         return ((int)x % 2);
    //     };
    //     int p                = pow(2, bits);
    //     std::string filename = "testdata/coeffcient/o3/" + std::to_string(bits) + "bit_mod_2.bin";

    //     std::vector<std::complex<double>> hermite_coefficients;

    //     std::ifstream check_file(filename);

    //     if (check_file.good()) {
    //         check_file.close();
    
    //         try {
    //             hermite_coefficients = LoadCoefficients(filename);
    //         }
    //         catch (const std::exception& e) {
    
    
    //             hermite_coefficients = hermiteInterpOrder3Coeff([f](double x) -> double { return 0.5 * f(x); }, p);
    //         }
    //     }
    //     else {
    
    
    //         check_file.close();
    //         hermite_coefficients = hermiteInterpOrder3Coeff([f](double x) -> double { return 0.5 * f(x); }, p);
    //         SaveCoefficients(filename, hermite_coefficients);
    //     }
    // }
}