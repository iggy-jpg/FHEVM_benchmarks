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

#include <iostream>
#include <unordered_set>
#include <chrono>
using namespace lbcrypto;
Ciphertext<DCRTPoly> MultByInteger(const ConstCiphertext<DCRTPoly>& ciphertext, const int64_t constant) {
    auto result = ciphertext->Clone();
    for (auto& c : result->GetElements())
        c *= static_cast<typename DCRTPoly::Integer>(constant);
    return result;
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

// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
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

std::vector<std::vector<double>> read_vectors_from_csv(const std::string& filename, int N) {
    std::vector<std::vector<double>> sift;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return sift;
    }

    std::string line;
    int count = 0;

    while (std::getline(file, line) && count < N) {
        std::vector<double> vector;
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss, value, ',')) {
            vector.push_back(std::stod(value));
        }

        sift.push_back(vector);
        count++;
    }

    file.close();
    return sift;
}

int extractBits(double num, int start, int length) {
    return (int(num) >> start) & ((1 << length) - 1);
}

std::vector<double> splitNumber(double x, int w, int Bg) {
    if (x < 0 || x >= (1LL << (w * Bg))) {
        std::cerr << "Input " << x << " out of range for precision " << (1LL << (w * Bg)) << "\n";
    }
    std::vector<double> result;
    for (int i = 0; i < w; ++i) {
        int start_bit = i * Bg;
        int extracted = extractBits(x, start_bit, Bg);
        result.push_back(extracted);
    }
    return result;
}

std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix) {
    if (matrix.empty() || matrix[0].empty()) {
        return {};
    }

    int rows = matrix.size();
    int cols = matrix[0].size();

    std::vector<std::vector<double>> transposed(cols, std::vector<double>(rows));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            transposed[j][i] = matrix[i][j];
        }
    }

    return transposed;
}

std::vector<float> read_fvec(const std::string& filename, size_t index) {
    std::ifstream input(filename, std::ios::binary);
    if (!input.is_open()) {
        throw std::runtime_error("Error: cannot open fvecs file " + filename);
    }

    int dim;
    input.read((char*)&dim, 4);
    input.seekg(0, std::ios::end);
    size_t file_size = input.tellg();
    size_t num_vecs  = file_size / ((dim + 1) * 4);
    if (index >= num_vecs) {
        throw std::runtime_error("Index out of range in fvecs file");
    }
    input.seekg(index * (dim + 1) * 4, std::ios::beg);

    int d;
    input.read((char*)&d, 4);
    std::vector<float> vec(d);
    input.read((char*)vec.data(), d * 4);
    return vec;
}

std::vector<std::vector<double>> quantization(std::vector<double> vec, uint32_t precision, uint32_t Bg) {
    int w      = precision / Bg;
    int length = vec.size();
    std::vector<std::vector<double>> quant_plain(w, std::vector<double>(length));
    for (int i = 0; i < length; i++) {
        auto split_vi     = splitNumber(vec[i], w, Bg);
        quant_plain[0][i] = split_vi[0];
        quant_plain[1][i] = split_vi[1];
    }
    return quant_plain;
}

void compare_topk(const std::vector<int>& plain_dist, const std::vector<int>& res_dist, int K) {
    int n = plain_dist.size();
    if (K > n)
        K = n;
    std::vector<int> plain_idx(n), res_idx(n);
    for (int i = 0; i < n; i++) {
        plain_idx[i] = i;
        res_idx[i]   = i;
    }

    std::sort(plain_idx.begin(), plain_idx.end(), [&](int a, int b) { return plain_dist[a] > plain_dist[b]; });

    std::sort(res_idx.begin(), res_idx.end(), [&](int a, int b) { return res_dist[a] > res_dist[b]; });

    std::unordered_set<int> plain_topk(plain_idx.begin(), plain_idx.begin() + K);
    std::unordered_set<int> res_topk(res_idx.begin(), res_idx.begin() + K);

    int match = 0;
    for (int id : plain_topk) {
        if (res_topk.count(id)) {
            match++;
        }
    }

    // std::cout << "Top-" << K << " res: " << std::endl;
    // std::cout << "Plain Top-" << K << " indices: ";
    // for (int i = 0; i < K; i++)
    //     std::cout << plain_idx[i] << " ";
    // std::cout << std::endl;

    // std::cout << "Enc Top-" << K << " indices: ";
    // for (int i = 0; i < K; i++)
    //     std::cout << res_idx[i] << " ";
    // std::cout << std::endl;

    std::cout << "Top-" << K << " match num: " << match << " / " << K << std::endl;
}

std::vector<Ciphertext<lbcrypto::DCRTPoly>> HomomorphicCarryMult32(
    Ciphertext<lbcrypto::DCRTPoly>& ctxt0, Ciphertext<lbcrypto::DCRTPoly>& ctxt1, Ciphertext<lbcrypto::DCRTPoly>& ctxt2,
    Ciphertext<lbcrypto::DCRTPoly>& ctxt3, Ciphertext<lbcrypto::DCRTPoly>& ctxt_rshift_2k,
    std::function<double(double)> mod2k, size_t p, int k, int depth, std::vector<uint32_t> levelBudget,
    lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey, int numSlots) {
    auto cc           = ctxt0->GetCryptoContext();
    auto algo         = cc->GetScheme();
    int modLevel      = depth - levelBudget[1];
    int hermite_order = 1;
    std::cout << " hermite_order " << hermite_order << std::endl;
    Plaintext result0, result1, result2, result3, result4, result5, result6, result7;
    // Step 1: compute raw multiplications M0~M6 with Rescale after each EvalMult
    std::cout << " raw multiplications " << std::endl;
    std::chrono::system_clock::time_point start, end;
    start   = std::chrono::system_clock::now();
    auto M0 = cc->EvalMult(ctxt0, ctxt0);
    cc->RescaleInPlace(M0);
    end              = std::chrono::system_clock::now();
    double totaltime = 0;
    totaltime += double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << numSlots << " slots single mul time:  " << totaltime << std::endl;

    auto M1 = cc->EvalMult(ctxt0, ctxt1);
    cc->RescaleInPlace(M1);
    M1 = cc->EvalMult(M1, 2);
    cc->RescaleInPlace(M1);

    auto M2a = cc->EvalMult(ctxt1, ctxt1);
    cc->RescaleInPlace(M2a);
    auto M2b = cc->EvalMult(ctxt0, ctxt2);
    cc->RescaleInPlace(M2b);
    M2b = cc->EvalMult(M2b, 2);
    cc->RescaleInPlace(M2b);
    auto M2 = cc->EvalAdd(M2a, M2b);

    auto M3a = cc->EvalMult(ctxt0, ctxt3);
    cc->RescaleInPlace(M3a);
    M3a = cc->EvalMult(M3a, 2);
    cc->RescaleInPlace(M3a);
    auto M3b = cc->EvalMult(ctxt1, ctxt2);
    cc->RescaleInPlace(M3b);
    M3b = cc->EvalMult(M3b, 2);
    cc->RescaleInPlace(M3b);
    auto M3 = cc->EvalAdd(M3a, M3b);

    auto M4a = cc->EvalMult(ctxt2, ctxt2);
    cc->RescaleInPlace(M4a);
    auto M4b = cc->EvalMult(ctxt1, ctxt3);
    cc->RescaleInPlace(M4b);
    M4b = cc->EvalMult(M4b, 2);
    cc->RescaleInPlace(M4b);
    auto M4 = cc->EvalAdd(M4a, M4b);

    auto M5 = cc->EvalMult(ctxt2, ctxt3);
    cc->RescaleInPlace(M5);
    M5 = cc->EvalMult(M5, 2);
    cc->RescaleInPlace(M5);

    auto M6 = cc->EvalMult(ctxt3, ctxt3);
    cc->RescaleInPlace(M6);

    cc->Decrypt(privateKey, M0, &result0);
    cc->Decrypt(privateKey, M1, &result1);
    cc->Decrypt(privateKey, M2, &result2);
    cc->Decrypt(privateKey, M3, &result3);
    cc->Decrypt(privateKey, M4, &result4);
    cc->Decrypt(privateKey, M5, &result5);
    cc->Decrypt(privateKey, M6, &result6);
    // cc->Decrypt(keyPair.secretKey, cres[7], &result7);
    result0->SetLength(numSlots);
    result1->SetLength(numSlots);
    result2->SetLength(numSlots);
    result3->SetLength(numSlots);
    result4->SetLength(numSlots);
    result5->SetLength(numSlots);
    result6->SetLength(numSlots);
    // result7->SetLength(numSlots);
    auto finalResult0 = result0->GetCKKSPackedValue();
    auto finalResult1 = result1->GetCKKSPackedValue();
    auto finalResult2 = result2->GetCKKSPackedValue();
    auto finalResult3 = result3->GetCKKSPackedValue();
    auto finalResult4 = result4->GetCKKSPackedValue();
    auto finalResult5 = result5->GetCKKSPackedValue();
    auto finalResult6 = result6->GetCKKSPackedValue();
    // auto finalResult7 = result7->GetCKKSPackedValue();
    std::cout << "M0: " << finalResult0[0].real() << ", M1: " << finalResult1[0].real()
              << ", M2: " << finalResult2[0].real() << ", M3: " << finalResult3[0].real()
              << ", M4: " << finalResult4[0].real() << ", M5  " << finalResult5[0].real()
              << ", M6: " << finalResult6[0].real() << std::endl;
    auto f_mod16 = [](double x) -> double {
        // x - 16 * round(x/16)
        
        
        return x - 16.0 * std::floor((x + 8.0) / 16.0);
    };
    std::vector<std::complex<double>> coe_mod16 =
        hermiteInterpOrder1Coeffshift([f_mod16](double x) -> double { return 0.5 * f_mod16(x); }, p);
    // Step 2: EvalFuncBootstrap on each term to get mod2^k
    auto ModBoot = [&](Ciphertext<DCRTPoly> c) {
        auto boot = c->Clone();
        algo->LevelReduceInternalInPlace(boot, modLevel - boot->GetLevel());
        return cc->EvalFuncBootstrap_kx_parallel_opt_sparse(boot, p, 1, coe_mod16, 1, false);
    };
    // Step 3: Carry Propagation per level
    auto CarryStep = [&](Ciphertext<DCRTPoly> lower, Ciphertext<DCRTPoly> higher_mod) {
        auto high_part = cc->EvalSub(lower, higher_mod);
        algo->AdjustLevelsInPlace(high_part, ctxt_rshift_2k);
        auto high_div = cc->EvalMult(high_part, ctxt_rshift_2k);
        cc->RescaleInPlace(high_div);
        return high_div;
    };

    std::cout << "ModBoot" << std::endl;
    auto M0_mod = ModBoot(M0);
    auto M1_mod = ModBoot(M1);
    auto M2_mod = ModBoot(M2);
    auto M3_mod = ModBoot(M3);
    auto M4_mod = ModBoot(M4);
    auto M5_mod = ModBoot(M5);
    auto M6_mod = ModBoot(M6);

    cc->Decrypt(privateKey, M0_mod, &result0);
    cc->Decrypt(privateKey, M1_mod, &result1);
    cc->Decrypt(privateKey, M2_mod, &result2);
    cc->Decrypt(privateKey, M3_mod, &result3);
    cc->Decrypt(privateKey, M4_mod, &result4);
    cc->Decrypt(privateKey, M5_mod, &result5);
    cc->Decrypt(privateKey, M6_mod, &result6);
    // cc->Decrypt(keyPair.secretKey, cres[7], &result7);
    result0->SetLength(numSlots);
    result1->SetLength(numSlots);
    result2->SetLength(numSlots);
    result3->SetLength(numSlots);
    result4->SetLength(numSlots);
    result5->SetLength(numSlots);
    result6->SetLength(numSlots);
    // result7->SetLength(numSlots);
    finalResult0 = result0->GetCKKSPackedValue();
    finalResult1 = result1->GetCKKSPackedValue();
    finalResult2 = result2->GetCKKSPackedValue();
    finalResult3 = result3->GetCKKSPackedValue();
    finalResult4 = result4->GetCKKSPackedValue();
    finalResult5 = result5->GetCKKSPackedValue();
    finalResult6 = result6->GetCKKSPackedValue();
    // auto finalResult7 = result7->GetCKKSPackedValue();
    std::cout << "M0_mod: " << finalResult0[0].real() << ", M1_mod: " << finalResult1[0].real()
              << ", M2_mod: " << finalResult2[0].real() << ", M3_mod: " << finalResult3[0].real()
              << ", M4_mod: " << finalResult4[0].real() << ", M5_mod: " << finalResult5[0].real()
              << ", M6_mod: " << finalResult6[0].real() << std::endl;

    std::cout << "CarryStep" << std::endl;
    auto M1_high = CarryStep(M1, M1_mod);
    auto M2_high = CarryStep(M2, M2_mod);
    auto M3_high = CarryStep(M3, M3_mod);
    auto M4_high = CarryStep(M4, M4_mod);
    auto M5_high = CarryStep(M5, M5_mod);
    auto M6_high = CarryStep(M6, M6_mod);
    // Step 3: Carry Propagation per level (including Mi_high at each stage)
    auto t1 = cc->EvalAdd(CarryStep(M0, M0_mod), M1_mod);
    auto t2 = cc->EvalAdd(cc->EvalAdd(CarryStep(t1, ModBoot(t1)), M2_mod), M1_high);
    auto t3 = cc->EvalAdd(cc->EvalAdd(CarryStep(t2, ModBoot(t2)), M3_mod), M2_high);
    auto t4 = cc->EvalAdd(cc->EvalAdd(CarryStep(t3, ModBoot(t3)), M4_mod), M3_high);
    auto t5 = cc->EvalAdd(cc->EvalAdd(CarryStep(t4, ModBoot(t4)), M5_mod), M4_high);
    auto t6 = cc->EvalAdd(cc->EvalAdd(CarryStep(t5, ModBoot(t5)), M6_mod), M5_high);
    cc->Decrypt(privateKey, t1, &result0);
    cc->Decrypt(privateKey, t2, &result1);
    cc->Decrypt(privateKey, t3, &result2);
    cc->Decrypt(privateKey, t4, &result3);
    cc->Decrypt(privateKey, t5, &result4);
    cc->Decrypt(privateKey, t6, &result5);
    result0->SetLength(numSlots);
    result1->SetLength(numSlots);
    result2->SetLength(numSlots);
    result3->SetLength(numSlots);
    result4->SetLength(numSlots);
    result5->SetLength(numSlots);
    finalResult0 = result0->GetCKKSPackedValue();
    finalResult1 = result1->GetCKKSPackedValue();
    finalResult2 = result2->GetCKKSPackedValue();
    finalResult3 = result3->GetCKKSPackedValue();
    finalResult4 = result4->GetCKKSPackedValue();
    finalResult5 = result5->GetCKKSPackedValue();
    std::cout << "t1: " << finalResult0[0].real() << ", t2: " << finalResult1[0].real()
              << ", t3: " << finalResult2[0].real() << ", t4: " << finalResult3[0].real()
              << ", t5: " << finalResult4[0].real() << ", t6: " << finalResult5[0].real() << std::endl;
    auto mres7 = cc->EvalAdd(CarryStep(t6, ModBoot(t6)), M6_high);
    std::cout << "t7 level: " << mres7->GetLevel() << std::endl;
    std::cout << "CarryStep finish" << std::endl;
    // Step 4: Return all segments mres0 ~ mres7
    std::vector<Ciphertext<DCRTPoly>> mres = {M0_mod,      ModBoot(t1), ModBoot(t2), ModBoot(t3),
                                              ModBoot(t4), ModBoot(t5), ModBoot(t6), mres7};

    std::cout << "Carry finish" << std::endl;
    return mres;
}

std::vector<Ciphertext<lbcrypto::DCRTPoly>> HomomorphicCarryMult16(
    Ciphertext<lbcrypto::DCRTPoly>& ctxt0, Ciphertext<lbcrypto::DCRTPoly>& ctxt1,
    Ciphertext<lbcrypto::DCRTPoly>& ctxt_rshift_2k, std::function<double(double)> mod2k, size_t p, int k, int depth,
    std::vector<uint32_t> levelBudget, lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey, int numSlots) {
    auto cc           = ctxt0->GetCryptoContext();
    auto algo         = cc->GetScheme();
    int modLevel      = depth - levelBudget[1];
    int hermite_order = 1;
    std::cout << " hermite_order " << hermite_order << std::endl;
    Plaintext result0, result1, result2, result3, result4;
    // Step 1: compute raw multiplications M0~M6 with Rescale after each EvalMult
    std::cout << " raw multiplications " << std::endl;
    std::chrono::system_clock::time_point start, end;
    start   = std::chrono::system_clock::now();
    auto M0 = cc->EvalMult(ctxt0, ctxt0);
    cc->RescaleInPlace(M0);
    end              = std::chrono::system_clock::now();
    double totaltime = 0;
    totaltime += double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << numSlots << " slots single mul time:  " << totaltime << std::endl;

    auto M1 = cc->EvalMult(ctxt0, ctxt1);
    cc->RescaleInPlace(M1);
    M1 = cc->EvalMult(M1, 2);
    cc->RescaleInPlace(M1);

    auto M2 = cc->EvalMult(ctxt1, ctxt1);
    cc->RescaleInPlace(M2);

    // Step 2: EvalFuncBootstrap on each term to get mod2^k
    auto ModBoot = [&](Ciphertext<DCRTPoly> c) {
        auto boot = c->Clone();
        algo->LevelReduceInternalInPlace(boot, modLevel - boot->GetLevel());
        return cc->EvalFuncBootstrap(boot, mod2k, p, hermite_order);
    };
    // Step 3: Carry Propagation per level
    auto CarryStep = [&](Ciphertext<DCRTPoly> lower, Ciphertext<DCRTPoly> higher_mod) {
        auto high_part = cc->EvalSub(lower, higher_mod);
        // algo->AdjustLevelsInPlace(high_part, ctxt_rshift_2k);
        auto high_div = cc->EvalMult(high_part, 0.0625);
        cc->RescaleInPlace(high_div);
        return high_div;
    };

    std::cout << "ModBoot" << std::endl;
    auto M0_mod = ModBoot(M0);
    auto M1_mod = ModBoot(M1);
    auto M2_mod = ModBoot(M2);

    cc->Decrypt(privateKey, M0, &result0);
    cc->Decrypt(privateKey, M1, &result1);
    cc->Decrypt(privateKey, M2, &result2);
    result0->SetLength(numSlots);
    result1->SetLength(numSlots);
    result2->SetLength(numSlots);
    auto finalResult0 = result0->GetCKKSPackedValue();
    auto finalResult1 = result1->GetCKKSPackedValue();
    auto finalResult2 = result2->GetCKKSPackedValue();
    std::cout << "M0_mod: " << finalResult0[0].real() << ", M1_mod: " << finalResult1[0].real()
              << ", M2_mod: " << finalResult2[0].real() << std::endl;

    // std::cout << "CarryStep" << std::endl;
    std::vector<Ciphertext<DCRTPoly>> tmod;
    auto M1_high = CarryStep(M1, M1_mod);
    auto M2_high = CarryStep(M2, M2_mod);
    // Step 3: Carry Propagation per level (including Mi_high at each stage)
    auto t1 = cc->EvalAdd(CarryStep(M0, M0_mod), M1_mod);
    tmod.push_back(ModBoot(t1));
    auto t2 = cc->EvalAdd(cc->EvalAdd(CarryStep(t1, tmod[0]), M2_mod), M1_high);
    tmod.push_back(ModBoot(t2));
    auto mres3 = cc->EvalAdd(CarryStep(t2, tmod[1]), M2_high);

    // Step 4: Return all segments mres0 ~ mres7
    std::vector<Ciphertext<DCRTPoly>> mres = {M0_mod, tmod[0], tmod[1], mres3};
    cc->Decrypt(privateKey, M0_mod, &result0);
    cc->Decrypt(privateKey, t2, &result1);
    cc->Decrypt(privateKey, mres3, &result2);
    // cc->Decrypt(privateKey, t4, &result3);
    // cc->Decrypt(privateKey, t5, &result4);
    // cc->Decrypt(privateKey, t6, &result5);
    result0->SetLength(numSlots);
    result1->SetLength(numSlots);
    result2->SetLength(numSlots);
    finalResult0 = result0->GetCKKSPackedValue();
    finalResult1 = result1->GetCKKSPackedValue();
    finalResult2 = result2->GetCKKSPackedValue();
    std::cout << "t1: " << finalResult0[0].real() << ", t2: " << finalResult1[0].real() << std::endl;
    // std::cout << "t1: " << finalResult0[0].real() << ", t2: " << finalResult1[0].real()
    //           << ", t3: " << finalResult2[0].real() << ", t4: " << finalResult3[0].real()
    //           << ", t5: " << finalResult4[0].real() << ", t6: " << finalResult5[0].real() << std::endl;
    // std::cout << "Carry finish" << std::endl;
    return mres;
}

std::vector<Ciphertext<lbcrypto::DCRTPoly>> HomomorphicCarry(
    std::vector<Ciphertext<lbcrypto::DCRTPoly>>& ctxt, Ciphertext<lbcrypto::DCRTPoly>& ctxt_rshift_2Bg,
    std::function<double(double)> mod2Bg, size_t seg_after, size_t p, int Bg, int depth,
    std::vector<uint32_t> levelBudget, lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey, int numSlots) {
    auto cc           = ctxt[0]->GetCryptoContext();
    auto algo         = cc->GetScheme();
    int modLevel      = depth - levelBudget[1];
    int hermite_order = 1;
    std::vector<Plaintext> result;

    auto ModBoot = [&](Ciphertext<DCRTPoly> c) {
        auto boot = c->Clone();
        algo->LevelReduceInternalInPlace(boot, modLevel - boot->GetLevel());
        return cc->EvalFuncBootstrap(boot, mod2Bg, p, hermite_order);
    };

    // Step 3: Carry Propagation per level
    auto CarryStep = [&](Ciphertext<DCRTPoly> lower, Ciphertext<DCRTPoly> higher_mod) {
        auto high_part = cc->EvalSub(lower, higher_mod);
        algo->AdjustLevelsInPlace(high_part, ctxt_rshift_2Bg);
        auto high_div = cc->EvalMult(high_part, ctxt_rshift_2Bg);
        cc->RescaleInPlace(high_div);
        return high_div;
    };

    std::cout << "ctxt.size(): " << ctxt.size() << std::endl;
    std::cout << "seg_after: " << seg_after << std::endl;
    std::cout << "ModBoot" << std::endl;
    std::vector<Ciphertext<lbcrypto::DCRTPoly>> M_mod, M_high, t, t_mod;

    for (size_t i = 0; i < ctxt.size(); i++)
        M_mod.push_back(ModBoot(ctxt[i]));

    std::cout << "CarryStep" << std::endl;
    for (size_t i = 1; i < ctxt.size(); i++)
        M_high.push_back(CarryStep(ctxt[i], M_mod[i]));

    t.push_back(cc->EvalAdd(CarryStep(ctxt[0], M_mod[0]), M_mod[1]));

    for (size_t i = 1; i < ctxt.size() - 1; i++) {
        std::cout << "i:" << i << std::endl;
        t_mod.push_back(ModBoot(t[i - 1]));
        t.push_back(cc->EvalAdd(cc->EvalAdd(CarryStep(t[i - 1], t_mod[i - 1]), M_mod[i + 1]), M_high[i - 1]));
    }
    t_mod.push_back(ModBoot(t[ctxt.size() - 2]));
    std::cout << "t_mod size:" << t_mod.size() << ", t size:" << t.size() << std::endl;
    std::cout << "t_mod level:" << t_mod[t_mod.size() - 1]->GetLevel() << std::endl;
    // for (int i = 1; i < seg_after; i++) {
    //     Plaintext res;
    //     cc->Decrypt(privateKey, t[i], &result0);
    // }
    // result.push_back(cc->Decrypt(privateKey, t1, &result0));
    // cc->Decrypt(privateKey, t1, &result0);
    // cc->Decrypt(privateKey, t2, &result1);
    // cc->Decrypt(privateKey, t3, &result2);
    // result0->SetLength(numSlots);
    // result1->SetLength(numSlots);
    // result2->SetLength(numSlots);
    // auto finalResult0 = result0->GetCKKSPackedValue();
    // auto finalResult1 = result1->GetCKKSPackedValue();
    // auto finalResult2 = result2->GetCKKSPackedValue();
    // std::cout << "t1: " << finalResult0[0].real() << ", t2: " << finalResult1[0].real()
    //           << ", t3: " << finalResult2[0].real() << std::endl;
    // std::cout << "t4 level: " << mres4->GetLevel() << std::endl;
    // std::cout << "CarryStep finish" << std::endl;
    // Step 4: Return all segments mres0 ~ mres7
    std::vector<Ciphertext<DCRTPoly>> mres;
    mres.push_back(M_mod[0]);
    for (size_t i = 0; i < t.size(); i++) {
        mres.push_back(ModBoot(t[i]));
    }
    if (seg_after > ctxt.size()) {
        auto mres_last = cc->EvalAdd(CarryStep(t[ctxt.size() - 2], t_mod[ctxt.size() - 2]), M_high[ctxt.size() - 2]);
        mres.push_back(mres_last);
    }

    std::cout << "Carry finish" << std::endl;
    return mres;
}

std::vector<Ciphertext<lbcrypto::DCRTPoly>> HomomorphicCarryMult32_opt(
    Ciphertext<lbcrypto::DCRTPoly>& ctxt0, Ciphertext<lbcrypto::DCRTPoly>& ctxt1, Ciphertext<lbcrypto::DCRTPoly>& ctxt2,
    Ciphertext<lbcrypto::DCRTPoly>& ctxt3, Ciphertext<lbcrypto::DCRTPoly>& ctxt_rshift_2k,
    std::function<double(double)> mod2k, size_t p, int k, int depth, std::vector<uint32_t> levelBudget,
    lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey, int numSlots, std::vector<std::complex<double>> coe_mod) {
    auto cc   = ctxt0->GetCryptoContext();
    auto algo = cc->GetScheme();

    
    int modLevel = depth - levelBudget[1];

    std::cout << "=== Starting Optimized HomomorphicCarryMult32 ===" << std::endl;
    auto start_total = std::chrono::system_clock::now();

    // ----------------------------------------------------------------------
    // Step 1: Compute Raw Multiplications (M0 ~ M6)
    
    // ----------------------------------------------------------------------

    // M0 = c0 * c0
    auto M0 = cc->EvalMult(ctxt0, ctxt0);
    cc->RescaleInPlace(M0);

    // M1 = 2 * c0 * c1
    auto M1 = cc->EvalMult(ctxt0, ctxt1);
    cc->RescaleInPlace(M1);
    cc->EvalAddInPlace(M1, M1);  // double it

    // M2 = c1*c1 + 2*c0*c2
    auto M2a = cc->EvalMult(ctxt1, ctxt1);
    cc->RescaleInPlace(M2a);
    auto M2b = cc->EvalMult(ctxt0, ctxt2);
    cc->RescaleInPlace(M2b);
    cc->EvalAddInPlace(M2b, M2b);  // double
    auto M2 = cc->EvalAdd(M2a, M2b);

    // M3 = 2*c0*c3 + 2*c1*c2
    auto M3a = cc->EvalMult(ctxt0, ctxt3);
    cc->RescaleInPlace(M3a);
    cc->EvalAddInPlace(M3a, M3a);
    auto M3b = cc->EvalMult(ctxt1, ctxt2);
    cc->RescaleInPlace(M3b);
    cc->EvalAddInPlace(M3b, M3b);
    auto M3 = cc->EvalAdd(M3a, M3b);

    // M4 = c2*c2 + 2*c1*c3
    auto M4a = cc->EvalMult(ctxt2, ctxt2);
    cc->RescaleInPlace(M4a);
    auto M4b = cc->EvalMult(ctxt1, ctxt3);
    cc->RescaleInPlace(M4b);
    cc->EvalAddInPlace(M4b, M4b);
    auto M4 = cc->EvalAdd(M4a, M4b);

    // M5 = 2*c2*c3
    auto M5 = cc->EvalMult(ctxt2, ctxt3);
    cc->RescaleInPlace(M5);
    cc->EvalAddInPlace(M5, M5);

    // M6 = c3*c3
    auto M6 = cc->EvalMult(ctxt3, ctxt3);
    cc->RescaleInPlace(M6);

    
    std::vector<Ciphertext<DCRTPoly>> M = {M0, M1, M2, M3, M4, M5, M6};

    // ----------------------------------------------------------------------
    // Step 2: Parallel Modular Reduction (Get M_i mod 2^k)
    // ----------------------------------------------------------------------

    
    auto ModBoot = [&](Ciphertext<DCRTPoly> c) {
        auto boot = c->Clone();
        
        if (boot->GetLevel() < (size_t)modLevel) {
            algo->LevelReduceInternalInPlace(boot, modLevel - boot->GetLevel());
        }
        return cc->EvalFuncBootstrap_kx_parallel_opt_sparse(boot, p, 1, coe_mod, 1, false);
    };

    std::vector<Ciphertext<DCRTPoly>> M_mod(7);
    
    for (int i = 0; i < 7; ++i) {
        M_mod[i] = ModBoot(M[i]);
    }

    // ----------------------------------------------------------------------
    // Step 3: Pre-compute Static Carries (High Parts of M_i)
    // ----------------------------------------------------------------------

    
    
    auto CarryStep = [&](Ciphertext<DCRTPoly> raw, Ciphertext<DCRTPoly> modded) {
        auto diff = cc->EvalSub(raw, modded);
        algo->AdjustLevelsInPlace(diff, ctxt_rshift_2k);  
        auto carry = cc->EvalMult(diff, ctxt_rshift_2k);
        cc->RescaleInPlace(carry);
        return carry;
    };

    std::vector<Ciphertext<DCRTPoly>> M_high(7);
    
    
    for (int i = 0; i < 7; ++i) {
        M_high[i] = CarryStep(M[i], M_mod[i]);
    }

    // ----------------------------------------------------------------------
    
    // ----------------------------------------------------------------------
    
    

    std::vector<Ciphertext<DCRTPoly>> final_res(8);  // 0~7
    final_res[0] = M_mod[0];                         

    
    auto current_carry = M_high[0];

    for (int i = 1; i < 7; ++i) {

        
        
        
        

        // t_sum = M[i]_mod + Carry_from_prev_stage
        
        
        

        
        auto total_carry_in = (i == 1) ? M_high[0] : cc->EvalAdd(M_high[i - 1], current_carry);

        
        auto t = cc->EvalAdd(M_mod[i], total_carry_in);

        
        final_res[i] = ModBoot(t);

        
        
        if (i < 7) {  
            current_carry = CarryStep(t, final_res[i]);
        }
    }

    
    // result 7 = M6_high + carry_from_t6
    final_res[7] = cc->EvalAdd(M_high[6], current_carry);

    auto end_total  = std::chrono::system_clock::now();
    double total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total).count();
    std::cout << "HomomorphicCarryMult32 Finished in " << total_ms << " ms" << std::endl;

    return final_res;
}

void Probe(std::string label, Ciphertext<DCRTPoly> ctxt, double expected, lbcrypto::PrivateKey<DCRTPoly>& sk,
           int slot_idx = 0) {
    Plaintext p;
    ctxt->GetCryptoContext()->Decrypt(sk, ctxt, &p);
    p->SetLength(slot_idx + 1);
    double got = p->GetCKKSPackedValue()[slot_idx].real();

    double diff = std::abs(expected - got);
    std::cout << "  [" << label << "] Exp: " << std::fixed << std::setprecision(2) << expected << " | Got: " << got
              << " | Diff: " << diff;

    if (diff > 1.0)
        std::cout << " <--- ERROR HERE!";
    std::cout << std::endl;
}

std::vector<Ciphertext<DCRTPoly>> HomomorphicCarryMult4BitDynamic(
    const std::vector<Ciphertext<DCRTPoly>>& chunksA, const std::vector<Ciphertext<DCRTPoly>>& chunksB, size_t p,
    int depth, const std::vector<uint32_t>& levelBudget, const std::vector<std::complex<double>>& coe_mod16) {
    if (chunksA.empty() || chunksB.empty()) {
        OPENFHE_THROW("4-bit dynamic multiplication requires non-empty inputs");
    }

    auto cc             = chunksA[0]->GetCryptoContext();
    auto algo           = cc->GetScheme();
    size_t raw_count    = chunksA.size() + chunksB.size() - 1;
    int modLevel        = depth - levelBudget[1];
    constexpr double iB = 1.0 / 16.0;

    // mod 16 is periodic over every p that is a multiple of 16.  The raw
    // convolution coefficient may be wider than bootstrap_bits, but FBT only
    // needs its residue class modulo p to recover the same mod-16 digit.
    size_t carryChunkCount = std::min(chunksA.size(), chunksB.size());
    if (p != 256 || carryChunkCount > 64) {
        OPENFHE_THROW("4-bit carry multiplication expects p=256 and at most 64 input chunks");
    }

    auto BootMod16 = [&](const Ciphertext<DCRTPoly>& c) {
        auto boot = c->Clone();
        if (boot->GetLevel() < static_cast<size_t>(modLevel)) {
            algo->LevelReduceInternalInPlace(boot, modLevel - boot->GetLevel());
        }
        return cc->EvalFuncBootstrap_kx_parallel_opt_sparse(boot, p, 1, coe_mod16, 1, false);
    };

    // Refresh the carry as independent 4-bit digits. This keeps every
    // functional bootstrap in the p=256 domain even when the complete carry
    // is larger than 255. Up to 64 input chunks need at most three hex digits
    // because the carry is below 15 * 64 = 960.
    size_t carryDigits = (carryChunkCount <= 1) ? 1 : ((carryChunkCount <= 16) ? 2 : 3);
    auto BootSegmentedCarry = [&](const Ciphertext<DCRTPoly>& c) {
        std::vector<Ciphertext<DCRTPoly>> digits;
        digits.reserve(carryDigits);
        auto remaining = c->Clone();

        for (size_t digitIndex = 0; digitIndex < carryDigits; ++digitIndex) {
            auto digit = BootMod16(remaining);
            digits.push_back(digit);
            if (digitIndex + 1 == carryDigits)
                break;

            auto aligned = remaining->Clone();
            auto digitForSplit = digit->Clone();
            algo->AdjustLevelsInPlace(aligned, digitForSplit);
            auto high = cc->EvalSub(aligned, digitForSplit);
            remaining = cc->EvalMult(high, iB);
            cc->ModReduceInPlace(remaining);
        }

        auto refreshed = digits[0]->Clone();
        for (size_t digitIndex = 1; digitIndex < digits.size(); ++digitIndex) {
            auto term = MultByInteger(digits[digitIndex], int64_t{1} << (4 * digitIndex));
            algo->AdjustLevelsInPlace(refreshed, term);
            cc->EvalAddInPlace(refreshed, term);
        }
        return refreshed;
    };

    // Base-16 schoolbook convolution.  When both references are the same
    // vector (the square benchmark), compute each symmetric cross term once
    // and double it with an addition.
    bool square = &chunksA == &chunksB;
    std::vector<Ciphertext<DCRTPoly>> raw(raw_count);
    for (size_t k = 0; k < raw_count; ++k) {
        bool first = true;
        for (size_t i = 0; i < chunksA.size(); ++i) {
            if (k < i)
                continue;
            size_t j = k - i;
            if (j >= chunksB.size())
                continue;
            if (square && i > j)
                continue;

            auto term = cc->EvalMult(chunksA[i], chunksB[j]);
            cc->ModReduceInPlace(term);
            if (square && i != j)
                cc->EvalAddInPlace(term, term);

            if (first) {
                raw[k] = term;
                first  = false;
            }
            else {
                cc->EvalAddInPlace(raw[k], term);
            }
        }
    }

    auto CarryStep = [&](Ciphertext<DCRTPoly> value, Ciphertext<DCRTPoly> digit) {
        auto aligned = value->Clone();
        algo->AdjustLevelsInPlace(aligned, digit);
        auto high = cc->EvalSub(aligned, digit);
        auto carry = cc->EvalMult(high, iB);
        cc->ModReduceInPlace(carry);
        return carry;
    };

    std::vector<Ciphertext<DCRTPoly>> result(raw_count + 1);
    Ciphertext<DCRTPoly> carry;
    for (size_t k = 0; k < raw_count; ++k) {
        auto current = raw[k];
        if (k != 0) {
            if (carry->GetLevel() + 2 >= static_cast<size_t>(modLevel)) {
                std::cout << "[CarryRefresh] chunk " << k << ", input level " << carry->GetLevel() << std::endl;
                carry = BootSegmentedCarry(carry);
            }
            algo->AdjustLevelsInPlace(current, carry);
            cc->EvalAddInPlace(current, carry);
        }

        result[k] = BootMod16(current);
        carry = CarryStep(current, result[k]);
    }

    result[raw_count] = carry;
    return result;
}

// Add two little-endian base-16 ciphertext vectors and normalize every chunk.
// result[k] is the k-th base-16 chunk and result.back() is the final carry.
std::vector<Ciphertext<DCRTPoly>> HomomorphicCarryAdd4BitDynamic(
    const std::vector<Ciphertext<DCRTPoly>>& chunksA, const std::vector<Ciphertext<DCRTPoly>>& chunksB, size_t p,
    int depth, const std::vector<uint32_t>& levelBudget, const std::vector<std::complex<double>>& coe_mod16,
    const Ciphertext<DCRTPoly>& ctxt_rshift_16) {
    if (chunksA.empty() || chunksB.empty()) {
        OPENFHE_THROW("4-bit dynamic addition requires non-empty inputs");
    }

    auto cc             = chunksA[0]->GetCryptoContext();
    auto algo           = cc->GetScheme();
    size_t raw_count    = std::max(chunksA.size(), chunksB.size());
    int modLevel        = depth - levelBudget[1];

    auto BootMod16 = [&](const Ciphertext<DCRTPoly>& c) {
        auto boot = c->Clone();
        if (boot->GetLevel() < static_cast<size_t>(modLevel)) {
            algo->LevelReduceInternalInPlace(boot, modLevel - boot->GetLevel());
        }
        return cc->EvalFuncBootstrap_kx_parallel_opt_sparse(boot, p, 1, coe_mod16, 1, false);
    };

    // The add/subtract and encrypted right shift below need two levels. Refresh
    // the carry before it enters the levels reserved by functional bootstrap.
    auto RefreshCarryIfNeeded = [&](const Ciphertext<DCRTPoly>& c, size_t chunk) {
        constexpr size_t levelsNeeded = 2;
        size_t level                  = c->GetLevel();
        size_t levelsLeft = level < static_cast<size_t>(modLevel) ? static_cast<size_t>(modLevel) - level : 0;
        if (levelsLeft <= levelsNeeded) {
            std::cout << "[CarryRefresh] chunk " << chunk << ", input level " << level << std::endl;
            // An addition carry is 0 or 1, so mod 16 is an identity on it.
            return BootMod16(c);
        }
        return c->Clone();
    };

    std::vector<Ciphertext<DCRTPoly>> result(raw_count + 1);
    Ciphertext<DCRTPoly> carry;

    for (size_t k = 0; k < raw_count; ++k) {
        Ciphertext<DCRTPoly> current;
        if (k < chunksA.size() && k < chunksB.size()) {
            current = cc->EvalAdd(chunksA[k], chunksB[k]);
        }
        else if (k < chunksA.size()) {
            current = chunksA[k]->Clone();
        }
        else {
            current = chunksB[k]->Clone();
        }

        // Add the carry generated by the previous, less-significant chunk.
        if (k != 0) {
            carry = RefreshCarryIfNeeded(carry, k);
            algo->AdjustLevelsInPlace(current, carry);
            cc->EvalAddInPlace(current, carry);
        }

        // digit = current mod 16
        result[k] = BootMod16(current);

        // carry = (current - digit) / 16. As in the fixed carry example,
        // multiply by an encrypted 1/16 so the CRT basis is explicitly
        // aligned before multiplication.
        auto current_aligned = current->Clone();
        algo->AdjustLevelsInPlace(current_aligned, result[k]);
        auto high = cc->EvalSub(current_aligned, result[k]);
        auto shift = ctxt_rshift_16->Clone();
        algo->AdjustLevelsInPlace(high, shift);
        carry = cc->EvalMult(high, shift);
        cc->RescaleInPlace(carry);
    }

    result[raw_count] = carry;
    return result;
}

std::vector<Ciphertext<DCRTPoly>> HomomorphicCarryMult8BitDynamic_Debug(
    const std::vector<Ciphertext<DCRTPoly>>& chunksA, const std::vector<Ciphertext<DCRTPoly>>& chunksB,
    
    const std::vector<std::vector<double>>& plainA, const std::vector<std::vector<double>>& plainB, size_t p, int depth,
    std::vector<uint32_t> levelBudget, lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey,
    std::vector<std::complex<double>> coe_mod) {
    auto cc             = chunksA[0]->GetCryptoContext();
    auto algo           = cc->GetScheme();
    size_t total_chunks = chunksA.size() + chunksB.size() - 1;

    
    int DEBUG_SLOT = 0;
    std::cout << "\n>>> DEBUG TRACE START (Slot " << DEBUG_SLOT << ") <<<" << std::endl;

    double c_inv16  = 1.0 / 16.0;
    double c_inv256 = 1.0 / 256.0;
    int modLevel    = depth - levelBudget[1];

    auto Boot16 = [&](Ciphertext<DCRTPoly> c) {
        auto boot = c->Clone();
        if (boot->GetLevel() < (size_t)modLevel)
            algo->LevelReduceInternalInPlace(boot, modLevel - boot->GetLevel());
        return cc->EvalFuncBootstrap_kx_parallel_opt_sparse(boot, p, 1, coe_mod, 1, false);
    };

    // ---------------------------------------------------------
    
    // ---------------------------------------------------------
    std::vector<Ciphertext<DCRTPoly>> M_raw(total_chunks);
    std::vector<double> M_raw_expected(total_chunks, 0.0);

    for (size_t k = 0; k < total_chunks; ++k) {
        
        double expected_val = 0;
        for (size_t i = 0; i < chunksA.size(); ++i) {
            size_t j = k - i;
            if (k >= i && j < chunksB.size()) {
                expected_val += plainA[i][DEBUG_SLOT] * plainB[j][DEBUG_SLOT];
            }
        }
        M_raw_expected[k] = expected_val;

        
        Ciphertext<DCRTPoly> sum_k;
        bool first = true;
        for (size_t i = 0; i < chunksA.size(); ++i) {
            size_t j = k - i;
            if (k >= i && j < chunksB.size()) {
                auto term = cc->EvalMult(chunksA[i], chunksB[j]);
                cc->ModReduceInPlace(term);  
                if (first) {
                    sum_k = term;
                    first = false;
                }
                else {
                    cc->EvalAddInPlace(sum_k, term);
                }
            }
        }
        M_raw[k] = sum_k;

        
        // std::cout << "Step 1: Convolution Chunk " << k << std::endl;
        // Probe("M_raw", M_raw[k], expected_val, privateKey, DEBUG_SLOT);
    }

    // ---------------------------------------------------------
    
    // ---------------------------------------------------------
    std::vector<Ciphertext<DCRTPoly>> final_results(total_chunks + 1);
    Ciphertext<DCRTPoly> current_carry;
    double current_carry_expected = 0.0;

    for (size_t i = 0; i < total_chunks; ++i) {
        std::cout << "\nStep 2: Processing Chunk " << i << std::endl;

        
        Ciphertext<DCRTPoly> current_val = M_raw[i];
        double current_val_exp           = M_raw_expected[i];

        if (i > 0) {
            algo->AdjustLevelsInPlace(current_val, current_carry);
            current_val = cc->EvalAdd(current_val, current_carry);
            current_val_exp += current_carry_expected;

            
            // Probe("Val+Carry", current_val, current_val_exp, privateKey, DEBUG_SLOT);
        }

        // ==========================================
        
        // ==========================================

        // A. r_low (x % 16)
        auto r_low = Boot16(current_val);
        
        double exp_r_low = std::fmod(current_val_exp, 16.0);
        
        if (exp_r_low < 0)
            exp_r_low += 16.0;
        
        exp_r_low = std::round(exp_r_low);

        // Probe("A. r_low", r_low, exp_r_low, privateKey, DEBUG_SLOT);

        // B. x_shifted ((x - r_low) / 16)
        auto x_raw_aligned = current_val->Clone();  
        auto diff_low      = cc->EvalSub(x_raw_aligned, r_low);
        auto x_shifted     = cc->EvalMult(diff_low, c_inv16);
        cc->ModReduceInPlace(x_shifted);  

        double exp_shifted = (current_val_exp - exp_r_low) / 16.0;
        // Probe("B. x_shifted", x_shifted, exp_shifted, privateKey, DEBUG_SLOT);

        // C. r_high (x_shifted % 16)
        auto r_high       = Boot16(x_shifted);
        double exp_r_high = std::fmod(exp_shifted, 16.0);
        if (exp_r_high < 0)
            exp_r_high += 16.0;
        exp_r_high = std::round(exp_r_high);

        // Probe("C. r_high", r_high, exp_r_high, privateKey, DEBUG_SLOT);

        // D. Result (r_low + 16 * r_high)
        auto term_high = cc->EvalMult(r_high, 16);
        cc->ModReduceInPlace(term_high);  

        auto r_low_deep = r_low->Clone();
        algo->AdjustLevelsInPlace(r_low_deep, term_high);
        auto res_mod256 = cc->EvalAdd(r_low_deep, term_high);

        double exp_res = exp_r_low + 16.0 * exp_r_high;
        // Probe("D. res_256", res_mod256, exp_res, privateKey, DEBUG_SLOT);

        // E. Carry Out ((x - res) / 256)
        auto x_raw_deep = current_val->Clone();  
        algo->AdjustLevelsInPlace(x_raw_deep, res_mod256);
        auto diff_total = cc->EvalSub(x_raw_deep, res_mod256);
        auto carry_out  = cc->EvalMult(diff_total, c_inv256);
        cc->ModReduceInPlace(carry_out);  

        double exp_carry = (current_val_exp - exp_res) / 256.0;
        // Probe("E. carry", carry_out, exp_carry, privateKey, DEBUG_SLOT);

        
        final_results[i]       = res_mod256;
        current_carry          = carry_out;
        current_carry_expected = exp_carry;
    }

    final_results[total_chunks] = current_carry;
    return final_results;
}

void DebugDecrypt(const std::string& label, Ciphertext<lbcrypto::DCRTPoly> ctxt,
                  const lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& sk, double expected = -99999.0) {
    auto cc = ctxt->GetCryptoContext();
    Plaintext result;
    cc->Decrypt(sk, ctxt, &result);
    double val     = result->GetCKKSPackedValue()[0].real();
    double int_val = std::round(val);

    std::cout << " [DEBUG] " << label << ": " << int_val;
    if (expected != -99999.0) {
        std::cout << " (Expected: " << expected << ")";
        if (std::abs(int_val - expected) > 0.5)
            std::cout << " [FAIL]";
        else
            std::cout << " [OK]";
    }
    std::cout << std::endl;
}  // degree = 128

std::vector<Ciphertext<lbcrypto::DCRTPoly>> Homomorphic_lazy_carry(
    std::vector<std::vector<Ciphertext<lbcrypto::DCRTPoly>>>& ctxt, size_t seg_after, size_t p, int Bg, int depth,
    std::vector<uint32_t> levelBudget, lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey, int numSlots, int vec_size,
    double expect_S0, double expect_S1, double expect_S2) {
    auto cc      = ctxt[0][0]->GetCryptoContext();
    auto algo    = cc->GetScheme();
    int dim      = ctxt.size();
    int modLevel = depth - levelBudget[1];
    std::vector<std::complex<double>> coe_1, coe_2, coe_3, coe_4, coe_rs_4bit;
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

    coe_1       = hermiteInterpOrder1Coeffshift([f_1](double x) -> double { return 0.5 * f_1(x); }, p);
    coe_2       = hermiteInterpOrder1Coeffshift([f_2](double x) -> double { return 0.5 * f_2(x); }, p);
    coe_3       = hermiteInterpOrder1Coeffshift([f_3](double x) -> double { return 0.5 * f_3(x); }, p);
    coe_4       = hermiteInterpOrder1Coeffshift([f_4](double x) -> double { return 0.5 * f_4(x); }, p);
    coe_rs_4bit = hermiteInterpOrder1Coeffshift([f_rs_4bit](double x) -> double { return 0.5 * f_rs_4bit(x); }, p);
    std::vector<std::vector<std::complex<double>>> coe_1234 = {coe_1, coe_2, coe_3, coe_4};
    std::vector<std::vector<std::complex<double>>> coe_1_rs = {coe_1, coe_rs_4bit};

    std::vector<std::vector<Ciphertext<lbcrypto::DCRTPoly>>> csquare(3,
                                                                     std::vector<Ciphertext<lbcrypto::DCRTPoly>>(dim));
    std::vector<Ciphertext<lbcrypto::DCRTPoly>> csquare_sum(3);
    for (int i = 0; i < dim; i++) {
        csquare[0][i] = cc->EvalMult(ctxt[i][0], ctxt[i][0]);
        cc->RescaleInPlace(csquare[0][i]);
        csquare[1][i] = cc->EvalMult(ctxt[i][0], ctxt[i][1]);
        cc->RescaleInPlace(csquare[1][i]);
        csquare[2][i] = cc->EvalMult(ctxt[i][1], ctxt[i][1]);
        cc->RescaleInPlace(csquare[2][i]);
    }
    csquare_sum[0] = cc->EvalAddMany(csquare[0]);
    csquare_sum[1] = cc->EvalAddMany(csquare[1]);
    csquare_sum[2] = cc->EvalAddMany(csquare[2]);
    csquare_sum[1] = MultByInteger(csquare_sum[1], 2);
    std::cout << "\n--- [DEBUG] Aggregation Check ---" << std::endl;
    DebugDecrypt("S0 (Sum x0*x0)", csquare_sum[0], privateKey, expect_S0);
    DebugDecrypt("S1 (Sum 2*x0*x1)", csquare_sum[1], privateKey, expect_S1);
    DebugDecrypt("S2 (Sum x1*x1)", csquare_sum[2], privateKey, expect_S2);
    int opt_depth = 6;
    bool use_img  = false;
    // decompose as 16bit (up to 14bits), each part is 4bit
    auto MVBoot = [&](Ciphertext<DCRTPoly> c, std::vector<std::vector<std::complex<double>>> coe) {
        auto boot = c->Clone();
        algo->LevelReduceInternalInPlace(boot, modLevel - boot->GetLevel());
        return cc->EvalMVFuncBootstrap_kx_parallel_opt(boot, p, 1, coe, opt_depth, use_img);
    };

    //
    // auto csquare_pack  = cc->EvalAdd(csquare_sum[0], cc->EvalRotate());
    auto c1_square_seg = MVBoot(csquare_sum[0], coe_1234);
    auto c2_square_seg = MVBoot(csquare_sum[1], coe_1234);
    auto c3_square_seg = MVBoot(csquare_sum[2], coe_1234);

    // --- Carry Propagation ---
    std::cout << "\n--- [DEBUG] Carry Propagation Check ---" << std::endl;

    // seg_0
    auto mres0 = c1_square_seg[0];
    DebugDecrypt("MRes0 (Seg 0)", mres0, privateKey);

    // seg_1 = (c1[1] + c2[0]) % Bg
    auto temp1 = cc->EvalAdd(c1_square_seg[1], c2_square_seg[0]);
    DebugDecrypt("Temp1 (Pre-Carry Seg1)", temp1, privateKey);  
    // auto temp2 = cc->EvalMVFuncBootstrap_kx_parallel_opt(temp1, p, 1, coe_1_rs, opt_depth, use_img);

    std::cout << "temp1 noise_deg:" << temp1->GetNoiseScaleDeg() << std::endl;
    auto temp2 = MVBoot(temp1, coe_1_rs);
    std::cout << "temp2[0]:" << temp2[0]->GetLevel() << std::endl;
    auto mres1 = temp2[0];
    DebugDecrypt("MRes1 (Seg 1)", mres1, privateKey);

    // seg_2 = (c1[2] + c2[1] + c3[0] + carry_from_seg1) % Bg
    auto temp3 = cc->EvalAddMany({c1_square_seg[2], c2_square_seg[1], c3_square_seg[0], temp2[1]});
    DebugDecrypt("Temp3 (Pre-Carry Seg2)", temp3, privateKey);
    // auto temp4 = cc->EvalMVFuncBootstrap_kx_parallel_opt(temp3, p, 1, coe_1_rs, opt_depth, use_img);

    std::cout << "temp3 noise_deg:" << temp3->GetNoiseScaleDeg() << std::endl;
    auto temp4 = MVBoot(temp3, coe_1_rs);
    auto mres2 = temp4[0];
    DebugDecrypt("MRes2 (Seg 2)", mres2, privateKey);

    // seg_3
    auto temp5 = cc->EvalAddMany({c1_square_seg[3], c2_square_seg[2], c3_square_seg[1], temp4[1]});

    std::cout << "temp5 noise_deg:" << temp5->GetNoiseScaleDeg() << std::endl;
    DebugDecrypt("Temp5 (Pre-Carry Seg3)", temp5, privateKey);
    // auto temp6 = cc->EvalMVFuncBootstrap_kx_parallel_opt(temp5, p, 1, coe_1_rs, opt_depth, use_img);
    auto temp6 = MVBoot(temp5, coe_1_rs);
    auto mres3 = temp6[0];
    DebugDecrypt("MRes3 (Seg 3)", mres3, privateKey);

    // seg_4
    auto temp7 = cc->EvalAddMany({c2_square_seg[3], c3_square_seg[2], temp6[1]});
    DebugDecrypt("Temp7 (Pre-Carry Seg4)", temp7, privateKey);
    // auto temp8 = cc->EvalMVFuncBootstrap_kx_parallel_opt(temp7, p, 1, coe_1_rs, opt_depth, use_img);
    auto temp8 = MVBoot(temp7, coe_1_rs);
    auto mres4 = temp8[0];
    DebugDecrypt("MRes4 (Seg 4)", mres4, privateKey);

    // seg_5
    auto temp9 = cc->EvalAdd(c3_square_seg[3], temp8[1]);
    DebugDecrypt("Temp9 (Pre-Carry Seg5)", temp9, privateKey);
    // auto temp10 = cc->EvalMVFuncBootstrap_kx_parallel_opt(temp9, p, 1, coe_1_rs, opt_depth, use_img);
    auto temp10 = MVBoot(temp9, coe_1_rs);
    auto mres5  = temp10[0];
    DebugDecrypt("MRes5 (Seg 5)", mres5, privateKey);

    return std::vector<Ciphertext<lbcrypto::DCRTPoly>>{mres0, mres1, mres2, mres3, mres4, mres5};
}

std::vector<Ciphertext<lbcrypto::DCRTPoly>> Homomorphic_lazy_carry_iterate(
    std::vector<std::vector<Ciphertext<lbcrypto::DCRTPoly>>>& ctxt, Ciphertext<lbcrypto::DCRTPoly> ctxt_rshift_16,
    std::vector<std::complex<double>> coe_mod16, size_t seg_after, size_t p, int Bg, int depth,
    std::vector<uint32_t> levelBudget, lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey, int numSlots, int vec_size,
    double expect_S0, double expect_S1, double expect_S2) {
    auto cc   = ctxt[0][0]->GetCryptoContext();
    auto algo = cc->GetScheme();
    int dim   = ctxt.size();

    std::vector<std::vector<Ciphertext<lbcrypto::DCRTPoly>>> csquare(3,
                                                                     std::vector<Ciphertext<lbcrypto::DCRTPoly>>(dim));
    std::vector<Ciphertext<lbcrypto::DCRTPoly>> csquare_sum(3);

    for (int i = 0; i < dim; i++) {
        csquare[0][i] = cc->EvalMult(ctxt[i][0], ctxt[i][0]);
        cc->RescaleInPlace(csquare[0][i]);
        csquare[1][i] = cc->EvalMult(ctxt[i][0], ctxt[i][1]);
        cc->RescaleInPlace(csquare[1][i]);
        csquare[2][i] = cc->EvalMult(ctxt[i][1], ctxt[i][1]);
        cc->RescaleInPlace(csquare[2][i]);
    }

    
    csquare_sum[0] = cc->EvalAddMany(csquare[0]);
    csquare_sum[1] = cc->EvalAddMany(csquare[1]);
    csquare_sum[2] = cc->EvalAddMany(csquare[2]);
    csquare_sum[1] = MultByInteger(csquare_sum[1], 2);

    // std::cout << "\n--- [DEBUG] Aggregation Check ---" << std::endl;
    // DebugDecrypt("S0 (Sum x0*x0)", csquare_sum[0], privateKey, expect_S0);
    // DebugDecrypt("S1 (Sum 2*x0*x1)", csquare_sum[1], privateKey, expect_S1);
    // DebugDecrypt("S2 (Sum x1*x1)", csquare_sum[2], privateKey, expect_S2);

    // =================================================================================
    
    // =================================================================================
    int modLevel  = depth - levelBudget[1];
    int opt_depth = 6;
    bool use_img  = false;

    
    auto BootMod16 = [&](Ciphertext<DCRTPoly> c) {
        // std::cout << "BootMod16 start" << std::endl;
        auto boot = c->Clone();
        

        // std::cout << "level before boot adjust:" << boot->GetLevel() << std::endl;
        if (boot->GetLevel() < size_t(modLevel)) {
            algo->LevelReduceInternalInPlace(boot, modLevel - boot->GetLevel());
        }
        

        // std::cout << "level after boot adjust:" << boot->GetLevel() << std::endl;
        auto res_vec = cc->EvalFuncBootstrap_kx_parallel_opt_sparse(boot, p, 1, coe_mod16, opt_depth, use_img);
        // std::cout << "BootMod16 finish" << std::endl;
        return res_vec;
    };

    // std::cout << "level ctxt_rshift_16:" << ctxt_rshift_16->GetLevel() << std::endl;
    
    
    auto CarryStep = [&](Ciphertext<DCRTPoly> lower, Ciphertext<DCRTPoly> higher_mod) {
        auto high_part = cc->EvalSub(lower, higher_mod);

        // std::cout << "CarryStep start" << std::endl;
        // std::cout << "level high_part:" << high_part->GetLevel() << std::endl;
        
        algo->AdjustLevelsInPlace(high_part, ctxt_rshift_16);
        // std::cout << "level high_part adjust:" << high_part->GetLevel() << std::endl;
        // std::cout << "level high_part adjust:" << ctxt_rshift_16->GetLevel() << std::endl;

        
        auto high_div = cc->EvalMult(high_part, ctxt_rshift_16);
        cc->RescaleInPlace(high_div);
        // std::cout << "CarryStep finish" << std::endl;
        return high_div;
    };

    // =================================================================================
    
    // =================================================================================
    std::cout << "\n--- Iterative Carry Propagation ---" << std::endl;

    std::vector<Ciphertext<lbcrypto::DCRTPoly>> results;

    
    Ciphertext<DCRTPoly> accumulator = csquare_sum[0];

    
    for (int i = 0; i < 5; i++) {
        std::cout << "Processing Seg " << i << "..." << std::endl;

        
        // MRes = Accumulator % 16
        auto mres = BootMod16(accumulator);
        results.push_back(mres);
        // DebugDecrypt("MRes" + std::to_string(i), mres, privateKey);

        
        if (i == 5)
            break;

        
        // Carry = (Accumulator - MRes) / 16

        // std::cout << "CarryStep " << i << "..." << std::endl;
        auto carry = CarryStep(accumulator, mres);

        
        
        if (i == 0) {
            
            accumulator = cc->EvalAdd(carry, csquare_sum[1]);
            // DebugDecrypt("Accumulator (Pre-Seg " + std::to_string(i + 1) + ")", accumulator, privateKey);
        }
        else if (i == 1) {
            
            accumulator = cc->EvalAdd(carry, csquare_sum[2]);
            // DebugDecrypt("Accumulator (Pre-Seg " + std::to_string(i + 1) + ")", accumulator, privateKey);
        }
        else {
            
            accumulator = carry;
            // DebugDecrypt("Accumulator (Pre-Seg " + std::to_string(i + 1) + ")", accumulator, privateKey);
        }
    }

    return results;
}

std::vector<std::vector<Ciphertext<lbcrypto::DCRTPoly>>> Homdist(
    std::vector<std::vector<Ciphertext<lbcrypto::DCRTPoly>>>& ev0,
    std::vector<std::vector<Ciphertext<lbcrypto::DCRTPoly>>>& ev1) {
    auto cc  = ev1[0][0]->GetCryptoContext();
    auto dim = ev1.size();  //degree
    auto w   = ev1[0].size();
    std::vector<std::vector<Ciphertext<lbcrypto::DCRTPoly>>> sub(dim, std::vector<Ciphertext<lbcrypto::DCRTPoly>>(w));
    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < w; j++) {
            sub[i][j] = cc->EvalSub(ev0[i][j], ev1[i][j]);
        }
    }
    return sub;
}

void CarrymulExample32bits(int bits, int p, int k) {
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

    /*usint depth = levelBudget[0] + levelBudget[1] + 12 + std::log2(p);*/
    // usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 4 + bits + 8;

    usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 4 + 6 + 3 + 8 + 4;
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    cc->Enable(FBTS);
    
    auto f_mod = [](double x) -> double {
        return x - 16.0 * std::floor((x + 8.0) / 16.0);  
    };
    
    std::vector<std::complex<double>> coe_mod =
        hermiteInterpOrder1Coeffshift([f_mod](double x) -> double { return 0.5 * f_mod(x); }, p);

    usint ringDim = cc->GetRingDimension();
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;
    std::cout << "8 bit shift " << std::endl;
    std::cout << "FBS Precision " << bits << " bits" << std::endl;
    std::cout << "Seg precision " << k << std::endl << std::endl;

    std::cout << "depth " << depth << std::endl << std::endl;
    // Precomputations for bootstrapping
    int numSlots = ringDim / 2;
    cc->EvalFuncBootstrapSetup(levelBudget, bsgsDim, numSlots, bits);

    // Key Generation
    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(-(pow(2, k) - 1), pow(2, k) - 1);
    // Input vector [0, 1, ..., p-2, p-1]
    std::vector<std::complex<double>> x0, x1, x2, x3;
    std::vector<double> x_rshift_16(numSlots, double(1.0 / 16));
    // 8bit squared x0+x1*2^k
    for (int i = 0; i < numSlots; i++) {
        // x.push_back((std::complex<double>(i % p, p - i % p - 1)));
        x0.push_back((std::complex<double>(message(engine))));
        x1.push_back((std::complex<double>(message(engine))));
        x2.push_back((std::complex<double>(message(engine))));
        x3.push_back((std::complex<double>(message(engine))));
    }
    Plaintext ptxt0          = cc->MakeCKKSPackedPlaintext(x0, 1);
    Plaintext ptxt1          = cc->MakeCKKSPackedPlaintext(x1, 1);
    Plaintext ptxt2          = cc->MakeCKKSPackedPlaintext(x2, 1);
    Plaintext ptxt3          = cc->MakeCKKSPackedPlaintext(x3, 1);
    Plaintext ptxt_rshift_16 = cc->MakeCKKSPackedPlaintext(x_rshift_16);
    ptxt0->SetLength(numSlots);
    ptxt1->SetLength(numSlots);
    ptxt2->SetLength(numSlots);
    ptxt3->SetLength(numSlots);
    ptxt_rshift_16->SetLength(numSlots);

    Ciphertext<DCRTPoly> ctxt0          = cc->Encrypt(keyPair.publicKey, ptxt0);
    Ciphertext<DCRTPoly> ctxt1          = cc->Encrypt(keyPair.publicKey, ptxt1);
    Ciphertext<DCRTPoly> ctxt2          = cc->Encrypt(keyPair.publicKey, ptxt2);
    Ciphertext<DCRTPoly> ctxt3          = cc->Encrypt(keyPair.publicKey, ptxt3);
    Ciphertext<DCRTPoly> ctxt_rshift_16 = cc->Encrypt(keyPair.publicKey, ptxt_rshift_16);

    auto mod16 = [](double x) -> double {
        return (double)((int)x % 16);
    };
    std::chrono::system_clock::time_point start, end;
    start     = std::chrono::system_clock::now();
    auto cres = HomomorphicCarryMult32_opt(ctxt0, ctxt1, ctxt2, ctxt3, ctxt_rshift_16, mod16, p, k, depth, levelBudget,
                                           keyPair.secretKey, numSlots, coe_mod);
    end       = std::chrono::system_clock::now();

    Plaintext result0, result1, result2, result3, result4, result5, result6, result7;

    cc->Decrypt(keyPair.secretKey, cres[0], &result0);
    cc->Decrypt(keyPair.secretKey, cres[1], &result1);
    cc->Decrypt(keyPair.secretKey, cres[2], &result2);
    cc->Decrypt(keyPair.secretKey, cres[3], &result3);
    cc->Decrypt(keyPair.secretKey, cres[4], &result4);
    cc->Decrypt(keyPair.secretKey, cres[5], &result5);
    cc->Decrypt(keyPair.secretKey, cres[6], &result6);
    cc->Decrypt(keyPair.secretKey, cres[7], &result7);
    result0->SetLength(numSlots);
    result1->SetLength(numSlots);
    result2->SetLength(numSlots);
    result3->SetLength(numSlots);
    result4->SetLength(numSlots);
    result5->SetLength(numSlots);
    result6->SetLength(numSlots);
    result7->SetLength(numSlots);
    auto finalResult0 = result0->GetCKKSPackedValue();
    auto finalResult1 = result1->GetCKKSPackedValue();
    auto finalResult2 = result2->GetCKKSPackedValue();
    auto finalResult3 = result3->GetCKKSPackedValue();
    auto finalResult4 = result4->GetCKKSPackedValue();
    auto finalResult5 = result5->GetCKKSPackedValue();
    auto finalResult6 = result6->GetCKKSPackedValue();
    auto finalResult7 = result7->GetCKKSPackedValue();
    int wrong         = 0;
    std::cout << "x0: " << x0[0].real() << ", x1: " << x1[0].real() << ", cres0: " << finalResult0[0].real()
              << ", cres1: " << finalResult1[0].real() << ", cres2: " << finalResult2[0].real()
              << ", cres3: " << finalResult3[0].real() << std::endl;
    for (int i = 0; i < numSlots; i++) {
        auto M = int(x0[i].real()) + int(x1[i].real()) * pow(2, k) + int(x2[i].real()) * pow(2, 2 * k) +
                 int(x3[i].real()) * pow(2, 3 * k);
        auto res = int(finalResult0[i].real()) + int(finalResult1[i].real()) * pow(2, k) +
                   int(finalResult2[i].real()) * pow(2, 2 * k) + int(finalResult3[i].real()) * pow(2, 3 * k) +
                   int(finalResult4[i].real()) * pow(2, 4 * k) + int(finalResult5[i].real()) * pow(2, 5 * k) +
                   int(finalResult6[i].real()) * pow(2, 6 * k) + int(finalResult7[i].real()) * pow(2, 7 * k);

        if (M * M != res) {
            std::cout << "M*M : " << M * M << ", res:" << res << std::endl;
            std::cout << "x0: " << x0[i].real() << ", x1: " << x1[i].real() << ", x2: " << x2[i].real()
                      << ", x3: " << x3[i].real() << ", cres0: " << finalResult0[i].real()
                      << ", cres1: " << finalResult1[i].real() << ", cres2: " << finalResult2[i].real()
                      << ", cres3: " << finalResult3[i].real() << ", cres4: " << finalResult4[i].real()
                      << ", cres5: " << finalResult5[i].real() << ", cres6: " << finalResult6[i].real()
                      << ", cres7: " << finalResult7[i].real() << std::endl;
            wrong++;
        }
    }
    std::cout << "error num:" << wrong << std::endl;
    double totaltime = 0;
    totaltime += double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / double(numSlots));
    std::cout << "shift time amortized:" << totaltime << "ms" << std::endl;
}

void CarrymulExample16bits(int bits, int p, int k) {
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

    /*usint depth = levelBudget[0] + levelBudget[1] + 12 + std::log2(p);*/
    usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 4 + bits + 8;
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    cc->Enable(FBTS);

    usint ringDim = cc->GetRingDimension();
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;
    std::cout << "8 bit shift " << std::endl;
    std::cout << "FBS Precision " << bits << " bits" << std::endl;
    std::cout << "Seg precision " << k << std::endl << std::endl;

    std::cout << "depth " << depth << std::endl << std::endl;
    // Precomputations for bootstrapping
    int numSlots = ringDim / 2;
    cc->EvalFuncBootstrapSetup(levelBudget, bsgsDim, numSlots, bits);

    // Key Generation
    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, pow(2, k) - 1);
    // Input vector [0, 1, ..., p-2, p-1]
    std::vector<double> x0, x1, y0, y1;
    std::vector<double> x_rshift_16(numSlots, double(1.0 / 16));
    // 8bit squared x0+x1*2^k
    for (int i = 0; i < numSlots; i++) {
        // x.push_back((std::complex<double>(i % p, p - i % p - 1)));
        // x0.push_back((message(engine)));
        // x1.push_back((message(engine)));
        // y0.push_back((message(engine)));
        // y1.push_back((message(engine)));
        x0.push_back(13);
        x1.push_back(3);
        y0.push_back(7);
        y1.push_back(0);
    }
    Plaintext ptxt0          = cc->MakeCKKSPackedPlaintext(x0, 1);
    Plaintext ptxt1          = cc->MakeCKKSPackedPlaintext(x1, 1);
    Plaintext ptxt_rshift_16 = cc->MakeCKKSPackedPlaintext(x_rshift_16);
    Plaintext ptxt2          = cc->MakeCKKSPackedPlaintext(y0, 1);
    Plaintext ptxt3          = cc->MakeCKKSPackedPlaintext(y1, 1);
    // ptxt0->SetLength(numSlots);
    // ptxt1->SetLength(numSlots);
    // ptxt_rshift_16->SetLength(numSlots);

    Ciphertext<DCRTPoly> ctxt0          = cc->Encrypt(keyPair.publicKey, ptxt0);
    Ciphertext<DCRTPoly> ctxt1          = cc->Encrypt(keyPair.publicKey, ptxt1);
    Ciphertext<DCRTPoly> ctxt_rshift_16 = cc->Encrypt(keyPair.publicKey, ptxt_rshift_16);
    Ciphertext<DCRTPoly> ctxt2          = cc->Encrypt(keyPair.publicKey, ptxt2);
    Ciphertext<DCRTPoly> ctxt3          = cc->Encrypt(keyPair.publicKey, ptxt3);

    auto mod16 = [](double x) -> double {
        return (double)((int)x % 16);
    };
    auto csub0 = cc->EvalSub(ctxt0, ctxt2);
    auto csub1 = cc->EvalSub(ctxt1, ctxt3);
    std::chrono::system_clock::time_point start, end;
    start     = std::chrono::system_clock::now();
    auto cres = HomomorphicCarryMult16(csub0, csub1, ctxt_rshift_16, mod16, p, k, depth, levelBudget, keyPair.secretKey,
                                       numSlots);
    end       = std::chrono::system_clock::now();

    Plaintext result0, result1, result2, result3;

    cc->Decrypt(keyPair.secretKey, cres[0], &result0);
    cc->Decrypt(keyPair.secretKey, cres[1], &result1);
    cc->Decrypt(keyPair.secretKey, cres[2], &result2);
    cc->Decrypt(keyPair.secretKey, cres[3], &result3);
    result0->SetLength(numSlots);
    result1->SetLength(numSlots);
    result2->SetLength(numSlots);
    result3->SetLength(numSlots);
    auto finalResult0 = result0->GetCKKSPackedValue();
    auto finalResult1 = result1->GetCKKSPackedValue();
    auto finalResult2 = result2->GetCKKSPackedValue();
    auto finalResult3 = result3->GetCKKSPackedValue();
    int wrong         = 0;
    // std::cout << "x0: " << x0[0].real() << ", x1: " << x1[0].real() << " y0: " << y0[0].real()
    //           << ", y1: " << y1[0].real() << ", cres0: " << finalResult0[0].real()
    //           << ", cres1: " << finalResult1[0].real() << ", cres2: " << finalResult2[0].real()
    //           << ", cres3: " << finalResult3[0].real() << std::endl;
    // for (int i = 0; i < numSlots; i++) {
    //     auto M   = int(x0[i].real() - y0[i].real()) + int(x1[i].real() - y1[i].real()) * pow(2, k);
    //     auto res = int(finalResult0[i].real()) + int(finalResult1[i].real()) * pow(2, k) +
    //                int(finalResult2[i].real()) * pow(2, 2 * k) + int(finalResult3[i].real()) * pow(2, 3 * k);

    //     if (M * M != res) {
    //         std::cout << "M*M : " << M * M << ", res:" << res << std::endl;
    //         std::cout << "x0: " << x0[i].real() << ", x1: " << x1[i].real() << " y0: " << y0[i].real()
    //                   << ", y1: " << y1[i].real() << ", cres0: " << finalResult0[i].real()
    //                   << ", cres1: " << finalResult1[i].real() << ", cres2: " << finalResult2[i].real()
    //                   << ", cres3: " << finalResult3[i].real() << std::endl;
    //         wrong++;
    //     }
    // }
    std::cout << "x0: " << x0[0] << ", x1: " << x1[0] << " y0: " << y0[0] << ", y1: " << y1[0]
              << ", cres0: " << finalResult0[0].real() << ", cres1: " << finalResult1[0].real()
              << ", cres2: " << finalResult2[0].real() << ", cres3: " << finalResult3[0].real() << std::endl;
    for (int i = 0; i < numSlots; i++) {
        auto M0  = int(x0[i]) + int(x1[i]) * pow(2, k);
        auto M1  = int(y0[i]) + int(y1[i]) * pow(2, k);
        auto M   = M0 - M1;
        auto res = int(finalResult0[i].real()) + int(finalResult1[i].real()) * pow(2, k) +
                   int(finalResult2[i].real()) * pow(2, 2 * k) + int(finalResult3[i].real()) * pow(2, 3 * k);

        if (M * M != res) {
            std::cout << "M*M : " << M * M << ", res:" << res << std::endl;
            std::cout << "x0: " << x0[i] << ", x1: " << x1[i] << " y0: " << y0[i] << ", y1: " << y1[i]
                      << ", cres0: " << finalResult0[i].real() << ", cres1: " << finalResult1[i].real()
                      << ", cres2: " << finalResult2[i].real() << ", cres3: " << finalResult3[i].real() << std::endl;
            wrong++;
        }
    }
    std::cout << "level:" << depth - cres[0]->GetLevel() << std::endl;
    std::cout << "error num:" << wrong << std::endl;
    double totaltime = 0;
    totaltime += double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / double(numSlots));
    std::cout << "shift time:" << totaltime << "ms" << std::endl;
}

void LazyCarry_dist(int bits, int p) {
    CCParams<CryptoContextCKKSRNS> parameters;

    std::cout << "SIFT1M LazyCarry_dist test" << std::endl << std::endl;
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
    int plain_bits                = 8;
    int Bg                        = 4;
    int w                         = plain_bits / Bg;  // number of segments
    std::cout << "w =" << w << std::endl;
    /*usint depth = levelBudget[0] + levelBudget[1] + 12 + std::log2(p);*/
    usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 4 + 6 + 3 + 8 + 4;
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    cc->Enable(FBTS);

    usint ringDim = cc->GetRingDimension();
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;
    std::cout << "8 bit shift " << std::endl;
    std::cout << "FBS Precision " << bits << " bits" << std::endl;

    std::cout << "depth " << depth << std::endl << std::endl;
    // Precomputations for bootstrapping
    int numSlots = ringDim / 2;

    cc->EvalFuncBootstrapSetup(levelBudget, bsgsDim, numSlots, bits);

    // Key Generation
    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);
    
    // auto f_mod16 = [](double x) -> double {
    //     return (int)x % 16;
    // };
    auto f = [](double x) -> double {
        return x - 16.0 * std::floor((x + 8.0) / 16.0);
    };
    std::vector<std::complex<double>> hermite_coefficients;
    std::string filename = "testdata/coeffcient/" + std::to_string(bits) + "carry_mod_16.bin";

    std::ifstream check_file(filename);

    if (check_file.good()) {
        check_file.close();
        std::cout << "检测到缓存文件，正在加载..." << std::endl;
        try {
            hermite_coefficients = LoadCoefficients(filename);
        }
        catch (const std::exception& e) {
            std::cerr << "加载失败: " << e.what() << "，将重新计算。" << std::endl;
            
            hermite_coefficients = hermiteInterpOrder1Coeffshift([f](double x) -> double { return 0.5 * f(x); }, p);
        }
    }
    else {
        
        std::cout << "未检测到缓存文件，正在计算..." << std::endl;
        check_file.close();
        hermite_coefficients = hermiteInterpOrder1Coeffshift([f](double x) -> double { return 0.5 * f(x); }, p);
        SaveCoefficients(filename, hermite_coefficients);
    }

    int vec_size = 2048;
    int degree;
    std::vector<std::vector<double>> read_sift, db_feature;

    std::string filename_index = "./index/parsed/sift1m_centroids.csv";
    read_sift                  = read_vectors_from_csv(filename_index, vec_size);
    // std::string filename = "../sift_base_100k.csv";
    // read_sift = read_vectors_from_csv(filename, vec_size);

    // std::cout << "read_sift " << read_sift[0] << std::endl;
    degree     = 128;
    db_feature = transpose(read_sift);

    std::cout << "degree =" << degree << std::endl;
    // std::cout << "ciphertext num = " << ct_num << std::endl;
    // double S0_max, S1_max, S1_min, S2_max, dismax;
    // S0_max = 0;
    // S1_max = 0;
    // S1_min = 0;
    // S2_max = 0;
    // dismax = 0;
    // degree * w * vecsize
    std::vector<std::vector<Ciphertext<DCRTPoly>>> edb_feature(degree), enc_input(degree);
    std::vector<std::vector<std::vector<double>>> discrete_db, discrete_query;
    for (int i = 0; i < degree; i++) {
        discrete_db.push_back(quantization(db_feature[i], plain_bits, Bg));

        for (int j = 0; j < vec_size; j++) {
            double quant_res = 0;
            for (int k = 0; k < 2; k++) {
                quant_res += discrete_db[i][k][j] * pow(2, Bg * k);
            }
            if (quant_res != db_feature[i][j])
                std::cout << "db_feature: " << db_feature[i][j] << ", quant_res: " << quant_res << " ->  "
                          << discrete_db[i][0][j] << ", " << discrete_db[i][1][j] << std::endl;
        }
        for (int j = 0; j < w; j++) {
            Plaintext plain_feat = cc->MakeCKKSPackedPlaintext(discrete_db[i][j], 1);
            edb_feature[i].push_back(cc->Encrypt(keyPair.publicKey, plain_feat));
        }
    }
    std::vector<double> x_rshift_16(numSlots, double(1.0 / 16));

    Plaintext ptxt_rshift_16            = cc->MakeCKKSPackedPlaintext(x_rshift_16);
    Ciphertext<DCRTPoly> ctxt_rshift_16 = cc->Encrypt(keyPair.publicKey, ptxt_rshift_16);
    std::vector<std::vector<Ciphertext<DCRTPoly>>> enc_dis;
    int query_num          = 2;
    std::string query_file = "../sift_query.fvecs";
    std::cout << "query_num:" << query_num << std::endl;
    double totaltime = 0.0;
    int seg_num      = 5;
    std::chrono::system_clock::time_point start, end;
    for (int query_id = 0; query_id < query_num; query_id++) {
        std::cout << "---------------- ct: " << query_id << " ----------------" << std::endl;

        discrete_query.clear();
        for (auto& vec : enc_input) {
            vec.clear();
        }

        std::vector<float> query_vec_f = read_fvec(query_file, query_id);  
        std::vector<double> query_vec(query_vec_f.begin(), query_vec_f.end());
        std::vector<double> vec_i(vec_size);

        double true_S0            = 0;
        double true_S1            = 0;
        double true_S2            = 0;
        double exact_integer_dist = 0;

        for (int i = 0; i < degree; i++) {
            std::fill(vec_i.begin(), vec_i.end(), query_vec[i]);

            
            discrete_query.push_back(quantization(vec_i, plain_bits, Bg));

            for (int j = 0; j < w; j++) {
                Plaintext plain_vec = cc->MakeCKKSPackedPlaintext(discrete_query[i][j], 1);
                enc_input[i].push_back(cc->Encrypt(keyPair.publicKey, plain_vec));
            }

            double db_low  = discrete_db[i][0][0];
            double db_high = discrete_db[i][1][0];

            double q_low  = discrete_query[i][0][0];
            double q_high = discrete_query[i][1][0];

            double diff_low  = q_low - db_low;
            double diff_high = q_high - db_high;

            true_S0 += diff_low * diff_low;
            true_S1 += 2 * diff_low * diff_high;
            true_S2 += diff_high * diff_high;
        }
        // S0_max             = S0_max > true_S0 ? S0_max : true_S0;
        // S1_max             = S1_max > true_S1 ? S1_max : true_S1;
        // S1_min             = S1_min < true_S1 ? S1_min : true_S1;
        // S2_max             = S2_max > true_S2 ? S2_max : true_S2;
        exact_integer_dist = true_S0 + true_S1 * 16.0 + true_S2 * 256.0;
        // dismax             = dismax > exact_integer_dist ? dismax : exact_integer_dist;
        // std::cout << "[TRUTH] Expected S0: " << true_S0 << std::endl;
        // std::cout << "[TRUTH] Expected S1: " << true_S1 << std::endl;
        // std::cout << "[TRUTH] Expected S2: " << true_S2 << std::endl;
        // std::cout << "[TRUTH] Exact Dist : " << exact_integer_dist << std::endl;

        start = std::chrono::system_clock::now();
        
        enc_dis = Homdist(enc_input, edb_feature);

        // std::cout << "Homdist finish  " << std::endl;
        // auto res = Homomorphic_lazy_carry(enc_dis, 6, p, Bg, depth, levelBudget, keyPair.secretKey, numSlots, vec_size,
        //                                   true_S0, true_S1, true_S2);
        auto res = Homomorphic_lazy_carry_iterate(enc_dis, ctxt_rshift_16, hermite_coefficients, seg_num, p, Bg, depth,
                                                  levelBudget, keyPair.secretKey, numSlots, vec_size, true_S0, true_S1,
                                                  true_S2);
        end      = std::chrono::system_clock::now();
        totaltime += double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
        std::vector<int64_t> final_segments(seg_num);
        for (int i = 0; i < seg_num; ++i) {
            Plaintext ptxt;
            cc->Decrypt(keyPair.secretKey, res[i], &ptxt);
            final_segments[i] = std::round(ptxt->GetCKKSPackedValue()[0].real());
        }

        double decrypted_dist = 0;
        for (int i = 0; i < seg_num; ++i) {
            decrypted_dist += final_segments[i] * std::pow(16.0, i);
        }

        std::cout << "\n---------------- Result Analysis ----------------" << std::endl;
        std::cout << "Exact Integer Dist  : " << (long long)exact_integer_dist << std::endl;
        std::cout << "Decrypted Dist      : " << (long long)decrypted_dist << std::endl;
        std::cout << "Diff                : " << std::abs(exact_integer_dist - decrypted_dist) << std::endl;

        std::cout << "Segments (Low->High): ";
        for (auto v : final_segments)
            std::cout << v << " ";
        std::cout << std::endl;
    }

    // std::cout << "S0 max: " << S0_max << std::endl;
    // std::cout << "S1 max: " << S1_max << std::endl;
    // std::cout << "S2 max: " << S2_max << std::endl;
    // std::cout << "S1 min: " << S1_min << std::endl;
    // std::cout << "dismax: " << dismax << std::endl;
    std::cout << "carry time per query: " << totaltime / query_num << " ms" << std::endl;
    std::cout << "EvalFuncBootstrap amortize: " << totaltime / (query_num * numSlots) * vec_size << " ms" << std::endl;

    cc->ClearEvalSumKeys();
    cc->ClearEvalMultKeys();
    cc->ClearEvalAutomorphismKeys();
    CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();

    return;
}

// -------------------------------------------------------------------------



// -------------------------------------------------------------------------
void CarrymulExampleDynamic(int total_bits, int segment_bits, int p) {
    std::cout << "\n===============================================================" << std::endl;
    std::cout << "Running Dynamic Carry Mul Test" << std::endl;
    std::cout << "Total Precision: " << total_bits << " bits" << std::endl;
    std::cout << "Segment Size:    " << segment_bits << " bits" << std::endl;
    std::cout << "Res Precision: " << total_bits * 2 << " bits" << std::endl;

    if (total_bits <= 0 || (segment_bits != 4 && segment_bits != 8) || p <= 0) {
        std::cerr << "Error: total_bits and p must be positive, and segment_bits must be 4 or 8." << std::endl;
        return;
    }

    int num_chunks = total_bits / segment_bits;
    if (total_bits % segment_bits != 0) {
        std::cerr << "Error: Total bits must be divisible by segment bits." << std::endl;
        return;
    }
    std::cout << "Number of Chunks: " << num_chunks << std::endl;

    // 1. Context Setup
    // ----------------
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecretKeyDist(SPARSE_TERNARY);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 16);  
    parameters.SetNumLargeDigits(3);
    parameters.SetKeySwitchTechnique(HYBRID);

    parameters.SetScalingModSize(59);
    parameters.SetFirstModSize(60);
    parameters.SetScalingTechnique(FIXEDMANUAL);

    std::vector<uint32_t> levelBudget = {4, 2};
    std::vector<uint32_t> bsgsDim     = {0, 0};

    
    
    usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 4 + 6 + 3 + 8 + 4;

    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);
    cc->Enable(FBTS);

    std::cout << "CKKS Ring Dimension: " << cc->GetRingDimension() << std::endl;
    std::cout << "Computed Depth:      " << depth << std::endl;

    // 2. Precomputations
    // ------------------
    int numSlots = cc->GetRingDimension() / 2;
    int bootstrap_bits = static_cast<int>(std::ceil(std::log2(static_cast<double>(p))));
    cc->EvalFuncBootstrapSetup(levelBudget, bsgsDim, numSlots, bootstrap_bits);

    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);

    // 3. Data Generation
    // ------------------
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    long long max_val = std::pow(2, segment_bits) - 1;
    std::uniform_int_distribution<int> message(0, max_val);
    std::vector<std::vector<double>> inputs(num_chunks);

    for (int k = 0; k < num_chunks; k++) {
        inputs[k].reserve(numSlots);
        for (int i = 0; i < numSlots; i++) {
            inputs[k].push_back((double)message(engine));
        }
    }

    // 4. Encryption
    // -------------
    std::vector<Ciphertext<DCRTPoly>> ciphertexts(num_chunks);
    for (int k = 0; k < num_chunks; k++) {
        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(inputs[k]);
        ciphertexts[k] = cc->Encrypt(keyPair.publicKey, ptxt);
    }

    // 5. Execution
    // ------------
    std::cout << "Executing HomomorphicCarryMult..." << std::endl;

    std::vector<Ciphertext<DCRTPoly>> cres;
    auto f_mod = [](double x) -> double {
        return int(x) % 16;
    };

    std::vector<std::complex<double>> coe_mod =
        hermiteInterpOrder1Coeff([f_mod](double x) -> double { return 0.5 * f_mod(x); }, p);
    auto start = std::chrono::system_clock::now();
    if (segment_bits == 8) {
        
        
        cres = HomomorphicCarryMult8BitDynamic_Debug(ciphertexts, ciphertexts, inputs, inputs, p, depth, levelBudget,
                                                     keyPair.secretKey, coe_mod);
    }
    else if (segment_bits == 4) {
        cres = HomomorphicCarryMult4BitDynamic(ciphertexts, ciphertexts, p, depth, levelBudget, coe_mod);
    }
    else {
        std::cerr << "Error: segment_bits currently supports only 4 or 8." << std::endl;
        return;
    }

    auto end = std::chrono::system_clock::now();

    std::cout << "cres level: " << depth - cres[0]->GetLevel() << std::endl;
    // 6. Verification & Decryption
    // ----------------------------
    
    size_t result_count = cres.size();
    std::vector<std::vector<double>> decrypted_results(result_count);

    for (size_t k = 0; k < result_count; k++) {
        Plaintext result_ptxt;
        cc->Decrypt(keyPair.secretKey, cres[k], &result_ptxt);
        result_ptxt->SetLength(numSlots);
        std::vector<std::complex<double>> complex_res = result_ptxt->GetCKKSPackedValue();
        for (const auto& val : complex_res) {
            decrypted_results[k].push_back(val.real());
        }
    }

    int wrong      = 0;
    double epsilon = 1e-2;  

    for (int i = 0; i < numSlots; i++) {
        
        double M = 0;
        for (int k = 0; k < num_chunks; k++) {
            M += inputs[k][i] * std::pow(2, k * segment_bits);
        }

        
        double Res = 0;
        for (size_t k = 0; k < result_count; k++) {
            Res += decrypted_results[k][i] * std::pow(2, k * segment_bits);
        }

        
        double Expected = M * M;

        
        double diff    = std::abs(Expected - Res);
        double rel_err = (Expected != 0) ? diff / Expected : diff;

        if (diff > epsilon && rel_err > 0.001) {
            if (wrong < 10) {
                std::cout << "Slot " << i << " Fail: "
                          << "Exp=" << std::fixed << std::setprecision(0) << Expected << ", Got=" << Res
                          << " (Diff=" << diff << ")" << std::endl;
            }
            wrong++;
        }
    }

    std::cout << "Total Errors: " << wrong << " / " << numSlots << std::endl;

    double totaltime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Total Time: " << totaltime << " ms" << std::endl;
    std::cout << "Avg Time per Slot: " << totaltime / numSlots << " ms" << std::endl;
}

void CarryaddExampleDynamic(int total_bits, int segment_bits, int p) {
    auto total_start = std::chrono::steady_clock::now();

    std::cout << "\n===============================================================" << std::endl;
    std::cout << "Running Dynamic Carry Add Test" << std::endl;
    std::cout << "Total Precision: " << total_bits << " bits" << std::endl;
    std::cout << "Segment Size:    " << segment_bits << " bits" << std::endl;

    if (total_bits <= 0 || segment_bits != 4 || total_bits % segment_bits != 0 || p <= 0) {
        std::cerr << "Error: addition test requires positive, 4-bit-aligned precision." << std::endl;
        return;
    }

    int num_chunks = total_bits / segment_bits;
    std::cout << "Number of Chunks: " << num_chunks << std::endl;

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecretKeyDist(SPARSE_TERNARY);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 16);
    parameters.SetNumLargeDigits(3);
    parameters.SetKeySwitchTechnique(HYBRID);
    parameters.SetScalingModSize(59);
    parameters.SetFirstModSize(60);
    parameters.SetScalingTechnique(FIXEDMANUAL);

    std::vector<uint32_t> levelBudget = {4, 2};
    std::vector<uint32_t> bsgsDim     = {0, 0};

    // Keep the same fixed depth as CarrymulExample32bits. Carry propagation
    // uses an encrypted 1/16 and functional bootstrapping at every segment,
    // so the context depth does not grow with the number of chunks.
    usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 4 + 6 + 3 + 8 + 4;
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);
    cc->Enable(FBTS);

    int numSlots       = cc->GetRingDimension() / 2;
    int bootstrap_bits = static_cast<int>(std::ceil(std::log2(static_cast<double>(p))));
    std::cout << "CKKS Ring Dimension: " << cc->GetRingDimension() << std::endl;
    std::cout << "Computed Depth:      " << depth << std::endl;
    std::cout << "Packed Slots:        " << numSlots << std::endl;

    cc->EvalFuncBootstrapSetup(levelBudget, bsgsDim, numSlots, bootstrap_bits);

    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);

    int base = 1 << segment_bits;
    std::mt19937 engine(0xC0FFEEu + static_cast<uint32_t>(total_bits));
    std::uniform_int_distribution<int> message(0, base - 1);
    std::vector<std::vector<double>> inputsA(num_chunks, std::vector<double>(numSlots));
    std::vector<std::vector<double>> inputsB(num_chunks, std::vector<double>(numSlots));
    std::vector<Ciphertext<DCRTPoly>> chunksA(num_chunks), chunksB(num_chunks);

    for (int k = 0; k < num_chunks; ++k) {
        for (int i = 0; i < numSlots; ++i) {
            inputsA[k][i] = message(engine);
            inputsB[k][i] = message(engine);
        }
        chunksA[k] = cc->Encrypt(keyPair.publicKey, cc->MakeCKKSPackedPlaintext(inputsA[k]));
        chunksB[k] = cc->Encrypt(keyPair.publicKey, cc->MakeCKKSPackedPlaintext(inputsB[k]));
    }

    std::vector<double> inverse_base(numSlots, 1.0 / base);
    auto ctxt_rshift_16 = cc->Encrypt(keyPair.publicKey, cc->MakeCKKSPackedPlaintext(inverse_base));

    auto f_mod16 = [](double x) -> double {
        return static_cast<int>(x) % 16;
    };
    std::vector<std::complex<double>> coe_mod16 =
        hermiteInterpOrder1Coeff([f_mod16](double x) -> double { return 0.5 * f_mod16(x); }, p);

    std::cout << "Executing HomomorphicCarryAdd4BitDynamic..." << std::endl;
    auto eval_start = std::chrono::steady_clock::now();
    auto cres =
        HomomorphicCarryAdd4BitDynamic(chunksA, chunksB, p, depth, levelBudget, coe_mod16, ctxt_rshift_16);
    auto eval_end = std::chrono::steady_clock::now();

    std::vector<std::vector<std::complex<double>>> decrypted(cres.size());
    for (size_t k = 0; k < cres.size(); ++k) {
        Plaintext result;
        cc->Decrypt(keyPair.secretKey, cres[k], &result);
        result->SetLength(numSlots);
        decrypted[k] = result->GetCKKSPackedValue();
    }

    size_t wrong_chunks = 0;
    size_t wrong_slots  = 0;
    int printed         = 0;
    for (int slot = 0; slot < numSlots; ++slot) {
        int carry       = 0;
        bool slot_wrong = false;
        for (int k = 0; k < num_chunks; ++k) {
            int sum      = static_cast<int>(inputsA[k][slot]) + static_cast<int>(inputsB[k][slot]) + carry;
            int expected = sum % base;
            carry        = sum / base;
            double got   = decrypted[k][slot].real();
            if (std::abs(got - expected) > 0.5) {
                ++wrong_chunks;
                slot_wrong = true;
                if (printed < 10) {
                    std::cout << "Slot " << slot << ", chunk " << k << " Fail: Exp=" << expected
                              << ", Got=" << got << std::endl;
                    ++printed;
                }
            }
        }

        double got_final_carry = decrypted[num_chunks][slot].real();
        if (std::abs(got_final_carry - carry) > 0.5) {
            ++wrong_chunks;
            slot_wrong = true;
            if (printed < 10) {
                std::cout << "Slot " << slot << ", final carry Fail: Exp=" << carry
                          << ", Got=" << got_final_carry << std::endl;
                ++printed;
            }
        }
        if (slot_wrong)
            ++wrong_slots;
    }

    auto total_end = std::chrono::steady_clock::now();
    double eval_ms = std::chrono::duration_cast<std::chrono::milliseconds>(eval_end - eval_start).count();
    double total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - total_start).count();

    std::cout << "Result Chunks: " << cres.size() << std::endl;
    std::cout << "Chunk Errors: " << wrong_chunks << " / "
              << static_cast<size_t>(numSlots) * static_cast<size_t>(num_chunks + 1) << std::endl;
    std::cout << "Slot Errors: " << wrong_slots << " / " << numSlots << std::endl;
    std::cout << "Carry Add Time: " << eval_ms << " ms" << std::endl;
    std::cout << "Avg Add Time per Slot: " << eval_ms / numSlots << " ms" << std::endl;
    std::cout << "End-to-End Time: " << total_ms << " ms" << std::endl;
}

int main(int argc, char* argv[]) {
    int total_bits   = (argc > 1) ? std::stoi(argv[1]) : 32;
    int segment_bits = (argc > 2) ? std::stoi(argv[2]) : 4;
    std::string operation = (argc > 3) ? argv[3] : "mul";

    if (total_bits <= 0 || (segment_bits != 4 && segment_bits != 8) || total_bits % segment_bits != 0 ||
        (operation != "add" && operation != "mul") || (operation == "add" && segment_bits != 4)) {
        std::cerr << "Usage: " << argv[0] << " [total_bits] [segment_bits: 4|8] [add|mul]" << std::endl;
        return 1;
    }

    // A 4-bit carry is refreshed as separate hexadecimal digits before it is
    // reconstructed, so every functional bootstrap stays in the p=256 domain.
    int num_chunks = total_bits / segment_bits;
    int base       = 1 << segment_bits;
    int bootstrap_bits;
    if (operation == "add") {
        long long max_bootstrap_input = 2LL * (base - 1) + 1;
        bootstrap_bits = static_cast<int>(std::ceil(std::log2(static_cast<double>(max_bootstrap_input) + 1.0)));
    }
    else if (segment_bits == 4) {
        bootstrap_bits = 8;
        if (num_chunks > 64) {
            std::cerr << "Error: segmented 4-bit carry currently supports at most 256-bit inputs." << std::endl;
            return 1;
        }
    }
    else {
        long long max_bootstrap_input = static_cast<long long>(num_chunks) * base * (base - 1);
        bootstrap_bits = static_cast<int>(std::ceil(std::log2(static_cast<double>(max_bootstrap_input) + 1.0)));
    }
    if (bootstrap_bits >= 31) {
        std::cerr << "Error: the required functional-bootstrap domain does not fit in int." << std::endl;
        return 1;
    }
    int p = 1 << bootstrap_bits;

    std::cout << "Functional-bootstrap domain: 2^" << bootstrap_bits << " = " << p << std::endl;

    // CarryExample(bits, p, k);
    // CarrymulExample16bits(bits, p, k);
    // CarrymulExample32bits(bits, p, k);
    // LazyCarry_dist(bits, p);

    if (operation == "add")
        CarryaddExampleDynamic(total_bits, segment_bits, p);
    else
        CarrymulExampleDynamic(total_bits, segment_bits, p);
}
