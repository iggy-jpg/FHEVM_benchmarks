#include "openfhe.h"
#include <chrono>
#include "comp.h"
#include "quant.h"
#include "orderapp.h"
#include "ordergen.h"
#include "sort.h"
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
namespace openfhe
{
    void output(Ciphertext<lbcrypto::DCRTPoly> &res, int slot, lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey);

    void inverse(Ciphertext<lbcrypto::DCRTPoly> &res);

    double error_estimate(std::vector<double> plain, Ciphertext<lbcrypto::DCRTPoly> &res,
                          lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey, size_t num_slots);

    double error_estimate_unbounded_modular(std::vector<double> plain, std::vector<std::vector<Ciphertext<lbcrypto::DCRTPoly>>> &res,
                                            lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey, int num_ct, int num_slots, int length);

    double error_estimate_modular(std::vector<double> plain, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &res,
                                  lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey, size_t num_slots);

    double error_estimate_unbounded(std::vector<double> plain, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &res,
                                    lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey);

    std::vector<int> dot_product(std::vector<int> &A, std::vector<int> &B);

    void cross_rotate(int N, int l, int T);

    std::vector<Ciphertext<lbcrypto::DCRTPoly>> gen_rotate(int N, int l, int T,
                                                           std::vector<Ciphertext<lbcrypto::DCRTPoly>> &big_vec,
                                                           std::vector<double> &test,
                                                           lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey);

    void test_gen_rotate(int plain_bits, int l, int num_slots);

    double Eval_Agg_4(int plain_bits, int num_slots);

    double Eval_Agg_1(int plain_bits, int num_slots);

    std::vector<std::vector<double>> read_vectors_from_csv(const std::string &filename, int N);

    double Euclid_distance(int plain_bits);

    Ciphertext<lbcrypto::DCRTPoly> homdis(
        Ciphertext<lbcrypto::DCRTPoly> &ct0, Ciphertext<lbcrypto::DCRTPoly> &ct1);

    std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> &matrix);
}