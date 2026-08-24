#include "openfhe.h"
#include <chrono>
#include "math/chebyshev.h"
#include <stdexcept>
namespace openfhe
{
    void bitonic_comp_unbounded_modular(int stage, int part, int slots, int length, std::vector<std::vector<Ciphertext<lbcrypto::DCRTPoly>>> &ct,
                                        Ciphertext<lbcrypto::DCRTPoly> &res, double precision, std::vector<double> &coefficients,
                                        lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey, int blocks);
    void bitonic_comp_topk(int stage, int slots, int topk_i, int k, Ciphertext<lbcrypto::DCRTPoly> &ct,
                           Ciphertext<lbcrypto::DCRTPoly> &res, double precision, std::vector<double> &coefficients,
                           lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey);

    void bitonic_comp(int stage, int part, int slots, Ciphertext<lbcrypto::DCRTPoly> &ct,
                      Ciphertext<lbcrypto::DCRTPoly> &res, double precision, std::vector<double> &coefficients,
                      lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey);

    void bitonic_comp_modular(int stage, int part, int slots, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ct,
                              Ciphertext<lbcrypto::DCRTPoly> &res, double precision, std::vector<double> &coefficients,
                              lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey);

    void bitonic_comp_unbounded(int stage, int part, int slots, int length, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ct,
                                Ciphertext<lbcrypto::DCRTPoly> &res, double precision, std::vector<double> &coefficients,
                                lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey);
}