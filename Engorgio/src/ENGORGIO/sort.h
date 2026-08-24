#include "openfhe.h"
#include <chrono>
#include "math/chebyshev.h"
#include <stdexcept>
namespace openfhe
{
    using namespace lbcrypto;

    // std::vector<double> coeff1 = {1.5, -0.5}; // 1.5 x  - 0.5 x ^ 3
    // std::vector<double> coeff3 = {2.1875, -2.1875, 1.3125, -0.3125};
    // std::vector<double> coeff5 = {2.4609375 / 2, -3.28125 / 2, 2.953125 / 2, -1.40625 / 2, 0.2734375 / 2};
    // -0.2095 x ^ 15 + 1.692 x ^ 13 + -5.999 x ^ 11 + 12.22 x ^ 9 + -15.71 x ^ 7 + 13.2 x ^ 5 + -7.332 x ^ 3 + 3.142 x

    std::vector<std::vector<std::vector<double>>> bitonic_test_plain(std::vector<double> vec_unsort,
                                                                     std::vector<std::vector<double>> &plain_sort,
                                                                     std::uint32_t pack_slots);

    std::vector<std::vector<std::vector<double>>> bitonic_test_plain_full_table_scan(std::vector<double> vec_unsort,
                                                                                     std::vector<std::vector<double>> &plain_sort,
                                                                                     std::uint32_t pack_slots, int test_stage);
    void bitonic_sort(int plain_bits, int num_slots);

    double bitonic_sort_query(int plain_bits, int num_slots, int length, double precision, usint multDepth, Ciphertext<lbcrypto::DCRTPoly> &ciphertext_unsort,
                              lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privatekey, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &sync_matrix);

    void bitonic_sort_small(int plain_bits, int num_slots);

    std::vector<std::vector<std::vector<double>>> bitonictopk_test_plain(std::vector<double> vec_unsort,
                                                                         std::vector<std::vector<double>> &plain_sort,
                                                                         std::uint32_t pack_slots, std::uint32_t k);
    double topk_sort_test(int plain_bits, int num_slots, int k);

    std::vector<double> getKSmallestBySort(const std::vector<double> &nums, int N, int k);

    double topk_sort(int plain_bits, int num_slots, int length, int k, double precision,
                     std::vector<double> topk_plain, usint multDepth, Ciphertext<lbcrypto::DCRTPoly> &ciphertext_unsort,
                     lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privatekey, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &sync_matrix);

    double bitonic_sort_full_table_scan(int plain_bits, int length);

    std::vector<double> bitonic_sort_modular(int plain_bits, int blocks, int num_slots);

    std::vector<double> bitonic_sort_modular_query(int plain_bits, int blocks, int num_slots);

    double sync_test_big(int num_column, int plain_bits, int num_slots);

    double sync_test_small(int plain_bits, int num_slots);
}