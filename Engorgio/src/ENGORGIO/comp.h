#include "openfhe.h"
#include <chrono>

namespace openfhe
{
    using namespace lbcrypto;
    Ciphertext<DCRTPoly> MultByInteger(const ConstCiphertext<DCRTPoly> &ciphertext, const int64_t constant);
    template <typename T>
    inline uint64_t CeilLog2(T x);

    void print_moduli_chain(const DCRTPoly &poly);
    extern std::vector<double> coeff1; // 1.5 x  - 0.5 x ^ 3
    extern std::vector<double> coeff3;
    extern std::vector<double> coeff5;
    extern std::vector<double> g;

    void EvalPower(std::vector<double> coefficients, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &power_basis,
                   Ciphertext<lbcrypto::DCRTPoly> &result);

    void poly_evaluate_power(Ciphertext<lbcrypto::DCRTPoly> &result, Ciphertext<lbcrypto::DCRTPoly> &x,
                             std::vector<double> &coefficients);

    void Homround(Ciphertext<lbcrypto::DCRTPoly> &cipher);

    void comp_greater_than(Ciphertext<lbcrypto::DCRTPoly> &ct1, Ciphertext<lbcrypto::DCRTPoly> &ct2, double precision,
                           int polyDegree, Ciphertext<lbcrypto::DCRTPoly> &res,
                           lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey);

    double CalculateApproximationError(const std::vector<std::complex<double>> &result,
                                       const std::vector<double> &expectedResult);

    int itboot(Ciphertext<lbcrypto::DCRTPoly> &res, lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey);

    void comp_equal(Ciphertext<lbcrypto::DCRTPoly> &ct1, Ciphertext<lbcrypto::DCRTPoly> &ct2, double precision,
                    int polyDegree, Ciphertext<lbcrypto::DCRTPoly> &res);

    void comp_greater_than_modular(std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ciphertext_a, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ciphertext_b, double precision,
                                   int polyDegree, Ciphertext<lbcrypto::DCRTPoly> &res,
                                   lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey);

    void comp_equal_modular(std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ciphertext_a, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ciphertext_b, double precision,
                            int polyDegree, Ciphertext<lbcrypto::DCRTPoly> &res);

    void comp_partial(Ciphertext<lbcrypto::DCRTPoly> &ct1, Ciphertext<lbcrypto::DCRTPoly> &ct2, double precision,
                      std::vector<double> &coefficients, Ciphertext<lbcrypto::DCRTPoly> &res,
                      lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey);

    void comp_partial_modular(std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ciphertext_a,
                              std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ciphertext_b, double precision,
                              std::vector<double> &coefficients, Ciphertext<lbcrypto::DCRTPoly> &comp_res,
                              lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey);

    void EvalSignExample(std::uint32_t plain_bits);

    void EvalequalExample(std::uint32_t plain_bits);

    void Eval_modular_greater_than(std::uint32_t plain_bits, std::uint32_t block);

    double Eval_modular_Strictly_greater_than(std::uint32_t plain_bits, std::uint32_t block);

    double Eval_modular_Strictly_equal(std::uint32_t plain_bits, std::uint32_t block);

    void Eval_modular_greater_test(std::uint32_t plain_bits, std::uint32_t block);

    double Eval_modular_Strictly_greater_test(std::uint32_t plain_bits, std::uint32_t block);

    double Eval_modular_Strictly_equal_test(std::uint32_t plain_bits, std::uint32_t block);
}