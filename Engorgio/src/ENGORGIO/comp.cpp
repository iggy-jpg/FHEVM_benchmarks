#include "utils.h"
#include <chrono>
#include "cryptocontext.h"

#include "key/privatekey.h"
#include "key/publickey.h"
#include "math/chebyshev.h"
#include "schemerns/rns-scheme.h"
#include "scheme/ckksrns/ckksrns-cryptoparameters.h"
namespace openfhe
{
    using namespace lbcrypto;
    std::vector<double> coeff1 = {1.5, -0.5}; // 1.5 x  - 0.5 x ^ 3
    std::vector<double> coeff3 = {2.1875, -2.1875, 1.3125, -0.3125};
    std::vector<double> coeff5 = {2.4609375 / 2, -3.28125 / 2, 2.953125 / 2, -1.40625 / 2, 0.2734375 / 2};
    std::vector<double> g = {double(4589 / 1024), double(-16577 / 1024), double(25614 / 1024), double(-12860 / 1024)};
    template <typename T>
    inline uint64_t CeilLog2(T x)
    {
        return static_cast<uint64_t>(std::ceil(std::log2(x)));
    }

    Ciphertext<DCRTPoly> MultByInteger(const ConstCiphertext<DCRTPoly> &ciphertext, const int64_t constant)
    {
        auto result = ciphertext->Clone();
        for (auto &c : result->GetElements())
            c *= static_cast<typename DCRTPoly::Integer>(constant);
        return result;
    }

    double CalculateApproximationError(const std::vector<std::complex<double>> &result,
                                       const std::vector<double> &expectedResult)
    {
        // using the infinity norm
        double maxError = 0;
        for (size_t i = 0; i < result.size(); ++i)
        {
            double error = std::abs(result[i].real() - expectedResult[i]);
            if (maxError < error)
                maxError = error;
        }

        return std::abs(std::log2(maxError));
    }

    void print_moduli_chain(const DCRTPoly &poly)
    {
        int num_primes = poly.GetNumOfElements();
        double total_bit_len = 0.0;
        for (int i = 0; i < num_primes; i++)
        {
            auto qi = poly.GetParams()->GetParams()[i]->GetModulus();
            std::cout << "q_" << i << ": " << qi.ConvertToInt() << ",  log q_" << i << ": "
                      << log(qi.ConvertToInt()) / log(2) << std::endl;
            total_bit_len += log(qi.ConvertToInt()) / log(2);
        }
        std::cout << "Total bit length: " << total_bit_len << std::endl;
    }

    void EvalPower(std::vector<double> coefficients, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &power_basis,
                   Ciphertext<lbcrypto::DCRTPoly> &result)
    {
        auto cc = power_basis[0]->GetCryptoContext();
        if (coefficients.size() == 1)
        {
            if (std::fabs(std::round(coefficients[0] * pow(2, 50))) > 1.)
            {
                result = power_basis[0];
                result = cc->EvalMult(result, coefficients[0]);
                return;
            }
            else
            {
                return;
            }
        }
        std::vector<double> quotient, remainder;
        uint64_t degree = 2 * coefficients.size() - 1;
        uint64_t m = CeilLog2(degree + 1);
        remainder.resize((1 << (m - 1)) / 2);
        quotient.resize(coefficients.size() - remainder.size());

        for (size_t i = 0; i < remainder.size(); i++)
        {
            remainder[i] = coefficients[i];
        }
        for (size_t i = 0; i < quotient.size(); i++)
        {
            quotient[i] = coefficients[i + remainder.size()];
        }
        Ciphertext<lbcrypto::DCRTPoly> cipher_quotient, cipher_remainder;
        EvalPower(quotient, power_basis, cipher_quotient);
        EvalPower(remainder, power_basis, cipher_remainder);
        result = cc->EvalMult(cipher_quotient, power_basis[m - 1]);
        result = cc->EvalAdd(result, cipher_remainder);
    }

    void poly_evaluate_power(Ciphertext<lbcrypto::DCRTPoly> &result, Ciphertext<lbcrypto::DCRTPoly> &x,
                             std::vector<double> &coefficients)
    {
        auto cc = x->GetCryptoContext();
        uint64_t degree = coefficients.size() * 2 - 1;
        uint64_t m = CeilLog2(degree + 1);

        std::vector<Ciphertext<lbcrypto::DCRTPoly>> power_basis(m);
        power_basis[0] = x;
        for (size_t i = 1; i < m; i++)
        {
            power_basis[i] = cc->EvalMult(power_basis[i - 1], power_basis[i - 1]);
        }

        EvalPower(coefficients, power_basis, result);
    }

    void Homround(Ciphertext<lbcrypto::DCRTPoly> &cipher)
    {
        auto cc = cipher->GetCryptoContext();
        cipher = MultByInteger(cipher, 2.0);
        cipher = cc->EvalAdd(cipher, -1.0);
        poly_evaluate_power(cipher, cipher, coeff1);
        poly_evaluate_power(cipher, cipher, coeff5);
        cipher = cc->EvalAdd(cipher, 0.5);
    }

    void comp_greater_than(Ciphertext<lbcrypto::DCRTPoly> &ct1, Ciphertext<lbcrypto::DCRTPoly> &ct2, double precision,
                           int polyDegree, Ciphertext<lbcrypto::DCRTPoly> &res,
                           lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        auto cc = ct1->GetCryptoContext();
        // ct1-ct2 >0 -> 1; ct1-ct2 ==0 -> 0.5; ct1-ct2 <0 -> 0
        auto ciphertext_sub = cc->EvalSub(ct1, ct2);
        auto ciphertext_sign = cc->EvalSign(ciphertext_sub, 3, -precision, precision, polyDegree);
        Homround(ciphertext_sign);
        // ct1-ct2 >0 -> 2; ct1-ct2 ==0 -> 1; ct1-ct2 <0 -> 0
        // comp(a==b)
        // ct2-ct1
        auto ciphertext_sub_neg = cc->EvalNegate(ciphertext_sign);
        auto result_sub_neg = cc->EvalAdd(ciphertext_sub_neg, 1.0);
        //(ct2-ct1)*(ct1-ct2)
        auto result_temp = cc->EvalMult(result_sub_neg, ciphertext_sign);
        auto result_equal = MultByInteger(result_temp, 2.0);

        // get equal ct1-ct2 ==0 -> 1, else 0
        // Homround(result_equal);

        res = cc->EvalSub(ciphertext_sign, result_equal);
        // int bootprecision = itboot(res, privateKey);
        // std::cout << "compres bootprecision : " << bootprecision << std::endl;
        // res = cc->EvalBootstrap(res, 2, bootprecision);
        // Homround(res);
    }
    void comp_equal(Ciphertext<lbcrypto::DCRTPoly> &ct1, Ciphertext<lbcrypto::DCRTPoly> &ct2, double precision,
                    int polyDegree, Ciphertext<lbcrypto::DCRTPoly> &res)
    {
        auto cc = ct1->GetCryptoContext();
        // ct1-ct2 >0 -> 1; ct1-ct2 ==0 -> 0.5; ct1-ct2 <0 -> 0
        auto ciphertext_sub = cc->EvalSub(ct1, ct2);
        auto ciphertext_sign = cc->EvalSign(ciphertext_sub, 3, -precision - 10, precision + 10, polyDegree);
        Homround(ciphertext_sign);
        // ct1-ct2 >0 -> 2; ct1-ct2 ==0 -> 1; ct1-ct2 <0 -> 0
        // auto ciphertext_greater_equal_2 = MultByInteger(ciphertext_sign, 2);
        // comp(a==b)
        // ct2-ct1
        auto ciphertext_sub_neg = cc->EvalNegate(ciphertext_sign);
        auto result_sub_neg = cc->EvalAdd(ciphertext_sub_neg, 1.0);
        //(ct2-ct1)*(ct1-ct2)
        auto result_equal = cc->EvalMult(result_sub_neg, ciphertext_sign);
        res = MultByInteger(result_equal, 4);
        // get equal ct1-ct2 ==0 -> 1, else 0
        // Homround(res);
    }
    void comp_greater_than_modular(std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ciphertext_a, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ciphertext_b, double precision,
                                   int polyDegree, Ciphertext<lbcrypto::DCRTPoly> &res,
                                   lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        auto cc = ciphertext_a[0]->GetCryptoContext();
        int block = ciphertext_a.size();
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_sub(block), result_sub(block), result_equal_comp(block),
            ciphertext_comp_greater(block);
        for (int i = block - 1; i >= 0; i--)
        {
            // comp(a>b) 1, 0.5, 0

            ciphertext_sub[i] = cc->EvalSub(ciphertext_a[i], ciphertext_b[i]);
            result_sub[i] = cc->EvalSign(ciphertext_sub[i], 3, -precision, precision, polyDegree);
            Homround(result_sub[i]);

            // comp(a==b) 1, 0

            auto ciphertext_sub_neg = cc->EvalNegate(result_sub[i]);
            auto result_sub_neg = cc->EvalAdd(ciphertext_sub_neg, 1.0);

            auto result_equal = cc->EvalMult(result_sub_neg, result_sub[i]);
            result_equal_comp[i] = MultByInteger(result_equal, 2.0);
            ciphertext_comp_greater[i] = cc->EvalSub(result_sub[i], result_equal_comp[i]);
            result_equal_comp[i] = MultByInteger(result_equal_comp[i], 2.0);
            // Homround(result_equal_comp[i]);
            // Homround(ciphertext_comp_greater[i]);
        }

        for (int i = 0; i < block; i++)
        {
            if (i == 0)
                res = ciphertext_comp_greater[i];
            else
            {
                auto ai_sub_bi = cc->EvalSub(res, ciphertext_comp_greater[i]);
                auto equal_mul_sub = cc->EvalMult(result_equal_comp[i], ai_sub_bi);
                res = cc->EvalAdd(equal_mul_sub, ciphertext_comp_greater[i]);
            }
        }
        Homround(res);
    }
    // check first boot precision
    int itboot(Ciphertext<lbcrypto::DCRTPoly> &res, lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        auto cc = res->GetCryptoContext();
        std::vector<double> error_vec;
        Plaintext plaintextDec;
        cc->Decrypt(privateKey, res, &plaintextDec);
        // std::cout << "Estimated precision: " << plaintextDec->GetLogPrecision() << std::endl;
        std::vector<double> expectedResult;
        auto dec = plaintextDec->GetCKKSPackedValue();
        for (int i = 0; i < int(dec.size()); i++)
            expectedResult.push_back(dec[i].real());
        auto boot_simple = cc->EvalBootstrap(res);
        Plaintext boot_1;
        cc->Decrypt(privateKey, boot_simple, &boot_1);
        int boot_precision = std::floor(CalculateApproximationError(boot_1->GetCKKSPackedValue(), expectedResult));

        return boot_precision;
    }

    // ct1==ct2?1:0
    void comp_equal_modular(std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ciphertext_a, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ciphertext_b, double precision,
                            int polyDegree, Ciphertext<lbcrypto::DCRTPoly> &res)
    {
        auto cc = ciphertext_a[0]->GetCryptoContext();
        int block = ciphertext_a.size();
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_sub(block), result_sub(block), result_equal_comp(block),
            ciphertext_comp_greater(block);
        for (int i = block - 1; i >= 0; i--)
        {
            // comp(a>b) 1, 0.5, 0

            ciphertext_sub[i] = cc->EvalSub(ciphertext_a[i], ciphertext_b[i]);
            result_sub[i] = cc->EvalSign(ciphertext_sub[i], 3, -precision - 10, precision + 10, polyDegree);
            Homround(result_sub[i]);

            // comp(a==b) 1, 0

            auto ciphertext_sub_neg = cc->EvalNegate(result_sub[i]);
            auto result_sub_neg = cc->EvalAdd(ciphertext_sub_neg, 1.0);

            auto result_equal = cc->EvalMult(result_sub_neg, result_sub[i]);
            result_equal_comp[i] = MultByInteger(result_equal, 4.0);
            ciphertext_comp_greater[i] = cc->EvalSub(result_sub[i], result_equal_comp[i]);
        }
        // mux
        for (int i = 0; i < block; i++)
        {

            if (i == 0)
                res = result_equal_comp[i];
            else
            {
                res = cc->EvalMult(res, result_equal_comp[i]);
            }
        }
        Homround(res);
    }

    // ct1-ct2 >0 -> res=1; ct1-ct2 ==0 -> res=0; ct1-ct2 <0 -> res=0
    void comp_partial(Ciphertext<lbcrypto::DCRTPoly> &ct1, Ciphertext<lbcrypto::DCRTPoly> &ct2, double precision,
                      std::vector<double> &coefficients, Ciphertext<lbcrypto::DCRTPoly> &res,
                      lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        auto cc = ct1->GetCryptoContext();
        auto ciphertext_sub = cc->EvalSub(ct1, ct2);
        // 11+7 =18level

        auto result_sub = cc->EvalChebyshevSeries(ciphertext_sub, coefficients, -precision - 10, precision + 10);

        Homround(result_sub);
        auto ciphertext_greater_equal_2 = MultByInteger(result_sub, 2.0);

        // comp(a==b)

        auto ciphertext_sub_neg = cc->EvalNegate(result_sub);
        auto result_sub_neg = cc->EvalAdd(ciphertext_sub_neg, 1.0);

        auto result_equal = cc->EvalMult(result_sub_neg, result_sub);
        auto result_equal_4 = MultByInteger(result_equal, 2.0);
        res = cc->EvalSub(result_sub, result_equal_4);
    }

    void comp_partial_modular(std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ciphertext_a,
                              std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ciphertext_b, double precision,
                              std::vector<double> &coefficients, Ciphertext<lbcrypto::DCRTPoly> &comp_res,
                              lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        auto cc = ciphertext_a[0]->GetCryptoContext();
        int block = ciphertext_a.size();
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_sub(block), result_sub(block), result_equal_comp(block),
            ciphertext_comp_greater(block);
        for (int i = block - 1; i >= 0; i--)
        {
            // comp(a>b) 1, 0.5, 0
            ciphertext_sub[i] = cc->EvalSub(ciphertext_a[i], ciphertext_b[i]);
            result_sub[i] = cc->EvalChebyshevSeries(ciphertext_sub[i], coefficients, -precision - 10, precision + 10);
            Homround(result_sub[i]);

            // comp(a==b) 1, 0

            auto ciphertext_sub_neg = cc->EvalNegate(result_sub[i]);
            auto result_sub_neg = cc->EvalAdd(ciphertext_sub_neg, 1.0);

            auto result_equal = cc->EvalMult(result_sub_neg, result_sub[i]);
            result_equal_comp[i] = MultByInteger(result_equal, 2.0);
            ciphertext_comp_greater[i] = cc->EvalSub(result_sub[i], result_equal_comp[i]);
            result_equal_comp[i] = MultByInteger(result_equal_comp[i], 2.0);
        }

        for (int i = 0; i < block; i++)
        {
            if (i == 0)
                comp_res = ciphertext_comp_greater[i];
            else
            {
                auto ai_sub_bi = cc->EvalSub(comp_res, ciphertext_comp_greater[i]);
                auto equal_mul_sub = cc->EvalMult(result_equal_comp[i], ai_sub_bi);
                comp_res = cc->EvalAdd(equal_mul_sub, ciphertext_comp_greater[i]);
            }
        }
        // Homround(comp_res);
    }

    void EvalSignExample(std::uint32_t plain_bits)
    {
        std::cout << "--------------------------------- EVAL SIGN FUNCTION ---------------------------------"
                  << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128
        usint scalingModSize = 78;
        usint firstModSize = 89;
#else
        usint scalingModSize = 59;
        usint firstModSize = 60;
#endif
        parameters.SetScalingModSize(scalingModSize);
        parameters.SetFirstModSize(firstModSize);
        std::uint32_t polyDegree = 59;
        std::uint32_t multDepth = 27;

        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);
        double precision = (1 << (plain_bits - 1) - 1);
        double lowerBound = -precision;
        double upperBound = precision;
        double bound = 3;
        auto keyPair = cc->KeyGen();
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message_1bit(-1, 1);
        std::uniform_int_distribution<int> message_part_1(-2, -1);
        std::uniform_int_distribution<int> message_part_2(1, 2);
        std::uniform_int_distribution<int> message_part_16_2(1, upperBound);

        usint ringDim = cc->GetRingDimension();
        int length = 128;
        std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl
                  << std::endl;
        std::vector<std::complex<double>> input0_1(length), input_4bit(length), input_16bit(length);
        for (int i = 0; i < length; i++)
        {
            input0_1[i] = double(message_1bit(engine));
        }
        for (int i = 0; i < length / 2; i++)
        {
            input_4bit[2 * i] = double(message_part_1(engine));
            input_4bit[2 * i + 1] = double(message_part_2(engine));
        }
        for (int i = 0; i < length / 2; i++)
        {
            input_16bit[2 * i] = double(message_part_16_2(engine));
            input_16bit[2 * i + 1] = double(message_part_16_2(engine));
        }
        size_t encodedLength = input0_1.size();
        Plaintext plaintext = cc->MakeCKKSPackedPlaintext(input0_1);
        Plaintext plaintext_4 = cc->MakeCKKSPackedPlaintext(input_4bit);
        Plaintext plaintext_16 = cc->MakeCKKSPackedPlaintext(input_16bit);
        auto ciphertext = cc->Encrypt(keyPair.publicKey, plaintext);
        auto ciphertext_4 = cc->Encrypt(keyPair.publicKey, plaintext_4);
        auto ciphertext_16 = cc->Encrypt(keyPair.publicKey, plaintext_16);

        auto result = cc->EvalSign(ciphertext, bound, lowerBound, upperBound, polyDegree);
    }

    void EvalequalExample(std::uint32_t plain_bits)
    {
        std::cout << "\nplain_bits equal comp: " << plain_bits << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;

        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128
        usint scalingModSize = 78;
        usint firstModSize = 89;
#else
        usint scalingModSize = 59;
        usint firstModSize = 60;
#endif
        SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
        parameters.SetSecretKeyDist(secretKeyDist);
        parameters.SetScalingModSize(scalingModSize);
        parameters.SetFirstModSize(firstModSize);
        std::vector<std::uint32_t> levelBudget = {4, 4};
        std::uint32_t polyDegree = 1007;

        std::uint32_t multDepth = 25;
        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(FHE);
        cc->Enable(ADVANCEDSHE);

        double precision = (1 << (plain_bits - 1) - 1);
        double lowerBound = -precision;
        double upperBound = precision;
        double bound = 3;
        auto keyPair = cc->KeyGen();
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;
        cc->EvalMultKeyGen(keyPair.secretKey);

        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);

        std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;
        std::vector<std::complex<double>> input_a(length), input_b(length);

        for (int i = 0; i < length / 2; i++)
        {
            input_a[i] = double(message(engine));
            input_b[i] = double(message(engine));
        }
        size_t encodedLength = input_a.size();
        Plaintext plaintext_a = cc->MakeCKKSPackedPlaintext(input_a);
        Plaintext plaintext_b = cc->MakeCKKSPackedPlaintext(input_b);
        auto ciphertext_a = cc->Encrypt(keyPair.publicKey, plaintext_a);
        auto ciphertext_b = cc->Encrypt(keyPair.publicKey, plaintext_b);
        // comp(a>b),  a-b

        std::cout << "number of levels remaining before comp: " << multDepth - ciphertext_a->GetLevel() << std::endl;
        std::chrono::system_clock::time_point start, end;
        start = std::chrono::system_clock::now();
        auto ciphertext_sub = cc->EvalSub(ciphertext_a, ciphertext_b);
        auto result_sub = cc->EvalSign(ciphertext_sub, bound, lowerBound, upperBound, polyDegree);
        std::cout << "number of levels remaining before Homround: " << multDepth - result_sub->GetLevel() << std::endl;
        Homround(result_sub);
        std::cout << "number of levels remaining after Homround: " << multDepth - result_sub->GetLevel() << std::endl;
        // auto ciphertext_greater_equal_2 = cc->EvalMult(result_sub, 2.0);

        // comp(a==b)

        auto ciphertext_sub_neg = cc->EvalNegate(result_sub);
        auto result_sub_neg = cc->EvalAdd(ciphertext_sub_neg, 1.0);

        auto result_equal = cc->EvalMult(result_sub_neg, result_sub);
        auto result_equal_4 = MultByInteger(result_equal, 2.0);
        // auto result_equal_4 = cc->EvalMult(result_equal, 4.0);
        // Homround(result_equal_4);

        // result_equal_4 = cc->EvalBootstrap(result_equal_4);

        auto ciphertext_comp_greater = cc->EvalSub(result_sub, result_equal_4);
        end = std::chrono::system_clock::now();
        std::cout << "number of levels remaining after equal comp: " << multDepth - result_equal_4->GetLevel() << std::endl;
        // auto ciphertext_comp_greater = cc->EvalSub(ciphertext_greater_equal_2, result_equal_4);
        // ciphertext_comp_greater      = cc->EvalMult(ciphertext_comp_greater, 0.5);
        // Homround(ciphertext_comp_greater);
        std::cout << "number of levels remaining after greater comp: " << multDepth - ciphertext_comp_greater->GetLevel()
                  << std::endl;
        auto ciphertext_res_greater = cc->EvalMult(ciphertext_a, ciphertext_comp_greater);
        auto ciphertext_res_equal = cc->EvalMult(ciphertext_a, result_equal_4);

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
        std::vector<std::complex<double>> Result_comp_equal = plaintextDec_result_equal->GetCKKSPackedValue();
        std::vector<std::complex<double>> Result_comp_greater = plaintextDec_result_greater->GetCKKSPackedValue();
        std::vector<std::complex<double>> Result_equal = plaintextDec_comp_res->GetCKKSPackedValue();
        std::vector<std::complex<double>> Result_greater = plaintextDec_greater_res->GetCKKSPackedValue();
        double err = 0;
        double err_greater = 0;
        double err_comp_equal = 0;
        double err_comp_greater = 0;
        int wrong = 0;
        std::vector<double> expectedOutput_Result_equal, expectedOutput_Result_greater, expectedOutput_Comp_equal,
            expectedOutput_Comp_greater;
        for (int i = 0; i < (int)encodedLength; i++)
        {
            expectedOutput_Result_equal.push_back(
                static_cast<double>(input_a[i].real() == input_b[i].real() ? input_a[i].real() : 0));
            expectedOutput_Result_greater.push_back(input_a[i].real() > input_b[i].real() ? input_a[i].real() : 0);
            expectedOutput_Comp_equal.push_back(input_a[i].real() == input_b[i].real() ? 1 : 0);
            expectedOutput_Comp_greater.push_back(input_a[i].real() > input_b[i].real() ? 1 : 0);

            if (std::fabs(expectedOutput_Comp_equal[i] - Result_comp_equal[i].real()) > 0.5 ||
                std::fabs(expectedOutput_Comp_greater[i] - Result_comp_greater[i].real()) > 0.5)
            {
                wrong++;
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

    // 1,0.5,0
    void Eval_modular_greater_than(std::uint32_t plain_bits, std::uint32_t block)
    {
        std::cout << plain_bits << " plain_bits " << block << " modular greater test " << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128
        usint scalingModSize = 78;
        usint firstModSize = 89;
#else
        usint scalingModSize = 59;
        usint firstModSize = 60;
#endif
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetScalingModSize(scalingModSize);
        parameters.SetFirstModSize(firstModSize);
        std::uint32_t polyDegree = 0;
        std::uint32_t multDepth = 20 + block;

        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);
        double precision = (1 << (plain_bits - 1) - 1);
        double lowerBound = -precision;
        double upperBound = precision;
        double bound = 3;
        auto keyPair = cc->KeyGen();
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;

        std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;
        std::cout << "using polyDegree " << polyDegree << std::endl;
        std::vector<std::vector<std::complex<double>>> input_a_mod, input_b_mod;
        std::vector<double> input_a_plain(length, 0), input_b_plain(length, 0);
        // a0,a1...  b0,b1...

        for (int i = 0; i < block; i++)
        {
            std::vector<std::complex<double>> a(length), b(length);
            for (int j = 0; j < length; j++)
            {
                a[j] = double(message(engine));
                b[j] = double(message(engine));
            }
            input_a_mod.push_back(a);
            input_b_mod.push_back(b);
        }

        for (int i = 0; i < length; i++)
        {
            for (int j = 0; j < block; j++)
            {
                input_a_plain[i] += input_a_mod[j][i].real() * std::pow(2, j * plain_bits);
                input_b_plain[i] += input_b_mod[j][i].real() * std::pow(2, j * plain_bits);
            }
        }
        double totalerr = 0;
        double totaltime = 0;
        size_t encodedLength = input_a_plain.size();

        for (int num_test = 0; num_test < 10; num_test++)
        {
            std::vector<Plaintext> plaintext_a, plaintext_b;
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_a, ciphertext_b;

            std::chrono::system_clock::time_point start, end;
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_sub(block), result_sub(block), result_equal_comp(block),
                ciphertext_comp_greater(block);
            for (int i = 0; i < block; i++)
            {
                Plaintext plain_a = cc->MakeCKKSPackedPlaintext(input_a_mod[i]);
                Plaintext plain_b = cc->MakeCKKSPackedPlaintext(input_b_mod[i]);
                ciphertext_a.push_back(cc->Encrypt(keyPair.publicKey, plain_a));
                ciphertext_b.push_back(cc->Encrypt(keyPair.publicKey, plain_b));
            }
            start = std::chrono::system_clock::now();
            for (int i = block - 1; i >= 0; i--)
            {
                // comp(a>b) 1, 0.5, 0

                ciphertext_sub[i] = cc->EvalSub(ciphertext_a[i], ciphertext_b[i]);
                result_sub[i] = cc->EvalSign(ciphertext_sub[i], bound, lowerBound, upperBound, polyDegree);
                Homround(result_sub[i]);

                // comp(a==b) 1, 0

                auto ciphertext_sub_neg = cc->EvalNegate(result_sub[i]);
                auto result_sub_neg = cc->EvalAdd(ciphertext_sub_neg, 1.0);

                auto result_equal = cc->EvalMult(result_sub_neg, result_sub[i]);
                result_equal_comp[i] = MultByInteger(result_equal, 4.0);
            }

            Ciphertext<lbcrypto::DCRTPoly> comp_res;
            for (int i = 0; i < block; i++)
            {
                if (i == 0)
                    comp_res = result_sub[i];
                else
                {
                    auto ai_sub_bi = cc->EvalSub(comp_res, result_sub[i]);
                    auto equal_mul_sub = cc->EvalMult(result_equal_comp[i], ai_sub_bi);
                    comp_res = cc->EvalAdd(equal_mul_sub, result_sub[i]);
                }
            }

            end = std::chrono::system_clock::now();
            Plaintext plaintextDec, plaintextDec_comp_res, plaintextDec_comp_res_equal;

            cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec_comp_res);

            plaintextDec_comp_res->SetLength(encodedLength);
            std::vector<double> expectedOutput_compres;
            std::vector<std::complex<double>> finalResult = plaintextDec_comp_res->GetCKKSPackedValue();
            double err = 0;
            int wrong = 0;

            for (int i = 0; i < int(encodedLength); i++)
            {
                if (input_a_plain[i] != input_b_plain[i])
                    expectedOutput_compres.push_back(static_cast<double>(input_a_plain[i] > input_b_plain[i]));
                else
                    expectedOutput_compres.push_back(0.5);
                if (fabs(expectedOutput_compres[i] - finalResult[i].real()) > 0.5)
                {
                    wrong++;
                    std::cout << "input_a:" << input_a_plain[i] << ", input_b:" << input_b_plain[i]
                              << ", res :" << finalResult[i] << std::endl;
                }

                err += std::fabs(expectedOutput_compres[i] - finalResult[i].real());
                totalerr += err / encodedLength;
            }
            totaltime +=
                double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / double(encodedLength));
            if (wrong > 0)
            {
                std::cout << "avg error: 2^" << std::log2(err / encodedLength) << std::endl;
                std::cout << "wrong number:" << wrong << std::endl;
            }
        }
        std::cout << encodedLength << "slots amortize time: " << double((totaltime / 10.0)) << "ms" << std::endl
                  << std::endl;
        std::cout << "avg error: 2^" << std::log2(totalerr / 10) << std::endl;
    }

    double Eval_modular_Strictly_greater_than(std::uint32_t plain_bits, std::uint32_t block)
    {
        std::cout << plain_bits * block << " plain_bits greater comparison test " << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128
        usint scalingModSize = 78;
        usint firstModSize = 89;
#else
        usint scalingModSize = 50;
        usint firstModSize = 60;
#endif
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetScalingModSize(scalingModSize);
        parameters.SetFirstModSize(firstModSize);
        parameters.SetRingDim(65536 * 2);
        std::uint32_t polyDegree = 119;
        std::uint32_t multDepth = 20 + block;

        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);
        double precision = (1 << (plain_bits - 1)) - 1;
        double lowerBound = -precision;
        double upperBound = precision;
        double bound = 3;
        auto keyPair = cc->KeyGen();
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;

        std::vector<std::vector<std::complex<double>>> input_a_mod, input_b_mod;
        std::vector<double> input_a_plain(length, 0), input_b_plain(length, 0);
        // a0,a1...  b0,b1...

        for (int i = 0; i < block; i++)
        {
            std::vector<std::complex<double>> a(length), b(length);
            for (int j = 0; j < length; j++)
            {
                a[j] = double(message(engine));
                b[j] = double(message(engine));
            }
            input_a_mod.push_back(a);
            input_b_mod.push_back(b);
        }

        for (int i = 0; i < length; i++)
        {
            for (int j = 0; j < block; j++)
            {
                input_a_plain[i] += input_a_mod[j][i].real() * std::pow(2, j * plain_bits);
                input_b_plain[i] += input_b_mod[j][i].real() * std::pow(2, j * plain_bits);
            }
        }

        double totalerr = 0;
        double err_max = 0;
        double totaltime = 0;
        size_t encodedLength = input_a_plain.size();
        double num_test = 1;
        for (int num = 0; num < int(num_test); num++)
        {
            std::vector<Plaintext> plaintext_a, plaintext_b;
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_a, ciphertext_b;

            std::chrono::system_clock::time_point start, end;
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_sub(block), result_sub(block), result_equal_comp(block),
                ciphertext_comp_greater(block);
            for (int i = 0; i < block; i++)
            {
                Plaintext plain_a = cc->MakeCKKSPackedPlaintext(input_a_mod[i]);
                Plaintext plain_b = cc->MakeCKKSPackedPlaintext(input_b_mod[i]);
                ciphertext_a.push_back(cc->Encrypt(keyPair.publicKey, plain_a));
                ciphertext_b.push_back(cc->Encrypt(keyPair.publicKey, plain_b));
            }
            start = std::chrono::system_clock::now();
            Ciphertext<lbcrypto::DCRTPoly> comp_res;
            comp_greater_than_modular(ciphertext_a, ciphertext_b, precision, polyDegree, comp_res, keyPair.secretKey);
            end = std::chrono::system_clock::now();
            Plaintext plaintextDec, plaintextDec_comp_res, plaintextDec_comp_res_equal;

            cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec_comp_res);

            plaintextDec_comp_res->SetLength(encodedLength);
            std::vector<double> expectedOutput_compres;
            std::vector<std::complex<double>> finalResult = plaintextDec_comp_res->GetCKKSPackedValue();
            double err = 0;
            int wrong = 0;

            for (int i = 0; i < int(encodedLength); i++)
            {
                expectedOutput_compres.push_back(static_cast<double>(input_a_plain[i] > input_b_plain[i]));

                if (fabs(expectedOutput_compres[i] - finalResult[i].real()) > 0.5)
                {
                    wrong++;
                    std::cout << "input_a:" << input_a_plain[i] << ", input_b:" << input_b_plain[i]
                              << ", res :" << finalResult[i] << std::endl;
                }
                err += std::fabs(expectedOutput_compres[i] - finalResult[i].real());
                err_max = err_max > std::fabs(expectedOutput_compres[i] - finalResult[i].real()) ? err_max : std::fabs(expectedOutput_compres[i] - finalResult[i].real());
            }
            totalerr += err / encodedLength;
            totaltime +=
                double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / double(encodedLength));

            if (wrong > 0)
            {
                std::cout << "avg error: 2^" << std::log2(err / encodedLength) << std::endl;
                std::cout << "wrong number:" << wrong << std::endl;
            }
        }

        std::cout << encodedLength << "slots amortize time: " << double((totaltime / num_test)) << "ms" << std::endl;
        std::cout << "avg error: 2^" << std::log2(totalerr / num_test) << std::endl;
        std::cout << "max error: 2^" << std::log2(err_max / num_test) << std::endl
                  << std::endl
                  << std::endl;
        cc->ClearEvalSumKeys();
        cc->ClearEvalMultKeys();
        cc->ClearEvalAutomorphismKeys();
        return double((totaltime / num_test));
    }

    double Eval_modular_Strictly_equal(std::uint32_t plain_bits, std::uint32_t block)
    {
        std::cout << plain_bits * block << " plain_bits equality comparison test" << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;

        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128
        usint scalingModSize = 78;
        usint firstModSize = 89;
#else
        usint scalingModSize = 40;
        usint firstModSize = 60;
#endif
        SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
        parameters.SetScalingModSize(scalingModSize);
        parameters.SetFirstModSize(firstModSize);

        std::uint32_t polyDegree = 59;
        std::uint32_t multDepth = 20 + block;
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
        cc->Enable(ADVANCEDSHE);
        double precision = (1 << (plain_bits - 1) - 1);
        double lowerBound = -precision;
        double upperBound = precision;
        double bound = 3;
        auto keyPair = cc->KeyGen();
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;

        std::vector<std::vector<std::complex<double>>> input_a_mod, input_b_mod;
        std::vector<double> input_a_plain(length, 0), input_b_plain(length, 0);
        // a0,a1...  b0,b1...

        for (int i = 0; i < block; i++)
        {
            std::vector<std::complex<double>> a(length), b(length);
            for (int j = 0; j < length; j++)
            {
                a[j] = double(message(engine));
                b[j] = double(message(engine));
            }
            input_a_mod.push_back(a);
            input_b_mod.push_back(b);
        }

        for (int i = 0; i < length; i++)
        {
            for (int j = 0; j < block; j++)
            {
                input_a_plain[i] += input_a_mod[j][i].real() * std::pow(2, j * plain_bits);
                input_b_plain[i] += input_b_mod[j][i].real() * std::pow(2, j * plain_bits);
            }
        }
        double totalerr = 0;
        double err_max = 0;
        double totaltime = 0;
        size_t encodedLength = input_a_plain.size();
        double num_test = 1;
        for (int num = 0; num < int(num_test); num++)
        {
            std::vector<Plaintext> plaintext_a, plaintext_b;
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_a, ciphertext_b;

            std::chrono::system_clock::time_point start, end;
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_sub(block), result_sub(block), result_equal_comp(block),
                ciphertext_comp_greater(block);
            for (int i = 0; i < block; i++)
            {
                Plaintext plain_a = cc->MakeCKKSPackedPlaintext(input_a_mod[i]);
                Plaintext plain_b = cc->MakeCKKSPackedPlaintext(input_b_mod[i]);
                ciphertext_a.push_back(cc->Encrypt(keyPair.publicKey, plain_a));
                ciphertext_b.push_back(cc->Encrypt(keyPair.publicKey, plain_b));
            }
            start = std::chrono::system_clock::now();
            for (int i = block - 1; i >= 0; i--)
            {
                // comp(a>b) 1, 0.5, 0

                ciphertext_sub[i] = cc->EvalSub(ciphertext_a[i], ciphertext_b[i]);
                result_sub[i] = cc->EvalSign(ciphertext_sub[i], bound, lowerBound, upperBound, polyDegree);
                Homround(result_sub[i]);

                // comp(a==b) 1, 0

                auto ciphertext_sub_neg = cc->EvalNegate(result_sub[i]);
                auto result_sub_neg = cc->EvalAdd(ciphertext_sub_neg, 1.0);

                auto result_equal = cc->EvalMult(result_sub_neg, result_sub[i]);
                result_equal_comp[i] = MultByInteger(result_equal, 4.0);
                ciphertext_comp_greater[i] = cc->EvalSub(result_sub[i], result_equal_comp[i]);
            }
            // mux
            Ciphertext<lbcrypto::DCRTPoly> comp_res;
            for (int i = 0; i < block; i++)
            {

                if (i == 0)
                    comp_res = result_equal_comp[i];
                else
                {
                    comp_res = cc->EvalMult(comp_res, result_equal_comp[i]);
                }
            }
            Homround(comp_res);
            end = std::chrono::system_clock::now();
            Plaintext plaintextDec, plaintextDec_comp_res, plaintextDec_comp_res_equal;

            cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec_comp_res);

            plaintextDec_comp_res->SetLength(encodedLength);
            std::vector<double> expectedOutput_compres;
            std::vector<std::complex<double>> finalResult = plaintextDec_comp_res->GetCKKSPackedValue();
            double err = 0;

            int wrong = 0;

            for (int i = 0; i < int(encodedLength); i++)
            {
                expectedOutput_compres.push_back(static_cast<double>(input_a_plain[i] == input_b_plain[i]));

                if (fabs(expectedOutput_compres[i] - finalResult[i].real()) > 0.1)
                {
                    wrong++;
                    std::cout << "input_a:" << input_a_plain[i] << ", input_b:" << input_b_plain[i]
                              << ", res :" << finalResult[i] << std::endl;
                }
                err += std::fabs(expectedOutput_compres[i] - finalResult[i].real());
                err_max = err_max > std::fabs(expectedOutput_compres[i] - finalResult[i].real()) ? err_max : std::fabs(expectedOutput_compres[i] - finalResult[i].real());
            }
            totalerr += err / encodedLength;
            totaltime +=
                double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / double(encodedLength));

            if (wrong > 0)
            {
                std::cout << "avg error: 2^" << std::log2(err / encodedLength) << std::endl;
                std::cout << "wrong number:" << wrong << std::endl;
            }
        }
        std::cout << encodedLength << "slots amortize time: " << double((totaltime / num_test)) << "ms" << std::endl;
        std::cout << "avg error: 2^" << std::log2(totalerr / num_test) << std::endl;
        std::cout << "max error: 2^" << std::log2(err_max / num_test) << std::endl
                  << std::endl
                  << std::endl;
        cc->ClearEvalSumKeys();
        cc->ClearEvalMultKeys();
        cc->ClearEvalAutomorphismKeys();
        return double((totaltime / num_test));
    }

    void Eval_modular_greater_test(std::uint32_t plain_bits, std::uint32_t block)
    {
        std::cout << plain_bits << " plain_bits " << block << " modular greater test " << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;

        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128
        usint scalingModSize = 78;
        usint firstModSize = 89;
#else
        usint scalingModSize = 59;
        usint firstModSize = 60;
#endif

        parameters.SetScalingModSize(scalingModSize);
        parameters.SetFirstModSize(firstModSize);
        std::uint32_t polyDegree = 0;
        std::uint32_t multDepth = 23;
        // Choosing a higher degree yields better precision, but a longer runtime.
        if (plain_bits >= 8)
            polyDegree = 119;
        else
        {
            polyDegree = 27;
            multDepth = 20;
        }

        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);
        double precision = (1 << (plain_bits - 1) - 1);
        double lowerBound = -precision;
        double upperBound = precision;
        double bound = 3;
        auto keyPair = cc->KeyGen();
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;

        std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;
        std::cout << "using polyDegree " << polyDegree << std::endl;
        std::vector<std::vector<std::complex<double>>> input_a_mod, input_b_mod;
        std::vector<double> input_a_plain(length, 0), input_b_plain(length, 0);
        // a0,a1...  b0,b1...

        for (int i = 0; i < block; i++)
        {
            std::vector<std::complex<double>> a(length), b(length);
            for (int j = 0; j < length; j++)
            {
                a[j] = double(message(engine));
                b[j] = double(message(engine));
            }
            input_a_mod.push_back(a);
            input_b_mod.push_back(b);
        }

        for (int i = 0; i < length; i++)
        {
            for (int j = 0; j < block; j++)
            {
                input_a_plain[i] += input_a_mod[j][i].real() * std::pow(2, j * plain_bits);
                input_b_plain[i] += input_b_mod[j][i].real() * std::pow(2, j * plain_bits);
            }
        }
        double totalerr = 0;
        double totaltime = 0;
        size_t encodedLength = input_a_plain.size();

        for (int num_test = 0; num_test < 10; num_test++)
        {
            std::vector<Plaintext> plaintext_a, plaintext_b;
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_a, ciphertext_b;

            std::chrono::system_clock::time_point start, end;
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_sub(block), result_sub(block), result_equal_comp(block),
                ciphertext_comp_greater(block);
            for (int i = 0; i < block; i++)
            {
                Plaintext plain_a = cc->MakeCKKSPackedPlaintext(input_a_mod[i]);
                Plaintext plain_b = cc->MakeCKKSPackedPlaintext(input_b_mod[i]);
                ciphertext_a.push_back(cc->Encrypt(keyPair.publicKey, plain_a));
                ciphertext_b.push_back(cc->Encrypt(keyPair.publicKey, plain_b));
            }
            start = std::chrono::system_clock::now();
            for (int i = block - 1; i >= 0; i--)
            {
                // comp(a>b) 1, 0.5, 0

                ciphertext_sub[i] = cc->EvalSub(ciphertext_a[i], ciphertext_b[i]);
                result_sub[i] = cc->EvalSign(ciphertext_sub[i], bound, lowerBound, upperBound, polyDegree);
                Homround(result_sub[i]);

                // comp(a==b) 1, 0

                auto ciphertext_sub_neg = cc->EvalNegate(result_sub[i]);
                auto result_sub_neg = cc->EvalAdd(ciphertext_sub_neg, 1.0);

                auto result_equal = cc->EvalMult(result_sub_neg, result_sub[i]);
                result_equal_comp[i] = MultByInteger(result_equal, 4.0);
                // ciphertext_comp_greater[i] = cc->EvalSub(result_sub[i], result_equal_comp[i]);
            }

            Ciphertext<lbcrypto::DCRTPoly> comp_res;
            for (int i = 0; i < block; i++)
            {
                // std::cout << "mux:" << i << std::endl;
                if (i == 0)
                    comp_res = result_sub[i];
                else
                {
                    auto ai_sub_bi = cc->EvalSub(comp_res, result_sub[i]);
                    auto equal_mul_sub = cc->EvalMult(result_equal_comp[i], ai_sub_bi);
                    comp_res = cc->EvalAdd(equal_mul_sub, result_sub[i]);
                }
            }

            end = std::chrono::system_clock::now();
            Plaintext plaintextDec, plaintextDec_comp_res, plaintextDec_comp_res_equal;

            cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec_comp_res);

            plaintextDec_comp_res->SetLength(encodedLength);
            std::vector<double> expectedOutput_compres;
            std::vector<std::complex<double>> finalResult = plaintextDec_comp_res->GetCKKSPackedValue();
            double err = 0;
            int wrong = 0;

            for (int i = 0; i < int(encodedLength); i++)
            {
                if (input_a_plain[i] != input_b_plain[i])
                    expectedOutput_compres.push_back(static_cast<double>(input_a_plain[i] > input_b_plain[i]));
                else
                    expectedOutput_compres.push_back(0.5);
                if (fabs(expectedOutput_compres[i] - finalResult[i].real()) > 0.5)
                {
                    wrong++;
                    std::cout << "input_a:" << input_a_plain[i] << ", input_b:" << input_b_plain[i]
                              << ", res :" << finalResult[i] << std::endl;
                }

                err += std::fabs(expectedOutput_compres[i] - finalResult[i].real());
                totalerr += err / encodedLength;
            }
            totaltime +=
                double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / double(encodedLength));
            if (wrong > 0)
            {
                std::cout << "avg error: 2^" << std::log2(err / encodedLength) << std::endl;
                std::cout << "wrong number:" << wrong << std::endl;
            }
        }
        std::cout << encodedLength << "slots amortize time: " << double((totaltime / 10.0)) << "ms" << std::endl
                  << std::endl;
        std::cout << "avg error: 2^" << std::log2(totalerr / 10) << std::endl;
    }

    double Eval_modular_Strictly_greater_test(std::uint32_t plain_bits, std::uint32_t block)
    {
        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetSecurityLevel(HEStd_128_classic);
        // parameters.SetRingDim(1 << 16);
#if NATIVEINT == 128
        usint scalingModSize = 78;
        usint firstModSize = 89;
#else
        usint scalingModSize = 40;
        usint firstModSize = 60;
#endif
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetScalingModSize(scalingModSize);
        parameters.SetFirstModSize(firstModSize);
        parameters.SetRingDim(65536);
        std::uint32_t polyDegree = 59;
        std::uint32_t multDepth = 20 + block;

        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);
        double precision = (1 << (plain_bits - 1) - 1);
        double lowerBound = -precision;
        double upperBound = precision;
        double bound = 3;
        auto keyPair = cc->KeyGen();
        // We need to generate mult keys to run Chebyshev approximations.
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;

        std::vector<std::vector<std::complex<double>>> input_a_mod, input_b_mod;
        std::vector<double> input_a_plain(length, 0), input_b_plain(length, 0);
        // a0,a1...  b0,b1...

        for (int i = 0; i < block; i++)
        {
            std::vector<std::complex<double>> a(length), b(length);
            for (int j = 0; j < length; j++)
            {
                a[j] = double(message(engine));
                b[j] = double(message(engine));
            }
            input_a_mod.push_back(a);
            input_b_mod.push_back(b);
        }

        for (int i = 0; i < length; i++)
        {
            for (int j = 0; j < block; j++)
            {
                input_a_plain[i] += input_a_mod[j][i].real() * std::pow(2, j * plain_bits);
                input_b_plain[i] += input_b_mod[j][i].real() * std::pow(2, j * plain_bits);
            }
        }

        double totalerr = 0;
        double totaltime = 0;
        size_t encodedLength = input_a_plain.size();
        double num_test = 1;
        for (int num = 0; num < int(num_test); num++)
        {
            std::vector<Plaintext> plaintext_a, plaintext_b;
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_a, ciphertext_b;

            std::chrono::system_clock::time_point start, end;
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_sub(block), result_sub(block), result_equal_comp(block),
                ciphertext_comp_greater(block);
            for (int i = 0; i < block; i++)
            {
                Plaintext plain_a = cc->MakeCKKSPackedPlaintext(input_a_mod[i]);
                Plaintext plain_b = cc->MakeCKKSPackedPlaintext(input_b_mod[i]);
                ciphertext_a.push_back(cc->Encrypt(keyPair.publicKey, plain_a));
                ciphertext_b.push_back(cc->Encrypt(keyPair.publicKey, plain_b));
            }
            start = std::chrono::system_clock::now();
            for (int i = block - 1; i >= 0; i--)
            {
                // comp(a>b) 1, 0.5, 0

                ciphertext_sub[i] = cc->EvalSub(ciphertext_a[i], ciphertext_b[i]);
                result_sub[i] = cc->EvalSign(ciphertext_sub[i], bound, lowerBound, upperBound, polyDegree);
                Homround(result_sub[i]);

                // comp(a==b) 1, 0

                auto ciphertext_sub_neg = cc->EvalNegate(result_sub[i]);
                auto result_sub_neg = cc->EvalAdd(ciphertext_sub_neg, 1.0);

                auto result_equal = cc->EvalMult(result_sub_neg, result_sub[i]);
                result_equal_comp[i] = MultByInteger(result_equal, 2.0);
                ciphertext_comp_greater[i] = cc->EvalSub(result_sub[i], result_equal_comp[i]);
                result_equal_comp[i] = MultByInteger(result_equal_comp[i], 2.0);
            }

            Ciphertext<lbcrypto::DCRTPoly> comp_res;
            for (int i = 0; i < block; i++)
            {
                if (i == 0)
                    comp_res = ciphertext_comp_greater[i];
                else
                {
                    auto ai_sub_bi = cc->EvalSub(comp_res, ciphertext_comp_greater[i]);
                    auto equal_mul_sub = cc->EvalMult(result_equal_comp[i], ai_sub_bi);
                    comp_res = cc->EvalAdd(equal_mul_sub, ciphertext_comp_greater[i]);
                }
            }

            end = std::chrono::system_clock::now();
            Plaintext plaintextDec, plaintextDec_comp_res, plaintextDec_comp_res_equal;

            cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec_comp_res);

            plaintextDec_comp_res->SetLength(encodedLength);
            std::vector<double> expectedOutput_compres;
            std::vector<std::complex<double>> finalResult = plaintextDec_comp_res->GetCKKSPackedValue();
            double err = 0;
            int wrong = 0;

            for (int i = 0; i < int(encodedLength); i++)
            {
                expectedOutput_compres.push_back(static_cast<double>(input_a_plain[i] > input_b_plain[i]));

                if (fabs(expectedOutput_compres[i] - finalResult[i].real()) > 0.5)
                {
                    wrong++;
                }

                err += std::fabs(expectedOutput_compres[i] - finalResult[i].real());
                totalerr += err / encodedLength;
            }
            totaltime +=
                double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / double(encodedLength));
            // if (wrong > 0)
            // {
            //     std::cout << "avg error: 2^" << std::log2(err / encodedLength) << std::endl;
            //     std::cout << "wrong number:" << wrong << std::endl;
            // }
        }
        cc->ClearEvalSumKeys();
        cc->ClearEvalMultKeys();
        cc->ClearEvalAutomorphismKeys();
        CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
        return totaltime / num_test;
    }

    double Eval_modular_Strictly_equal_test(std::uint32_t plain_bits, std::uint32_t block)
    {
        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128
        usint scalingModSize = 78;
        usint firstModSize = 89;
#else
        usint scalingModSize = 40;
        usint firstModSize = 60;
#endif
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetScalingModSize(scalingModSize);
        parameters.SetFirstModSize(firstModSize);
        parameters.SetRingDim(65536 * 2);
        std::uint32_t polyDegree = 59;
        std::uint32_t multDepth = 20 + block;
        // Choosing a higher degree yields better precision, but a longer runtime.
        // if (plain_bits >= 8)
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
        cc->Enable(ADVANCEDSHE);
        double precision = (1 << (plain_bits - 1) - 1);
        double lowerBound = -precision;
        double upperBound = precision;
        double bound = 3;
        auto keyPair = cc->KeyGen();
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;

        std::vector<std::vector<std::complex<double>>> input_a_mod, input_b_mod;
        std::vector<double> input_a_plain(length, 0), input_b_plain(length, 0);
        // a0,a1...  b0,b1...

        for (int i = 0; i < block; i++)
        {
            std::vector<std::complex<double>> a(length), b(length);
            for (int j = 0; j < length; j++)
            {
                a[j] = double(message(engine));
                b[j] = double(message(engine));
            }
            input_a_mod.push_back(a);
            input_b_mod.push_back(b);
        }

        for (int i = 0; i < length; i++)
        {
            for (int j = 0; j < block; j++)
            {
                input_a_plain[i] += input_a_mod[j][i].real() * std::pow(2, j * plain_bits);
                input_b_plain[i] += input_b_mod[j][i].real() * std::pow(2, j * plain_bits);
            }
        }
        double totalerr = 0;
        double totaltime = 0;
        size_t encodedLength = input_a_plain.size();
        double num_test = 1;
        for (int num = 0; num < int(num_test); num++)
        {
            std::vector<Plaintext> plaintext_a, plaintext_b;
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_a, ciphertext_b;

            std::chrono::system_clock::time_point start, end;
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_sub(block), result_sub(block), result_equal_comp(block),
                ciphertext_comp_greater(block);
            for (int i = 0; i < block; i++)
            {
                Plaintext plain_a = cc->MakeCKKSPackedPlaintext(input_a_mod[i]);
                Plaintext plain_b = cc->MakeCKKSPackedPlaintext(input_b_mod[i]);
                ciphertext_a.push_back(cc->Encrypt(keyPair.publicKey, plain_a));
                ciphertext_b.push_back(cc->Encrypt(keyPair.publicKey, plain_b));
            }
            start = std::chrono::system_clock::now();
            for (int i = block - 1; i >= 0; i--)
            {
                // comp(a>b) 1, 0.5, 0

                ciphertext_sub[i] = cc->EvalSub(ciphertext_a[i], ciphertext_b[i]);
                result_sub[i] = cc->EvalSign(ciphertext_sub[i], bound, lowerBound, upperBound, polyDegree);
                Homround(result_sub[i]);

                // comp(a==b) 1, 0

                auto ciphertext_sub_neg = cc->EvalNegate(result_sub[i]);
                auto result_sub_neg = cc->EvalAdd(ciphertext_sub_neg, 1.0);

                auto result_equal = cc->EvalMult(result_sub_neg, result_sub[i]);
                result_equal_comp[i] = MultByInteger(result_equal, 4.0);
                ciphertext_comp_greater[i] = cc->EvalSub(result_sub[i], result_equal_comp[i]);
            }

            Ciphertext<lbcrypto::DCRTPoly> comp_res;
            for (int i = 0; i < block; i++)
            {
                // std::cout << "mux:" << i << std::endl;
                if (i == 0)
                    comp_res = result_equal_comp[i];
                else
                {
                    comp_res = cc->EvalMult(comp_res, result_equal_comp[i]);
                }
            }

            end = std::chrono::system_clock::now();
            Plaintext plaintextDec, plaintextDec_comp_res, plaintextDec_comp_res_equal;

            cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec_comp_res);

            plaintextDec_comp_res->SetLength(encodedLength);
            std::vector<double> expectedOutput_compres;
            std::vector<std::complex<double>> finalResult = plaintextDec_comp_res->GetCKKSPackedValue();
            double err = 0;
            int wrong = 0;

            for (int i = 0; i < int(encodedLength); i++)
            {
                expectedOutput_compres.push_back(static_cast<double>(input_a_plain[i] == input_b_plain[i]));

                if (fabs(expectedOutput_compres[i] - finalResult[i].real()) > 0.5)
                {
                    wrong++;
                }

                err += std::fabs(expectedOutput_compres[i] - finalResult[i].real());
                totalerr += err / encodedLength;
            }
            totaltime +=
                double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / double(encodedLength));

            if (wrong > 0)
            {
                // std::cout << "avg error: 2^" << std::log2(err / encodedLength) << std::endl;
                // std::cout << "wrong number:" << wrong << std::endl;
            }
        }
        cc->ClearEvalSumKeys();
        cc->ClearEvalMultKeys();
        cc->ClearEvalAutomorphismKeys();
        CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
        return totaltime / num_test;
    }
}