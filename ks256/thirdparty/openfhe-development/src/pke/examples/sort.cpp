#include "openfhe.h"
#include <chrono>
#include "math/chebyshev.h"
#include <stdexcept>
using namespace lbcrypto;
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

Ciphertext<DCRTPoly> MultByInteger(const ConstCiphertext<DCRTPoly>& ciphertext, const int64_t constant) {
    auto result = ciphertext->Clone();
    for (auto& c : result->GetElements())
        c *= static_cast<typename DCRTPoly::Integer>(constant);
    return result;
}
std::vector<double> coeff1 = {1.5, -0.5};  // 1.5 x  - 0.5 x ^ 3
std::vector<double> coeff3 = {2.1875, -2.1875, 1.3125, -0.3125};
std::vector<double> coeff5 = {2.4609375 / 2, -3.28125 / 2, 2.953125 / 2, -1.40625 / 2, 0.2734375 / 2};
// -0.2095 x ^ 15 + 1.692 x ^ 13 + -5.999 x ^ 11 + 12.22 x ^ 9 + -15.71 x ^ 7 + 13.2 x ^ 5 + -7.332 x ^ 3 + 3.142 x
void EvalPower(std::vector<double> coefficients, std::vector<Ciphertext<lbcrypto::DCRTPoly>>& power_basis,
               Ciphertext<lbcrypto::DCRTPoly>& result) {
    auto cc = power_basis[0]->GetCryptoContext();
    if (coefficients.size() == 1) {
        if (std::abs(std::round(coefficients[0] * pow(2, 50))) > 1.) {
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

double CalculateApproximationError(const std::vector<std::complex<double>>& result,
                                   const std::vector<double>& expectedResult) {
    // using the infinity norm
    double maxError = 0;
    for (size_t i = 0; i < result.size(); ++i) {
        double error = std::abs(result[i].real() - expectedResult[i]);
        if (maxError < error)
            maxError = error;
    }

    return std::abs(std::log2(maxError));
}
//check first boot precision
int itboot(Ciphertext<lbcrypto::DCRTPoly>& res, lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey) {
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

void Homround(Ciphertext<lbcrypto::DCRTPoly>& cipher) {
    auto cc = cipher->GetCryptoContext();
    cipher  = MultByInteger(cipher, 2);
    cipher  = cc->EvalAdd(cipher, -1.0);
    poly_evaluate_power(cipher, cipher, coeff1);
    // poly_evaluate_power(cipher, scale, cipher, coeff3, context, encoder, evaluator, relin_keys, decryptor);
    poly_evaluate_power(cipher, cipher, coeff5);
    cipher = cc->EvalAdd(cipher, 0.5);
}

void HomComp(Ciphertext<lbcrypto::DCRTPoly>& a, Ciphertext<lbcrypto::DCRTPoly>& b, double precision,
             uint32_t polyDegree) {
    auto cc             = a->GetCryptoContext();
    auto ciphertext_sub = cc->EvalSub(a, b);

    auto result_sub = cc->EvalSign(ciphertext_sub, 0, -precision, precision, polyDegree);
    Homround(result_sub);
}

std::vector<std::vector<std::vector<double>>> bitonic_test_plain(std::vector<double> vec_unsort,
                                                                 std::vector<std::vector<double>>& plain_sort,
                                                                 uint32_t pack_slots) {
    // vector<int> vec_unsort = {5, 2, 1, 7, 3, 8, 6, 4};
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<> dis(1, 10000);
    std::vector<std::vector<std::vector<double>>> martix;
    int len = vec_unsort.size();
    for (int stage = 0; pack_slots >> stage + 1 > 0; stage++) {
        for (int part = 0; stage - part >= 0; part++) {
            // printf("stage: %d part: %d\n", stage, part);
            std::vector<std::vector<double>> swap_martix;  // dimension = 2 or 3, depends on stage and part
            int Block = 1 << (stage + 1);
            int block = 1 << (stage + 1 - part);
            // printf("block: %d\n", block);
            std::vector<double> swap_vector(len, 0);
            std::vector<double> eliminate(len, 0), eliminate_rot(len, 0);
            //1level
            for (int i = 0; i < len; i += block) {
                for (int j = 0; j < block / 2; j++) {
                    eliminate[i + j] = 1;
                }
            }

            for (int i = 0; i < len / block; i++) {
                for (int j = 0; j < block / 2; j++) {
                    int swap_bit;
                    if ((part == 0 && i % 2 == 0) || (part > 0 && (i * block / Block) % 2 == 0))
                        swap_bit = vec_unsort[i * block + j] < vec_unsort[i * block + block / 2 + j] ? 1 : 0;
                    else
                        swap_bit = vec_unsort[i * block + j] > vec_unsort[i * block + block / 2 + j] ? 1 : 0;
                    swap_vector[i * block + j]             = swap_bit;
                    swap_vector[i * block + block / 2 + j] = swap_bit;
                }
            }

            swap_martix.push_back(swap_vector);
            std::vector<double> swap_vector_upper(len, 0);
            std::vector<double> swap_vector_below(len, 0);
            std::vector<double> swap_vector_tmp(len, 0);

            // if (block < len) {
            for (int i = 0; i < len / block; i++)
                for (int j = 0; j < block / 2; j++)
                    swap_vector_tmp[i * block + j] = 1;
            for (int i = 0; i < len; i++)
                swap_vector_below[i] = swap_vector_tmp[i] * (swap_vector_tmp[i] - swap_martix[0][i]);
            rotate_copy(swap_vector_below.begin(), swap_vector_below.begin() + len - block / 2, swap_vector_below.end(),
                        swap_vector_upper.begin());
            swap_martix.push_back(swap_vector_below);
            swap_martix.push_back(swap_vector_upper);
            // }
            // else {
            //     for (int i = 0; i < len; i++)
            //         swap_vector_tmp[i] = 1;
            //     for (int i = 0; i < len; i++)
            //         swap_vector_upper[i] = swap_vector_tmp[i] - swap_martix[0][i];
            //     swap_martix.push_back(swap_vector_upper);
            //     swap_martix.push_back(swap_vector_upper);
            // }
            std::vector<double> vec_sort_temp(len, 0);

            // if (block < len) {
            std::vector<double> vec_rot_upper(len, 0), vec_rot_below(len, 0);
            rotate_copy(vec_unsort.begin(), vec_unsort.begin() + block / 2, vec_unsort.end(), vec_rot_below.begin());
            rotate_copy(vec_unsort.begin(), vec_unsort.begin() + len - block / 2, vec_unsort.end(),
                        vec_rot_upper.begin());

            for (int i = 0; i < len; i++)
                vec_unsort[i] = swap_martix[0][i] * vec_unsort[i] + swap_martix[1][i] * vec_rot_below[i] +
                                swap_martix[2][i] * vec_rot_upper[i];
            // }
            // else {
            //     std::vector<double> vec_rot_below(len, 0);
            //     rotate_copy(vec_unsort.begin(), vec_unsort.begin() + block / 2, vec_unsort.end(),
            //                 vec_rot_below.begin());
            //     for (int i = 0; i < len; i++)
            //         vec_unsort[i] = swap_martix[0][i] * vec_unsort[i] + swap_martix[1][i] * vec_rot_below[i];
            // }
            plain_sort.push_back(vec_unsort);
            martix.push_back(swap_martix);
        }
    }
    return martix;
}

// ct1>ct2?1:0
void comp_greater_than(Ciphertext<lbcrypto::DCRTPoly>& ct1, Ciphertext<lbcrypto::DCRTPoly>& ct2, double precision,
                       std::vector<double> coefficients, Ciphertext<lbcrypto::DCRTPoly>& res,
                       lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey) {
    auto cc = ct1->GetCryptoContext();
    //ct1-ct2 >0 -> 1; ct1-ct2 ==0 -> 0.5; ct1-ct2 <0 -> 0
    auto ciphertext_sub  = cc->EvalSub(ct1, ct2);
    auto ciphertext_sign = cc->EvalChebyshevSeries(ciphertext_sub, coefficients, -precision, precision);
    Homround(ciphertext_sign);
    //ct1-ct2 >0 -> 2; ct1-ct2 ==0 -> 1; ct1-ct2 <0 -> 0
    auto ciphertext_greater_equal_2 = MultByInteger(ciphertext_sign, 2.0);
    //comp(a==b)
    //ct2-ct1
    auto ciphertext_sub_neg = cc->EvalNegate(ciphertext_sign);
    auto result_sub_neg     = cc->EvalAdd(ciphertext_sub_neg, 1.0);
    //(ct2-ct1)*(ct1-ct2)
    auto result_equal   = cc->EvalMult(result_sub_neg, ciphertext_sign);
    auto result_equal_4 = MultByInteger(result_equal, 4.0);
    //get equal ct1-ct2 ==0 -> 1, else 0
    Homround(result_equal_4);

    auto ciphertext_comp_greater = cc->EvalSub(ciphertext_greater_equal_2, result_equal_4);
    res                          = cc->EvalMult(ciphertext_comp_greater, 0.5);
    int bootprecision            = itboot(res, privateKey);
    std::cout << "compres bootprecision : " << bootprecision << std::endl;
    res = cc->EvalBootstrap(res, 2, bootprecision);
    // Homround(res);
}

// ct1==ct2?1:0
void comp_equal(Ciphertext<lbcrypto::DCRTPoly>& ct1, Ciphertext<lbcrypto::DCRTPoly>& ct2, double precision,
                std::vector<double> coefficients, Ciphertext<lbcrypto::DCRTPoly>& res) {
    auto cc = ct1->GetCryptoContext();
    //ct1-ct2 >0 -> 1; ct1-ct2 ==0 -> 0.5; ct1-ct2 <0 -> 0
    auto ciphertext_sub  = cc->EvalSub(ct1, ct2);
    auto ciphertext_sign = cc->EvalChebyshevSeries(ciphertext_sub, coefficients, -precision, precision);
    Homround(ciphertext_sign);
    //ct1-ct2 >0 -> 2; ct1-ct2 ==0 -> 1; ct1-ct2 <0 -> 0
    auto ciphertext_greater_equal_2 = MultByInteger(ciphertext_sign, 2);
    //comp(a==b)
    //ct2-ct1
    auto ciphertext_sub_neg = cc->EvalNegate(ciphertext_sign);
    auto result_sub_neg     = cc->EvalAdd(ciphertext_sub_neg, 1.0);
    //(ct2-ct1)*(ct1-ct2)
    auto result_equal = cc->EvalMult(result_sub_neg, ciphertext_sign);
    res               = MultByInteger(result_equal, 4);
    //get equal ct1-ct2 ==0 -> 1, else 0
    Homround(res);
}

//ct1-ct2 >0 -> res=1; ct1-ct2 ==0 -> res=0; ct1-ct2 <0 -> res=0
void comp_partial(Ciphertext<lbcrypto::DCRTPoly>& ct1, Ciphertext<lbcrypto::DCRTPoly>& ct2, double precision,
                  std::vector<double>& coefficients, Ciphertext<lbcrypto::DCRTPoly>& res,
                  lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey) {
    auto cc             = ct1->GetCryptoContext();
    auto ciphertext_sub = cc->EvalSub(ct1, ct2);
    //11+7 =18level

    auto result_sub = cc->EvalChebyshevSeries(ciphertext_sub, coefficients, -precision - 10, precision + 10);

    Homround(result_sub);
    auto ciphertext_greater_equal_2 = MultByInteger(result_sub, 2.0);

    //comp(a==b)

    auto ciphertext_sub_neg = cc->EvalNegate(result_sub);
    auto result_sub_neg     = cc->EvalAdd(ciphertext_sub_neg, 1.0);

    auto result_equal   = cc->EvalMult(result_sub_neg, result_sub);
    auto result_equal_4 = MultByInteger(result_equal, 2.0);
    res                 = cc->EvalSub(result_sub, result_equal_4);
}

void comp_partial_modular(std::vector<Ciphertext<lbcrypto::DCRTPoly>>& ciphertext_a,
                          std::vector<Ciphertext<lbcrypto::DCRTPoly>>& ciphertext_b, double precision,
                          std::vector<double>& coefficients, Ciphertext<lbcrypto::DCRTPoly>& comp_res,
                          lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey) {
    auto cc   = ciphertext_a[0]->GetCryptoContext();
    int block = ciphertext_a.size();
    std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_sub(block), result_sub(block), result_equal_comp(block),
        ciphertext_comp_greater(block);
    for (int i = block - 1; i >= 0; i--) {
        //comp(a>b) 1, 0.5, 0
        ciphertext_sub[i] = cc->EvalSub(ciphertext_a[i], ciphertext_b[i]);
        result_sub[i]     = cc->EvalChebyshevSeries(ciphertext_sub[i], coefficients, -precision - 10, precision + 10);
        Homround(result_sub[i]);

        //comp(a==b) 1, 0

        auto ciphertext_sub_neg = cc->EvalNegate(result_sub[i]);
        auto result_sub_neg     = cc->EvalAdd(ciphertext_sub_neg, 1.0);

        auto result_equal          = cc->EvalMult(result_sub_neg, result_sub[i]);
        result_equal_comp[i]       = MultByInteger(result_equal, 2.0);
        ciphertext_comp_greater[i] = cc->EvalSub(result_sub[i], result_equal_comp[i]);
        result_equal_comp[i]       = MultByInteger(result_equal_comp[i], 2.0);
    }

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
}

void output(Ciphertext<lbcrypto::DCRTPoly>& res, int slot, lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey) {
    Plaintext plaintextDec;
    auto cc = res->GetCryptoContext();
    cc->Decrypt(privateKey, res, &plaintextDec);
    std::vector<std::complex<double>> dec = plaintextDec->GetCKKSPackedValue();
    for (int i = 0; i < slot; i++) {
        std::cout << "dec: " << dec[i].real() << std::endl;
    }
}
//
void bitonic_comp(int stage, int part, int slots, Ciphertext<lbcrypto::DCRTPoly>& ct,
                  Ciphertext<lbcrypto::DCRTPoly>& res, double precision, std::vector<double>& coefficients,
                  lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey) {
    auto cc   = ct->GetCryptoContext();
    int Block = 1 << (stage + 1);
    int block = 1 << (stage + 1 - part);
    std::vector<double> mask(slots, 1), eliminate(slots, 0);
    // if (block < slots)
    for (int i = 0; i < slots; i += 2 * Block) {
        for (int j = 0; j < Block; j++) {
            mask[i + j] = -1;
        }
    }
    // else

    for (int i = 0; i < slots; i += block) {
        for (int j = 0; j < block / 2; j++) {
            eliminate[i + j] = 1;
        }
    }
    // for (int i = 0; i < 8; i++) {
    //     std::cout << "i:" << i << " mask: " << mask[i] << ", eliminate: " << eliminate[i] << std::endl;
    // }
    Plaintext plain_mask      = cc->MakeCKKSPackedPlaintext(mask);
    Plaintext plain_eliminate = cc->MakeCKKSPackedPlaintext(eliminate);
    auto ct_mask              = cc->EvalMult(ct, plain_mask);
    // std::cout << "\n EvalRotate: " << block / 2 << std::endl;
    auto ct_rot = cc->EvalRotate(ct_mask, block / 2);

    // std::cout << "\n comp_partial " << std::endl;
    ct_mask = cc->EvalMult(ct_mask, plain_eliminate);

    ct_rot = cc->EvalMult(ct_rot, plain_eliminate);

    // std::cout << " before sign  " << std::endl;
    // output(ct_rot, 8, privateKey);
    // std::cout << std::endl;
    Ciphertext<lbcrypto::DCRTPoly> ploy_res;
    // comp_greater_than(ct_mask, ct_rot, precision, coefficients, ploy_res, privateKey);
    // comp_partial(ct_mask, ct_rot, precision, coefficients, ploy_res, privateKey);
    comp_partial(ct_mask, ct_rot, precision, coefficients, res, privateKey);
    // res = cc->EvalMult(ploy_res, plain_eliminate);

    // std::cout << "compres  " << std::endl;
    // output(res, 8, privateKey);
    // std::cout << std::endl;

    // auto ct_rot = cc->EvalRotate(res, -Block / 2);
    // res         = cc->EvalAdd(res, ct_rot);
}

void bitonic_comp_modular(int stage, int part, int slots, std::vector<Ciphertext<lbcrypto::DCRTPoly>>& ct,
                          Ciphertext<lbcrypto::DCRTPoly>& res, double precision, std::vector<double>& coefficients,
                          lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey) {
    auto cc    = ct[0]->GetCryptoContext();
    int Block  = 1 << (stage + 1);
    int block  = 1 << (stage + 1 - part);
    int num_ct = ct.size();
    std::vector<double> mask(slots, 1), eliminate(slots, 0);
    for (int i = 0; i < slots; i += 2 * Block) {
        for (int j = 0; j < Block; j++) {
            mask[i + j] = -1;
        }
    }
    // else

    for (int i = 0; i < slots; i += block) {
        for (int j = 0; j < block / 2; j++) {
            eliminate[i + j] = 1;
        }
    }
    // for (int i = 0; i < 8; i++) {
    //     std::cout << "i:" << i << " mask: " << mask[i] << ", eliminate: " << eliminate[i] << std::endl;
    // }
    std::vector<Ciphertext<lbcrypto::DCRTPoly>> ct_mask(num_ct), ct_rot(num_ct);
    Plaintext plain_mask      = cc->MakeCKKSPackedPlaintext(mask);
    Plaintext plain_eliminate = cc->MakeCKKSPackedPlaintext(eliminate);
    for (int i = 0; i < num_ct; i++) {
        ct_mask[i] = cc->EvalMult(ct[i], plain_mask);
        // std::cout << "\n EvalRotate: " << block / 2 << std::endl;
        ct_rot[i] = cc->EvalRotate(ct_mask[i], block / 2);

        // std::cout << "\n comp_partial " << std::endl;
        ct_mask[i] = cc->EvalMult(ct_mask[i], plain_eliminate);

        ct_rot[i] = cc->EvalMult(ct_rot[i], plain_eliminate);
    }
    Ciphertext<lbcrypto::DCRTPoly> ploy_res;
    comp_partial_modular(ct_mask, ct_rot, precision, coefficients, ploy_res, privateKey);
    res = cc->EvalMult(ploy_res, plain_eliminate);
}

void inverse(Ciphertext<lbcrypto::DCRTPoly>& res) {
    auto cc  = res->GetCryptoContext();
    auto neg = cc->EvalNegate(res);

    res = cc->EvalAdd(1, res);
}

void bitonic_swap(int stage, int part, int slots, Ciphertext<lbcrypto::DCRTPoly>& compres,
                  Ciphertext<lbcrypto::DCRTPoly>& ct, Ciphertext<lbcrypto::DCRTPoly>& res,
                  lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey) {
    std::chrono::system_clock::time_point start, end;
    double time = 0.;
    start       = std::chrono::system_clock::now();
    auto cc     = ct->GetCryptoContext();
    // int Block = 1 << (stage + 1);
    int block = 1 << (stage + 1 - part);
    std::vector<Ciphertext<lbcrypto::DCRTPoly>> res_martrix;
    std::vector<double> eliminate(slots, 0), eliminate_rot(slots, 0);
    //1level
    for (int i = 0; i < slots; i += block) {
        for (int j = 0; j < block / 2; j++) {
            eliminate[i + j] = 1;
        }
    }
    std::rotate_copy(eliminate.begin(), eliminate.begin() + block / 2, eliminate.end(), eliminate_rot.begin());

    Plaintext plaintextDec;
    Plaintext plain_eliminate = cc->MakeCKKSPackedPlaintext(eliminate);
    end                       = std::chrono::system_clock::now();
    time                      = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    // std::cout << "preprocess time: " << time << "ms" << std::endl;
    // Plaintext plain_eliminate_rot = cc->MakeCKKSPackedPlaintext(eliminate_rot);
    // if (block < slots) {

    auto ct_upper = cc->EvalRotate(ct, -block / 2);
    auto ct_lower = cc->EvalRotate(ct, block / 2);

    // ct_upper = cc->EvalMult(ct_upper, plain_eliminate_rot);
    // ct_lower = cc->EvalMult(ct_lower, plain_eliminate);

    std::vector<std::complex<double>> enc;
    start            = std::chrono::system_clock::now();
    auto compres_neg = cc->EvalNegate(compres);
    auto compres_rot = cc->EvalRotate(compres, -block / 2);

    auto swap_vec_lower = cc->EvalAdd(compres_neg, plain_eliminate);

    auto swap_vec_upper = cc->EvalRotate(swap_vec_lower, -block / 2);

    // cc->Decrypt(privateKey, swap_vec_upper, &plaintextDec);
    // plaintextDec->SetLength(slots);
    // enc = plaintextDec->GetCKKSPackedValue();
    // for (int i = 0; i < 8; i++) {
    //     std::cout << " swap_vec_upper: " << enc[i].real() << std::endl;
    // }
    auto swap_vec_diagonal = cc->EvalAdd(compres, compres_rot);
    // cc->Decrypt(privateKey, swap_vec_diagonal, &plaintextDec);
    // plaintextDec->SetLength(slots);
    // enc = plaintextDec->GetCKKSPackedValue();
    // for (int i = 0; i < 8; i++) {
    //     std::cout << " swap_vec_diagonal: " << enc[i].real() << std::endl;
    // }
    // inverse(swap_vec_upper);
    // inverse(swap_vec_lower);
    end  = std::chrono::system_clock::now();
    time = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    // std::cout << "rotate time: " << time << "ms" << std::endl;

    start = std::chrono::system_clock::now();
    res_martrix.push_back(cc->EvalMult(ct_upper, swap_vec_upper));
    res_martrix.push_back(cc->EvalMult(ct, swap_vec_diagonal));
    res_martrix.push_back(cc->EvalMult(ct_lower, swap_vec_lower));
    // cc->Decrypt(privateKey, res_martrix[0], &plaintextDec);
    // plaintextDec->SetLength(slots);
    // enc = plaintextDec->GetCKKSPackedValue();
    // for (int i = 0; i < 8; i++) {
    //     std::cout << " res_martrix[0]: " << enc[i].real() << std::endl;
    // }
    // cc->Decrypt(privateKey, res_martrix[1], &plaintextDec);
    // plaintextDec->SetLength(slots);
    // enc = plaintextDec->GetCKKSPackedValue();
    // for (int i = 0; i < 8; i++) {
    //     std::cout << " res_martrix[1]: " << enc[i].real() << std::endl;
    // }
    // cc->Decrypt(privateKey, res_martrix[2], &plaintextDec);
    // plaintextDec->SetLength(slots);
    // enc = plaintextDec->GetCKKSPackedValue();
    // for (int i = 0; i < 8; i++) {
    //     std::cout << " res_martrix[2]: " << enc[i].real() << std::endl;
    // }
    res  = cc->EvalAdd(res_martrix[0], res_martrix[1]);
    res  = cc->EvalAdd(res, res_martrix[2]);
    end  = std::chrono::system_clock::now();
    time = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    // std::cout << "mul and add time: " << time << "ms" << std::endl;
    // }
    // else {
    //     auto ct_partial_upper  = cc->EvalMult(res, plain_eliminate);
    //     auto ct_upper          = cc->EvalRotate(res, -block / 2);
    //     auto swap_vec_upper    = cc->EvalRotate(compres, -block / 2);
    //     auto swap_vec_diagonal = cc->EvalAdd(compres, swap_vec_upper);
    //     res_martrix.push_back(cc->EvalMult(res, swap_vec_diagonal));
    //     inverse(swap_vec_diagonal);
    //     res_martrix.push_back(cc->EvalMult(ct_upper, swap_vec_diagonal));
    //     res = cc->EvalAddMany(res_martrix);
    // }
}

double error_estimate(std::vector<double> plain, Ciphertext<lbcrypto::DCRTPoly>& res,
                      lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey, size_t num_slots) {
    auto cc = res->GetCryptoContext();
    std::vector<double> error_vec;
    Plaintext plaintextDec;
    cc->Decrypt(privateKey, res, &plaintextDec);
    // std::cout << "Estimated precision: " << plaintextDec->GetLogPrecision() << std::endl;
    std::vector<std::complex<double>> Result = plaintextDec->GetCKKSPackedValue();
    // for (int i = 0; i < int(num_slots); i++) {
    //     std::cout << "Result dec: " << Result[i].real() << ", expect res: " << plain[i] << std::endl;
    // }
    double err = 0.;
    // std::cout << "plain.size:" << plain.size() << std::endl;
    for (int i = 0; i < int(plain.size()); i++) {
        err += fabs(plain[i] - Result[i].real());
        error_vec.push_back(fabs(plain[i] - plaintextDec->GetCKKSPackedValue()[i].real()));
        if (fabs(plain[i] - Result[i].real()) > 0.5)
            throw std::runtime_error(" err>0.5  ");
    }
    // int maxPosition = max_element(error_vec.begin(), error_vec.end()) - error_vec.begin();
    // std::cout << "error:  " << err / int(plain.size()) << " ~ 2^" << std::log2(err / int(plain.size())) << std::endl;
    return err / plain.size();
}

double error_estimate_modular(std::vector<double> plain, std::vector<Ciphertext<lbcrypto::DCRTPoly>>& res,
                              lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey, size_t num_slots) {
    auto cc = res[0]->GetCryptoContext();
    std::vector<double> error_vec;
    std::vector<double> merge(num_slots, 0);
    std::cout << "error_estimate_modular " << std::endl;
    Plaintext plaintextDec0, plaintextDec1;
    cc->Decrypt(privateKey, res[0], &plaintextDec0);
    std::cout << "Decrypt finish " << std::endl;
    std::vector<std::complex<double>> Result0 = plaintextDec0->GetCKKSPackedValue();
    cc->Decrypt(privateKey, res[1], &plaintextDec1);
    std::cout << "Decrypt finish " << std::endl;
    std::vector<std::complex<double>> Result1 = plaintextDec1->GetCKKSPackedValue();
    std::cout << "merge   " << std::endl;
    for (int j = 0; j < int(Result0.size()); j++) {
        merge[j] += round(Result0[j].real()) + round(Result1[j].real()) * pow(2, 8);
    }
    std::cout << "merge finish " << std::endl;

    // for (int i = 0; i < int(res.size()); i++) {
    //     Plaintext plaintextDec;
    //     std::cout << "Decrypt " << std::endl;
    //     cc->Decrypt(privateKey, res[i], &plaintextDec);
    //     std::cout << "Decrypt finish " << std::endl;
    //     std::vector<std::complex<double>> Result = plaintextDec->GetCKKSPackedValue();
    //     for (int j = 0; j < int(Result.size()); j++) {
    //         merge[j] += round(Result[j].real()) * pow(2, i);
    //     }
    //     std::cout << "merge finish " << std::endl;
    // }

    // std::cout << "Estimated precision: " << plaintextDec->GetLogPrecision() << std::endl;

    // for (int i = 0; i < int(num_slots); i++) {
    //     std::cout << "Result dec: " << Result[i].real() << ", expect res: " << plain[i] << std::endl;
    // }
    // double err = 0.;
    // std::cout << "merge.size:" << merge.size() << std::endl;
    // std::cout << "plain.size:" << plain.size() << std::endl;
    // for (int i = 0; i < int(plain.size()); i++) {
    //     err += fabs(plain[i] - merge[i]);
    //     error_vec.push_back(fabs(plain[i] - merge[i]));

    //     if (fabs(plain[i] - merge[i]) > 0.5)
    //         std::cout << "error! plain:  " << plain[i] << "  merge: " << merge[i] << std::endl;
    //     // throw std::runtime_error(" err>0.5  ");
    // }

    double err = 0.;
    std::cout << "merge.size:" << merge.size() << std::endl;
    std::cout << "plain.size:" << plain.size() << std::endl;
    for (int i = 0; i < int(plain.size()); i++) {
        err += fabs((int(plain[i]) >> 8) - Result1[i].real());
        error_vec.push_back(fabs(int(plain[i]) >> 8) - Result1[i].real());

        if (error_vec[i] > 0.5)
            std::cout << "error! plain:  " << plain[i] << "  merge: " << merge[i] << std::endl;
        // throw std::runtime_error(" err>0.5  ");
    }
    // int maxPosition = max_element(error_vec.begin(), error_vec.end()) - error_vec.begin();
    // std::cout << "error:  " << err / int(plain.size()) << " ~ 2^" << std::log2(err / int(plain.size())) << std::endl;
    return err / plain.size();
}

void bitonic_sort(int plain_bits, int num_slots) {
    std::cout << "\nsort: " << plain_bits << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits               = 78;
    usint firstMod               = 89;
#else
    usint dcrtBits = 59;
    usint firstMod = 60;
#endif
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetFirstModSize(firstMod);
    std::vector<uint32_t> levelBudget = {4, 4};
    // SecretKeyDist secretKeyDist       = UNIFORM_TERNARY;
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    // Choosing a higher degree yields better precision, but a longer runtime.
    uint32_t polyDegree = 495;  //

    // The multiplicative depth depends on the polynomial degree.
    // See the FUNCTION_EVALUATION.md file for a table mapping polynomial degrees to multiplicative depths.
    // uint32_t levelsAvailableAfterBootstrap = 35;  //35
    // uint32_t levelsbeforeBootstrap         = 15;
    uint32_t numIterations = 2;
    // usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist) +
    //                   (numIterations - 1) + levelsbeforeBootstrap;
    uint32_t levelsAvailableAfterBootstrap = 25;
    usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    std::cout << "polyDegree " << polyDegree << std::endl;
    std::cout << "CyclotomicOrder " << cc->GetCyclotomicOrder() << std::endl;
    std::cout << "RingDimension " << cc->GetRingDimension() << std::endl;
    std::cout << "Modulus " << cc->GetModulus() << std::endl;
    std::cout << "PlaintextModulus " << cc->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(FHE);
    cc->Enable(ADVANCEDSHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.

    double precision = (1 << (plain_bits - 1)) - 1;
    std::cout << "precision: " << precision << std::endl;

    double lowerBound = -precision - 10;
    double upperBound = precision + 10;
    // double bound      = 3;
    auto keyPair  = cc->KeyGen();
    usint ringDim = cc->GetRingDimension();
    int length    = ringDim / 2;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, precision);
    std::vector<double> input(length, 0);
    const std::vector<DCRTPoly>& ckks_pk = keyPair.publicKey->GetPublicElements();
    std::cout << "Moduli chain of pk: " << std::endl;
    print_moduli_chain(ckks_pk[0]);
    for (int i = 0; i < length; i++) {
        input[i] = message(engine);
    }
    // input[0] = 80;
    // input[1] = 10;
    // input[2] = 40;
    // input[3] = 60;
    // input[4] = 5;
    // input[5] = 81;
    // input[6] = 189;
    // input[7] = 165;
    for (int i = 0; i < 8; i++) {
        std::cout << input[i] << std::endl;
    }
    std::vector<std::vector<double>> plain_sort_vec;
    std::vector<std::vector<std::vector<double>>> matrix = bitonic_test_plain(input, plain_sort_vec, num_slots);
    // for (int i = 0; i < plain_sort_vec[0].size(); i++) {
    //     std::cout << plain_sort_vec[0][i] << ",";
    // }
    std::cout << "ringDim: " << ringDim << std::endl;
    cc->EvalMultKeyGen(keyPair.secretKey);
    std::vector<int> rotstep;
    for (int i = 1; i < num_slots; i *= 2) {
        rotstep.push_back(i);
        rotstep.push_back(-i);
    }

    cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
    cc->EvalBootstrapSetup(levelBudget);
    cc->EvalBootstrapKeyGen(keyPair.secretKey, length);

    std::vector<double> coefficients = EvalChebyshevCoefficients([](double x) -> double { return ((x > 0) ? 1 : 0); },
                                                                 lowerBound, upperBound, polyDegree);
    Plaintext plain                  = cc->MakeCKKSPackedPlaintext(input);
    auto ciphertext_unsort           = cc->Encrypt(keyPair.publicKey, plain);

    std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
    std::cout << "number of levels fresh: " << multDepth - ciphertext_unsort->GetLevel() << std::endl;
    start             = std::chrono::system_clock::now();
    double time_total = 0;
    for (int stage = 0; num_slots >> stage + 1 > 0; stage++) {
        for (int part = 0; stage - part >= 0; part++) {
            std::cout << " stage " << stage << " part " << part << std::endl;
            Ciphertext<lbcrypto::DCRTPoly> comp_res;
            Ciphertext<lbcrypto::DCRTPoly> sort_res;
            // std::cout << "bitonic_comp " << std::endl;
            start_comp = std::chrono::system_clock::now();
            bitonic_comp(stage, part, length, ciphertext_unsort, comp_res, precision, coefficients, keyPair.secretKey);
            end_comp = std::chrono::system_clock::now();
            // comp_res = cc->EvalBootstrap(comp_res);
            // Plaintext plaintextDec;
            // cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec);
            // std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
            // for (int i = 0; i < 2 * int(num_slots); i++) {
            //     std::cout << " comp_res: " << enc[i].real() << std::endl;
            // }
            // std::cout << "number of levels remaining before bitonic_swap: " << multDepth - ciphertext_unsort->GetLevel()
            //           << std::endl;
            std::cout << "bitonic_swap " << std::endl;
            start_swap = std::chrono::system_clock::now();
            bitonic_swap(stage, part, length, comp_res, ciphertext_unsort, sort_res, keyPair.secretKey);
            end_swap             = std::chrono::system_clock::now();
            double err_afterswap = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], sort_res,
                                                  keyPair.secretKey, 2 * num_slots);
            // std::cout << "number of levels remaining after bitonic_swap: " << multDepth - ciphertext_unsort->GetLevel()
            //           << std::endl;
            std::cout << "error after bitonic_swap is " << err_afterswap << " ~ 2^" << std::log2(err_afterswap)
                      << std::endl;

            if ((multDepth - sort_res->GetLevel()) < 23) {  //23
                sort_res = cc->EvalMult(sort_res, 1 / precision);
                std::cout << "boot " << std::endl;
                // int bootprecision = itboot(sort_res, keyPair.secretKey);
                // std::cout << " sort_res bootprecision: " << bootprecision << std::endl;
                start_boot        = std::chrono::system_clock::now();
                ciphertext_unsort = cc->EvalBootstrap(sort_res, numIterations, 8);
                // ciphertext_unsort = cc->EvalBootstrap(sort_res);
                ciphertext_unsort = cc->EvalMult(ciphertext_unsort, precision);
                end_boot          = std::chrono::system_clock::now();
                std::cout << "number of levels remaining  after boot: " << multDepth - ciphertext_unsort->GetLevel()
                          << std::endl;
            }
            else {
                start_boot        = std::chrono::system_clock::now();
                ciphertext_unsort = sort_res;
                end_boot          = std::chrono::system_clock::now();
                std::cout << "number of levels remaining : " << multDepth - ciphertext_unsort->GetLevel() << std::endl;
            }

            double err_afterBoot = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], ciphertext_unsort,
                                                  keyPair.secretKey, 2 * num_slots);

            std::cout << "error afterBoot in stage " << stage << " part " << part << " is " << err_afterBoot << " ~ 2^"
                      << std::log2(err_afterBoot) << std::endl;
            time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
            std::cout << "boot time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count())
                      << "ms" << std::endl;
            std::cout << "swap time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count())
                      << "ms" << std::endl;
            std::cout << "comp time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count())
                      << "ms" << std::endl
                      << std::endl;
        }

        std::cout << "stage " << stage << " sort total time: " << time_total << "ms" << std::endl;
        std::cout << "stage " << stage << " sort amortize time: " << time_total / double(length) * (1 << (stage + 1))
                  << "ms" << std::endl
                  << std::endl
                  << std::endl
                  << std::endl;
    }
    end = std::chrono::system_clock::now();
    Plaintext plaintextDec;
    cc->Decrypt(keyPair.secretKey, ciphertext_unsort, &plaintextDec);
    std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
    std::cout << " sortres " << std::endl;
    for (int i = 0; i < int(num_slots); i++) {
        std::cout << enc[i].real() << std::endl;
    }
    std::cout << length << "slots sort total time: " << time_total << "ms" << std::endl;
    std::cout << length << "slots sort amortize time: " << time_total / double(length) * num_slots << "ms" << std::endl
              << std::endl;
}

void bitonic_sort_small(int plain_bits, int num_slots) {
    std::cout << "\nsort: " << plain_bits << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits               = 78;
    usint firstMod               = 89;
#else
    usint dcrtBits = 59;
    usint firstMod = 60;
#endif
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetFirstModSize(firstMod);
    std::vector<uint32_t> levelBudget = {4, 4};
    // SecretKeyDist secretKeyDist       = UNIFORM_TERNARY;
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    // Choosing a higher degree yields better precision, but a longer runtime.
    uint32_t polyDegree = 119;  //

    // The multiplicative depth depends on the polynomial degree.
    // See the FUNCTION_EVALUATION.md file for a table mapping polynomial degrees to multiplicative depths.
    // uint32_t levelsAvailableAfterBootstrap = 35;  //35
    // uint32_t levelsbeforeBootstrap         = 15;
    uint32_t numIterations = 2;
    // usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist) +
    //                   (numIterations - 1) + levelsbeforeBootstrap;
    uint32_t levelsAvailableAfterBootstrap = 25;
    usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    std::cout << "polyDegree " << polyDegree << std::endl;
    std::cout << "CyclotomicOrder " << cc->GetCyclotomicOrder() << std::endl;
    std::cout << "RingDimension " << cc->GetRingDimension() << std::endl;
    std::cout << "Modulus " << cc->GetModulus() << std::endl;
    std::cout << "PlaintextModulus " << cc->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(FHE);
    cc->Enable(ADVANCEDSHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.

    double precision = (1 << (plain_bits - 1)) - 1;
    std::cout << "precision: " << precision << std::endl;

    double lowerBound = -precision - 10;
    double upperBound = precision + 10;
    // double bound      = 3;
    auto keyPair  = cc->KeyGen();
    usint ringDim = cc->GetRingDimension();
    int length    = ringDim / 2;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, precision);
    std::vector<double> input(length, 0);
    const std::vector<DCRTPoly>& ckks_pk = keyPair.publicKey->GetPublicElements();
    std::cout << "Moduli chain of pk: " << std::endl;
    print_moduli_chain(ckks_pk[0]);
    for (int i = 0; i < length; i++) {
        input[i] = message(engine);
    }
    // input[0] = 80;
    // input[1] = 10;
    // input[2] = 40;
    // input[3] = 60;
    // input[4] = 5;
    // input[5] = 81;
    // input[6] = 189;
    // input[7] = 165;
    for (int i = 0; i < 8; i++) {
        std::cout << input[i] << std::endl;
    }
    std::vector<std::vector<double>> plain_sort_vec;
    std::vector<std::vector<std::vector<double>>> matrix = bitonic_test_plain(input, plain_sort_vec, num_slots);
    // for (int i = 0; i < plain_sort_vec[0].size(); i++) {
    //     std::cout << plain_sort_vec[0][i] << ",";
    // }
    std::cout << "ringDim: " << ringDim << std::endl;
    cc->EvalMultKeyGen(keyPair.secretKey);
    std::vector<int> rotstep;
    for (int i = 1; i < num_slots; i *= 2) {
        rotstep.push_back(i);
        rotstep.push_back(-i);
    }

    cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
    cc->EvalBootstrapSetup(levelBudget);
    cc->EvalBootstrapKeyGen(keyPair.secretKey, length);

    std::vector<double> coefficients = EvalChebyshevCoefficients([](double x) -> double { return ((x > 0) ? 1 : 0); },
                                                                 lowerBound, upperBound, polyDegree);
    Plaintext plain                  = cc->MakeCKKSPackedPlaintext(input);
    auto ciphertext_unsort           = cc->Encrypt(keyPair.publicKey, plain);

    std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
    std::cout << "number of levels fresh: " << multDepth - ciphertext_unsort->GetLevel() << std::endl;
    start             = std::chrono::system_clock::now();
    double time_total = 0;
    for (int stage = 0; num_slots >> stage + 1 > 0; stage++) {
        for (int part = 0; stage - part >= 0; part++) {
            std::cout << " stage " << stage << " part " << part << std::endl;
            Ciphertext<lbcrypto::DCRTPoly> comp_res;
            Ciphertext<lbcrypto::DCRTPoly> sort_res;
            // std::cout << "bitonic_comp " << std::endl;
            start_comp = std::chrono::system_clock::now();
            bitonic_comp(stage, part, length, ciphertext_unsort, comp_res, precision, coefficients, keyPair.secretKey);
            end_comp = std::chrono::system_clock::now();
            // comp_res = cc->EvalBootstrap(comp_res);
            // Plaintext plaintextDec;
            // cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec);
            // std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
            // for (int i = 0; i < 2 * int(num_slots); i++) {
            //     std::cout << " comp_res: " << enc[i].real() << std::endl;
            // }
            // std::cout << "number of levels remaining before bitonic_swap: " << multDepth - ciphertext_unsort->GetLevel()
            //           << std::endl;
            std::cout << "bitonic_swap " << std::endl;
            start_swap = std::chrono::system_clock::now();
            bitonic_swap(stage, part, length, comp_res, ciphertext_unsort, sort_res, keyPair.secretKey);
            end_swap             = std::chrono::system_clock::now();
            double err_afterswap = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], sort_res,
                                                  keyPair.secretKey, 2 * num_slots);
            // std::cout << "number of levels remaining after bitonic_swap: " << multDepth - ciphertext_unsort->GetLevel()
            //           << std::endl;
            std::cout << "error after bitonic_swap is " << err_afterswap << " ~ 2^" << std::log2(err_afterswap)
                      << std::endl;

            if ((multDepth - sort_res->GetLevel()) < 23) {  //23
                sort_res = cc->EvalMult(sort_res, 1 / precision);
                std::cout << "boot " << std::endl;
                // int bootprecision = itboot(sort_res, keyPair.secretKey);
                // std::cout << " sort_res bootprecision: " << bootprecision << std::endl;
                start_boot        = std::chrono::system_clock::now();
                ciphertext_unsort = cc->EvalBootstrap(sort_res, numIterations, 8);
                // ciphertext_unsort = cc->EvalBootstrap(sort_res);
                ciphertext_unsort = cc->EvalMult(ciphertext_unsort, precision);
                end_boot          = std::chrono::system_clock::now();
                std::cout << "number of levels remaining  after boot: " << multDepth - ciphertext_unsort->GetLevel()
                          << std::endl;
            }
            else {
                start_boot        = std::chrono::system_clock::now();
                ciphertext_unsort = sort_res;
                end_boot          = std::chrono::system_clock::now();
                std::cout << "number of levels remaining : " << multDepth - ciphertext_unsort->GetLevel() << std::endl;
            }

            double err_afterBoot = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], ciphertext_unsort,
                                                  keyPair.secretKey, 2 * num_slots);

            std::cout << "error afterBoot in stage " << stage << " part " << part << " is " << err_afterBoot << " ~ 2^"
                      << std::log2(err_afterBoot) << std::endl;
            time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
            std::cout << "boot time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count())
                      << "ms" << std::endl;
            std::cout << "swap time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count())
                      << "ms" << std::endl;
            std::cout << "comp time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count())
                      << "ms" << std::endl
                      << std::endl;
        }

        std::cout << "stage " << stage << " sort total time: " << time_total << "ms" << std::endl;
        std::cout << "stage " << stage << " sort amortize time: " << time_total / double(length) * (1 << (stage + 1))
                  << "ms" << std::endl
                  << std::endl
                  << std::endl
                  << std::endl;
    }
    end = std::chrono::system_clock::now();
    Plaintext plaintextDec;
    cc->Decrypt(keyPair.secretKey, ciphertext_unsort, &plaintextDec);
    std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
    std::cout << " sortres " << std::endl;
    for (int i = 0; i < int(num_slots); i++) {
        std::cout << enc[i].real() << std::endl;
    }
    std::cout << length << "slots sort total time: " << time_total << "ms" << std::endl;
    std::cout << length << "slots sort amortize time: " << time_total / double(length) * num_slots << "ms" << std::endl
              << std::endl;
}

std::vector<int> dot_product(std::vector<int>& A, std::vector<int>& B) {
    std::vector<int> result(B.size());
    for (int i = 0; i < A.size(); i++) {
        result[i] += A[i] * B[i];
    }
    return result;
}
std::vector<int> vadd(std::vector<int>& A, std::vector<int>& B) {
    std::vector<int> result(B.size());
    for (int i = 0; i < A.size(); i++) {
        result[i] += A[i] + B[i];
    }
    return result;
}

void cross_rotate(int N, int l, int T) {
    std::vector<std::vector<int>> big_vec;
    std::vector<std::vector<int>> big_vec_rot;
    std::vector<int> test;
    // int N     = 8;
    int baby  = T % N;
    int giant = floor(T / N);
    for (int i = 0; i < l; i++) {
        std::vector<int> small_vec;
        for (int j = 0; j < N; j++) {
            small_vec.push_back(i * N + j);
            test.push_back(i * N + j);
        }
        big_vec.push_back(small_vec);
    }
    for (int i = 0; i < test.size(); i++)
        std::cout << test[i] << ", " << std::endl;
    std::vector<std::vector<int>> big_vec_prerot;
    for (int i = 0; i < l; i++) {
        big_vec_prerot.push_back(big_vec[(i + giant) % l]);
    }
    std::vector<int> m0(N, 0);
    std::vector<int> m1(N, 1);
    for (int i = 0; i < N - baby; i++) {
        m0[i] = 1;
        m1[i] = 0;
    }
    for (int i = 0; i < l - 1; i++) {
        std::vector<int> v0(N), v1(N), rm0(N), rm1(N);
        std::rotate_copy(big_vec_prerot[i].begin(), big_vec_prerot[i].begin() + baby, big_vec_prerot[i].end(),
                         v0.begin());
        std::rotate_copy(big_vec_prerot[i + 1].begin(), big_vec_prerot[i + 1].begin() + baby,
                         big_vec_prerot[i + 1].end(), v1.begin());
        rm0 = dot_product(v0, m0);
        rm1 = dot_product(v1, m1);
        big_vec_rot.push_back(vadd(rm0, rm1));
    }
    std::vector<int> v0(N), v1(N), rm0, rm1;
    std::rotate_copy(big_vec_prerot[l - 1].begin(), big_vec_prerot[l - 1].begin() + baby, big_vec_prerot[l - 1].end(),
                     v0.begin());
    std::rotate_copy(big_vec_prerot[0].begin(), big_vec_prerot[0].begin() + baby, big_vec_prerot[0].end(), v1.begin());
    rm0 = dot_product(v0, m0);
    rm1 = dot_product(v1, m1);
    big_vec_rot.push_back(vadd(rm0, rm1));
    std::rotate(test.begin(), test.begin() + T, test.end());
    for (int i = 0; i < l; i++) {
        for (int j = 0; j < N; j++) {
            // if (big_vec_rot[i][j] != test[i * l + j])
            std::cout << "i: " << i << ", j: " << j << ", big_vec_rot: " << big_vec_rot[i][j]
                      << ", test:" << test[i * N + j] << std::endl;
        }
    }
}
// N slots number, l vector number, T rotate step
std::vector<Ciphertext<lbcrypto::DCRTPoly>> gen_rotate(int N, int l, int T,
                                                       std::vector<Ciphertext<lbcrypto::DCRTPoly>>& big_vec,
                                                       std::vector<double>& test,
                                                       lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey) {
    std::vector<Ciphertext<lbcrypto::DCRTPoly>> big_vec_rot;

    // int N     = 8;
    int baby  = T % N;
    auto cc   = big_vec[0]->GetCryptoContext();
    int giant = floor(T / N);
    for (int i = 0; i < l; i++) {
        Plaintext plaintextDec;
        cc->Decrypt(privateKey, big_vec[i], &plaintextDec);
        std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
        for (int j = 0; j < 20; j++) {
            std::cout << enc[j].real() << std::endl;
        }
        for (int j = 0; j < N; j++) {
            if (round(enc[j].real()) != test[i * N + j])
                std::cout << "i: " << i << ", j: " << j << ", big_vec_rot: " << enc[j].real()
                          << ", test:" << test[i * N + j] << std::endl;
        }
    }
    std::rotate(test.begin(), test.begin() + T, test.end());
    std::vector<Ciphertext<lbcrypto::DCRTPoly>> big_vec_prerot;
    for (int i = 0; i < l; i++) {
        big_vec_prerot.push_back(big_vec[(i + giant) % l]);
    }
    std::vector<double> mask0(N, 0);
    std::vector<double> mask1(N, 1);
    for (int i = 0; i < N - baby; i++) {
        mask0[i] = 1;
        mask1[i] = 0;
    }
    Plaintext m0 = cc->MakeCKKSPackedPlaintext(mask0);
    Plaintext m1 = cc->MakeCKKSPackedPlaintext(mask1);
    for (int i = 0; i < l - 1; i++) {
        auto v0  = cc->EvalRotate(big_vec_prerot[i], baby);
        auto v1  = cc->EvalRotate(big_vec_prerot[i + 1], baby);
        auto rm0 = cc->EvalMult(v0, m0);
        auto rm1 = cc->EvalMult(v1, m1);
        big_vec_rot.push_back(cc->EvalAdd(rm0, rm1));
    }

    auto v0  = cc->EvalRotate(big_vec_prerot[l - 1], baby);
    auto v1  = cc->EvalRotate(big_vec_prerot[0], baby);
    auto rm0 = cc->EvalMult(v0, m0);
    auto rm1 = cc->EvalMult(v1, m1);
    big_vec_rot.push_back(cc->EvalAdd(rm0, rm1));
    Plaintext plaintextDec;
    for (int i = 0; i < l; i++) {
        cc->Decrypt(privateKey, big_vec_rot[i], &plaintextDec);
        std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
        std::cout << " sortres " << std::endl;
        for (int j = 0; j < 20; j++) {
            std::cout << enc[j].real() << std::endl;
        }
        for (int j = 0; j < N; j++) {
            if (round(enc[j].real()) != test[i * N + j])
                std::cout << "i: " << i << ", j: " << j << ", big_vec_rot: " << enc[i].real()
                          << ", test:" << test[i * N + j] << std::endl;
        }
    }
    std::cout << " finish " << std::endl;
}

void sync_test_big(int num_column, int plain_bits, int num_slots) {
    std::cout << "65536 slots " << num_column << " column sync test " << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits               = 78;
    usint firstMod               = 89;
#else
    usint dcrtBits = 40;
    usint firstMod = 60;
#endif
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetFirstModSize(firstMod);

    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetRingDim(2 * num_slots);
    usint multDepth = (log2(num_slots) + 1) * log2(num_slots) / 2 + 1;
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(FHE);
    cc->Enable(ADVANCEDSHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.

    double precision = (1 << (plain_bits - 1)) - 1;
    std::cout << "precision: " << precision << std::endl;

    double lowerBound = -precision - 10;
    double upperBound = precision + 10;
    // double bound      = 3;
    auto keyPair  = cc->KeyGen();
    usint ringDim = cc->GetRingDimension();
    int length    = ringDim / 2;
    std::cout << "ringDim: " << ringDim << std::endl;
    cc->EvalMultKeyGen(keyPair.secretKey);
    std::vector<int> rotstep;
    for (int i = 1; i < num_slots; i *= 2) {
        rotstep.push_back(i);
        rotstep.push_back(-i);
    }

    cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, precision);
    std::vector<std::vector<double>> column;
    std::vector<double> input(length, 0);
    for (int i = 0; i < length; i++) {
        input[i] = message(engine);
    }

    std::vector<std::vector<double>> plain_sort_vec;
    std::vector<std::vector<std::vector<double>>> matrix = bitonic_test_plain(input, plain_sort_vec, num_slots);
    std::vector<Ciphertext<lbcrypto::DCRTPoly>> swap_matrix;

    std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;

    double time_total = 0;
    for (int col = 0; col < num_column; col++) {
        Plaintext plain        = cc->MakeCKKSPackedPlaintext(input);
        auto ciphertext_unsort = cc->Encrypt(keyPair.publicKey, plain);
        int part_count         = 0;
        std::cout << "number of levels fresh: " << multDepth - ciphertext_unsort->GetLevel() << std::endl;
        for (int stage = 0; num_slots >> stage + 1 > 0; stage++) {
            for (int part = 0; stage - part >= 0; part++) {
                int block      = 1 << (stage + 1 - part);
                int part_count = part + ((1 + stage) * stage) / 2;
                std::vector<Ciphertext<lbcrypto::DCRTPoly>> swap_matrix;
                for (int i = 0; i < matrix[part_count].size(); i++) {
                    Plaintext plain_swap_vec                = cc->MakeCKKSPackedPlaintext(matrix[part_count][i]);
                    Ciphertext<lbcrypto::DCRTPoly> swap_vec = cc->Encrypt(keyPair.publicKey, plain_swap_vec);
                    swap_matrix.push_back(swap_vec);
                }
                start         = std::chrono::system_clock::now();
                auto ct_upper = cc->EvalRotate(ciphertext_unsort, -block / 2);
                auto ct_lower = cc->EvalRotate(ciphertext_unsort, block / 2);

                std::vector<Ciphertext<lbcrypto::DCRTPoly>> res_martrix;
                res_martrix.push_back(cc->EvalMult(ct_upper, swap_matrix[1]));
                res_martrix.push_back(cc->EvalMult(ciphertext_unsort, swap_matrix[0]));
                res_martrix.push_back(cc->EvalMult(ct_lower, swap_matrix[2]));

                auto temp         = cc->EvalAdd(res_martrix[0], res_martrix[1]);
                ciphertext_unsort = cc->EvalAdd(temp, res_martrix[2]);
                end               = std::chrono::system_clock::now();
                time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
            }
        }
    }

    std::cout << plain_bits << " bits " << num_column << " column sync time: " << time_total * num_slots / length
              << "ms" << std::endl;
}
//4096
void sync_test_small(int plain_bits, int num_slots) {
    std::cout << num_slots << " records sync test " << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecurityLevel(HEStd_128_classic);
    uint32_t scaleModSize = 50;
    parameters.SetScalingModSize(scaleModSize);

    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    // if (num_slots > 4096)
    //     parameters.SetRingDim(2 * num_slots);

    usint multDepth = (log2(num_slots) + 1) * log2(num_slots) / 2 + 1;
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(FHE);
    cc->Enable(ADVANCEDSHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.

    double precision = (1 << (plain_bits - 1)) - 1;

    double lowerBound = -precision - 10;
    double upperBound = precision + 10;
    // double bound      = 3;
    auto keyPair  = cc->KeyGen();
    usint ringDim = cc->GetRingDimension();
    int length    = ringDim / 2;
    std::cout << "ringDim: " << ringDim << std::endl;
    cc->EvalMultKeyGen(keyPair.secretKey);
    std::vector<int> rotstep;
    for (int i = 1; i < num_slots; i *= 2) {
        rotstep.push_back(i);
        rotstep.push_back(-i);
    }

    cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, precision);
    std::vector<std::vector<double>> column;
    std::vector<double> input(length, 0);
    for (int i = 0; i < length; i++) {
        input[i] = message(engine);
    }

    std::vector<std::vector<double>> plain_sort_vec;
    std::vector<std::vector<std::vector<double>>> matrix = bitonic_test_plain(input, plain_sort_vec, num_slots);
    std::vector<Ciphertext<lbcrypto::DCRTPoly>> swap_matrix;

    Plaintext plain        = cc->MakeCKKSPackedPlaintext(input);
    auto ciphertext_unsort = cc->Encrypt(keyPair.publicKey, plain);
    int part_count         = 0;
    std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
    std::cout << "number of levels fresh: " << multDepth - ciphertext_unsort->GetLevel() << std::endl;
    double time_total = 0;

    for (int stage = 0; num_slots >> stage + 1 > 0; stage++) {
        for (int part = 0; stage - part >= 0; part++) {
            int block      = 1 << (stage + 1 - part);
            int part_count = part + ((1 + stage) * stage) / 2;
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> swap_matrix;
            for (int i = 0; i < matrix[part_count].size(); i++) {
                Plaintext plain_swap_vec                = cc->MakeCKKSPackedPlaintext(matrix[part_count][i]);
                Ciphertext<lbcrypto::DCRTPoly> swap_vec = cc->Encrypt(keyPair.publicKey, plain_swap_vec);
                swap_matrix.push_back(swap_vec);
            }
            start         = std::chrono::system_clock::now();
            auto ct_upper = cc->EvalRotate(ciphertext_unsort, -block / 2);
            auto ct_lower = cc->EvalRotate(ciphertext_unsort, block / 2);

            std::vector<Ciphertext<lbcrypto::DCRTPoly>> res_martrix;
            res_martrix.push_back(cc->EvalMult(ct_upper, swap_matrix[1]));
            res_martrix.push_back(cc->EvalMult(ciphertext_unsort, swap_matrix[0]));
            res_martrix.push_back(cc->EvalMult(ct_lower, swap_matrix[2]));

            auto temp         = cc->EvalAdd(res_martrix[0], res_martrix[1]);
            ciphertext_unsort = cc->EvalAdd(temp, res_martrix[2]);
            end               = std::chrono::system_clock::now();
            time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
        }
    }

    std::cout << " sync time: " << time_total * num_slots / length << "ms" << std::endl;
}

void bitonic_sort_modular(int plain_bits, int blocks, int num_slots) {
    std::cout << "\nmodular sort: " << plain_bits << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits               = 78;
    usint firstMod               = 89;
#else
    usint dcrtBits = 59;
    usint firstMod = 60;
#endif
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetFirstModSize(firstMod);
    std::vector<uint32_t> levelBudget = {4, 4};
    // SecretKeyDist secretKeyDist       = UNIFORM_TERNARY;
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    // Choosing a higher degree yields better precision, but a longer runtime.
    uint32_t polyDegree = 247;  //

    // The multiplicative depth depends on the polynomial degree.
    // See the FUNCTION_EVALUATION.md file for a table mapping polynomial degrees to multiplicative depths.
    // uint32_t levelsAvailableAfterBootstrap = 35;  //35
    // uint32_t levelsbeforeBootstrap         = 15;
    uint32_t numIterations = 2;
    // usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist) +
    //                   (numIterations - 1) + levelsbeforeBootstrap;
    uint32_t levelsAvailableAfterBootstrap = 25;
    usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    std::cout << "CyclotomicOrder " << cc->GetCyclotomicOrder() << std::endl;
    std::cout << "RingDimension " << cc->GetRingDimension() << std::endl;
    // std::cout << "Modulus " << cc->GetModulus() << std::endl;
    std::cout << "PlaintextModulus " << cc->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(FHE);
    cc->Enable(ADVANCEDSHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.

    double precision = (1 << (plain_bits - 1)) - 1;
    std::cout << "precision: " << precision << std::endl;

    double lowerBound = -precision - 10;
    double upperBound = precision + 10;
    // double bound      = 3;
    auto keyPair  = cc->KeyGen();
    usint ringDim = cc->GetRingDimension();
    int length    = ringDim / 2;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, precision);
    std::vector<double> input_big(length, 0);
    std::vector<std::vector<double>> input_vec;
    for (int j = 0; j < blocks; j++) {
        std::vector<double> input_0(length, 0);
        for (int i = 0; i < length; i++) {
            input_0[i] = message(engine);
        }
        input_vec.push_back(input_0);
    }
    std::cout << "input_vec  " << std::endl;

    for (int i = 0; i < length; i++) {
        for (int j = 0; j < blocks; j++) {
            input_big[i] += double(int(input_vec[j][i]) << (plain_bits * j));
        }
        // std::cout << input_big[i] << std::endl;
    }
    std::cout << "input_big  " << std::endl;
    // input[0] = 80;
    // input[1] = 10;
    // input[2] = 40;
    // input[3] = 60;
    // input[4] = 5;
    // input[5] = 81;
    // input[6] = 189;
    // input[7] = 165;
    for (int i = 0; i < 16; i++) {
        std::cout << input_big[i] << std::endl;
    }
    std::vector<std::vector<double>> plain_sort_vec;
    std::vector<std::vector<std::vector<double>>> matrix = bitonic_test_plain(input_big, plain_sort_vec, length);

    std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_unsort;
    // for (int i = 0; i < plain_sort_vec[0].size(); i++) {
    //     std::cout << plain_sort_vec[0][i] << ",";
    // }
    std::cout << "ringDim: " << ringDim << std::endl;
    cc->EvalMultKeyGen(keyPair.secretKey);
    std::vector<int> rotstep;
    for (int i = 1; i < num_slots; i *= 2) {
        rotstep.push_back(i);
        rotstep.push_back(-i);
    }

    cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
    cc->EvalBootstrapSetup(levelBudget);
    cc->EvalBootstrapKeyGen(keyPair.secretKey, length);

    std::vector<double> coefficients = EvalChebyshevCoefficients([](double x) -> double { return ((x > 0) ? 1 : 0); },
                                                                 lowerBound, upperBound, polyDegree);
    std::vector<Plaintext> plain_big;
    // std::vector<Plaintext> plain_big;
    for (int i = 0; i < blocks; i++) {
        Plaintext plain_0        = cc->MakeCKKSPackedPlaintext(input_vec[i]);
        auto ciphertext_unsort_0 = cc->Encrypt(keyPair.publicKey, plain_0);
        ciphertext_unsort.push_back(ciphertext_unsort_0);
    }
    // double test_error = error_estimate_modular(input_big, ciphertext_unsort, keyPair.secretKey, length);
    // auto ciphertext_unsort_0 = cc->Encrypt(keyPair.publicKey, plain_0);
    // Plaintext plain_1        = cc->MakeCKKSPackedPlaintext(input_1);
    // auto ciphertext_unsort_1 = cc->Encrypt(keyPair.publicKey, plain_1);
    // ciphertext_unsort.push_back(ciphertext_unsort_0);
    // ciphertext_unsort.push_back(ciphertext_unsort_1);
    const std::vector<DCRTPoly>& ckks_pk = keyPair.publicKey->GetPublicElements();
    std::cout << "Moduli chain of pk: " << std::endl;
    print_moduli_chain(ckks_pk[0]);
    std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
    std::cout << "number of levels fresh: " << multDepth - ciphertext_unsort[0]->GetLevel() << std::endl;
    start             = std::chrono::system_clock::now();
    double time_total = 0;
    for (int stage = 0; num_slots >> stage + 1 > 0; stage++) {
        for (int part = 0; stage - part >= 0; part++) {
            std::cout << " stage " << stage << " part " << part << std::endl;
            Ciphertext<lbcrypto::DCRTPoly> comp_res;
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> sort_res_vec;
            // std::cout << "bitonic_comp " << std::endl;
            start_comp = std::chrono::system_clock::now();
            bitonic_comp_modular(stage, part, length, ciphertext_unsort, comp_res, precision, coefficients,
                                 keyPair.secretKey);
            end_comp = std::chrono::system_clock::now();
            std::cout << "bitonic_comp_modular finish" << std::endl;
            // comp_res = cc->EvalBootstrap(comp_res);
            // Plaintext plaintextDec;
            // cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec);
            // std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
            // for (int i = 0; i < 2 * int(num_slots); i++) {
            //     std::cout << " comp_res: " << enc[i].real() << std::endl;
            // }
            // std::cout << "number of levels remaining before bitonic_swap: " << multDepth - ciphertext_unsort->GetLevel()
            //           << std::endl;
            // std::cout << "bitonic_swap " << std::endl;
            start_swap = std::chrono::system_clock::now();

            for (int i = 0; i < ciphertext_unsort.size(); i++) {
                Ciphertext<lbcrypto::DCRTPoly> sort_res;
                bitonic_swap(stage, part, length, comp_res, ciphertext_unsort[i], sort_res, keyPair.secretKey);
                sort_res_vec.push_back(sort_res);
            }
            end_swap = std::chrono::system_clock::now();
            std::cout << "bitonic_swap finish" << std::endl;
            // double err_afterswap = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], sort_res[0],
            //                                       keyPair.secretKey, 2 * num_slots);
            // std::cout << "number of levels remaining after bitonic_swap: " << multDepth - ciphertext_unsort->GetLevel()
            //           << std::endl;
            // std::cout << "error after bitonic_swap is " << err_afterswap << " ~ 2^" << std::log2(err_afterswap)
            //           << std::endl;
            for (int i = 0; i < ciphertext_unsort.size(); i++) {
                if ((multDepth - sort_res_vec[i]->GetLevel()) < 23) {  //23

                    sort_res_vec[i] = cc->EvalMult(sort_res_vec[i], 1 / precision);
                    std::cout << "boot " << std::endl;
                    // int bootprecision = itboot(sort_res, keyPair.secretKey);
                    // std::cout << " sort_res bootprecision: " << bootprecision << std::endl;
                    start_boot = std::chrono::system_clock::now();

                    ciphertext_unsort[i] = cc->EvalBootstrap(sort_res_vec[i], numIterations, 8);
                    // ciphertext_unsort = cc->EvalBootstrap(sort_res);
                    ciphertext_unsort[i] = cc->EvalMult(ciphertext_unsort[i], precision);
                    end_boot             = std::chrono::system_clock::now();
                }
                else {
                    start_boot           = std::chrono::system_clock::now();
                    ciphertext_unsort[i] = sort_res_vec[i];
                    end_boot             = std::chrono::system_clock::now();
                }
            }

            std::cout << "number of levels remaining: " << multDepth - ciphertext_unsort[0]->GetLevel() << std::endl;
            double err_afterBoot = error_estimate_modular(plain_sort_vec[(1 + stage) * stage / 2 + part],
                                                          ciphertext_unsort, keyPair.secretKey, length);

            std::cout << "error afterBoot in stage " << stage << " part " << part << " is " << err_afterBoot << " ~ 2^"
                      << std::log2(err_afterBoot) << std::endl;
            time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
            std::cout << "boot time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count())
                      << "ms" << std::endl;
            std::cout << "swap time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count())
                      << "ms" << std::endl;
            std::cout << "comp time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count())
                      << "ms" << std::endl
                      << std::endl;
        }
        std::vector<double> plain_sort_res(num_slots, 0);
        for (int i = 0; i < ciphertext_unsort.size(); i++) {
            Plaintext plaintextDec;
            cc->Decrypt(keyPair.secretKey, ciphertext_unsort[i], &plaintextDec);
            std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
            for (int j = 0; j < num_slots; j++) {
                plain_sort_res[j] += enc[j].real() * pow(2, 8 * i);
            }
        }
        std::cout << " sortres " << std::endl;
        for (int i = 0; i < pow(2, stage); i++) {
            std::cout << "plain:" << plain_sort_vec[(1 + stage) * stage / 2 + stage][i]
                      << " ,sort: " << plain_sort_res[i] << std::endl;
        }
        std::cout << "stage " << stage << " sort total time: " << time_total << "ms" << std::endl;
        std::cout << "stage " << stage << " sort amortize time: " << time_total / double(length) * (1 << (stage + 1))
                  << "ms" << std::endl
                  << std::endl
                  << std::endl
                  << std::endl;
    }
    end = std::chrono::system_clock::now();

    std::cout << length << "slots sort total time: " << time_total << "ms" << std::endl;
    std::cout << length << "slots sort amortize time: " << time_total / double(length) * num_slots << "ms" << std::endl
              << std::endl;
}

void test_gen_rotate(int plain_bits, int l, int num_slots) {
    std::cout << "\ntest_gen_rotate: " << plain_bits << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits               = 78;
    usint firstMod               = 89;
#else
    usint dcrtBits = 59;
    usint firstMod = 60;
#endif
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetFirstModSize(firstMod);
    std::vector<uint32_t> levelBudget = {4, 4};
    // SecretKeyDist secretKeyDist       = UNIFORM_TERNARY;
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    // Choosing a higher degree yields better precision, but a longer runtime.
    uint32_t polyDegree                    = 495;  //
    uint32_t numIterations                 = 2;
    uint32_t levelsAvailableAfterBootstrap = 10;
    usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    std::cout << "CyclotomicOrder " << cc->GetCyclotomicOrder() << std::endl;
    std::cout << "RingDimension " << cc->GetRingDimension() << std::endl;
    std::cout << "PlaintextModulus " << cc->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(FHE);
    cc->Enable(ADVANCEDSHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.

    double precision = (1 << (plain_bits - 1)) - 1;
    std::cout << "precision: " << precision << std::endl;

    double lowerBound = -precision - 10;
    double upperBound = precision + 10;
    // double bound      = 3;
    auto keyPair  = cc->KeyGen();
    usint ringDim = cc->GetRingDimension();
    int length    = ringDim / 2;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, precision);
    std::vector<std::vector<double>> input_big;
    std::vector<double> test;
    std::vector<Plaintext> plain_big(l);
    std::vector<Ciphertext<lbcrypto::DCRTPoly>> cipher_big(l);

    for (int i = 0; i < l; i++) {
        std::vector<double> input_temp;
        for (int j = 0; j < length; j++) {
            input_temp.push_back(j);
            test.push_back(input_temp[j]);
        }
        input_big.push_back(input_temp);
        plain_big[i]  = cc->MakeCKKSPackedPlaintext(input_temp);
        cipher_big[i] = cc->Encrypt(keyPair.publicKey, plain_big[i]);
    }
    std::vector<std::vector<double>> plain_sort_vec;
    // for (int i = 0; i < plain_sort_vec[0].size(); i++) {
    //     std::cout << plain_sort_vec[0][i] << ",";
    // }
    std::cout << "ringDim: " << ringDim << std::endl;
    cc->EvalMultKeyGen(keyPair.secretKey);
    std::vector<int> rotstep;
    for (int i = 1; i < num_slots; i *= 2) {
        rotstep.push_back(i);
        rotstep.push_back(-i);
    }

    cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
    std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_rotate =
        gen_rotate(length, l, 16, cipher_big, test, keyPair.secretKey);
}

std::vector<std::vector<std::vector<double>>> bitonictopk_test_plain(std::vector<double> vec_unsort,
                                                                     std::vector<std::vector<double>>& plain_sort,
                                                                     uint32_t pack_slots, uint32_t k) {
    // vector<int> vec_unsort = {5, 2, 1, 7, 3, 8, 6, 4};
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<> dis(1, 10000);
    std::vector<std::vector<std::vector<double>>> martix;
    int len = vec_unsort.size();

    for (int stage = 0; int(std::log2(k) + 1) >> stage > 0; stage++) {
        for (int part = 0; stage - part >= 0; part++) {
            // printf("stage: %d part: %d\n", stage, part);
            std::vector<std::vector<double>> swap_martix;  // dimension = 2 or 3, depends on stage and part
            int Block = 1 << (stage + 1);
            int block = 1 << (stage + 1 - part);
            // printf("block: %d\n", block);
            std::vector<double> swap_vector(len, 0);
            std::vector<double> eliminate(len, 0), eliminate_rot(len, 0);
            //1level
            for (int i = 0; i < len; i += block) {
                for (int j = 0; j < block / 2; j++) {
                    eliminate[i + j] = 1;
                }
            }

            for (int i = 0; i < len / block; i++) {
                for (int j = 0; j < block / 2; j++) {
                    int swap_bit;
                    if ((part == 0 && i % 2 == 0) || (part > 0 && (i * block / Block) % 2 == 0))
                        swap_bit = vec_unsort[i * block + j] > vec_unsort[i * block + block / 2 + j] ? 1 : 0;
                    else
                        swap_bit = vec_unsort[i * block + j] < vec_unsort[i * block + block / 2 + j] ? 1 : 0;
                    swap_vector[i * block + j]             = swap_bit;
                    swap_vector[i * block + block / 2 + j] = swap_bit;
                }
            }

            swap_martix.push_back(swap_vector);
            std::vector<double> swap_vector_upper(len, 0);
            std::vector<double> swap_vector_below(len, 0);
            std::vector<double> swap_vector_tmp(len, 0);

            // if (block < len) {
            for (int i = 0; i < len / block; i++)
                for (int j = 0; j < block / 2; j++)
                    swap_vector_tmp[i * block + j] = 1;
            for (int i = 0; i < len; i++)
                swap_vector_below[i] = swap_vector_tmp[i] * (swap_vector_tmp[i] - swap_martix[0][i]);
            rotate_copy(swap_vector_below.begin(), swap_vector_below.begin() + len - block / 2, swap_vector_below.end(),
                        swap_vector_upper.begin());
            swap_martix.push_back(swap_vector_below);
            swap_martix.push_back(swap_vector_upper);
            // }
            // else {
            //     for (int i = 0; i < len; i++)
            //         swap_vector_tmp[i] = 1;
            //     for (int i = 0; i < len; i++)
            //         swap_vector_upper[i] = swap_vector_tmp[i] - swap_martix[0][i];
            //     swap_martix.push_back(swap_vector_upper);
            //     swap_martix.push_back(swap_vector_upper);
            // }
            std::vector<double> vec_sort_temp(len, 0);

            // if (block < len) {
            std::vector<double> vec_rot_upper(len, 0), vec_rot_below(len, 0);
            rotate_copy(vec_unsort.begin(), vec_unsort.begin() + block / 2, vec_unsort.end(), vec_rot_below.begin());
            rotate_copy(vec_unsort.begin(), vec_unsort.begin() + len - block / 2, vec_unsort.end(),
                        vec_rot_upper.begin());

            for (int i = 0; i < len; i++)
                vec_unsort[i] = swap_martix[0][i] * vec_unsort[i] + swap_martix[1][i] * vec_rot_below[i] +
                                swap_martix[2][i] * vec_rot_upper[i];
            // }
            // else {
            //     std::vector<double> vec_rot_below(len, 0);
            //     rotate_copy(vec_unsort.begin(), vec_unsort.begin() + block / 2, vec_unsort.end(),
            //                 vec_rot_below.begin());
            //     for (int i = 0; i < len; i++)
            //         vec_unsort[i] = swap_martix[0][i] * vec_unsort[i] + swap_martix[1][i] * vec_rot_below[i];
            // }
            plain_sort.push_back(vec_unsort);
            martix.push_back(swap_martix);
        }
    }
    std::cout << "partially sort:" << std::endl;
    for (int i = 0; i < len; i++) {
        std::cout << vec_unsort[i] << ",";
    }
    std::cout << std::endl;
    for (int i = 0; i < std::log2(pack_slots) - std::log2(k) - 1; i++) {
        std::cout << "times:" << i << std::endl;
        std::vector<std::vector<double>> swap_martix;  // dimension = 2 or 3, depends on stage and part
        int Block = (1 << int(std::log2(k) + 1)) + (pow(2, i + 1) - 1) * 2 * k;
        int block = 1 << (int(std::log2(k)) + 1);
        int step  = (1 << int(std::log2(k))) + (pow(2, i + 1) - 1) * 2 * k;
        // printf("block: %d\n", block);
        std::vector<double> swap_vector(len, 0);
        std::vector<double> eliminate(len, 0), eliminate_rot(len, 0);
        //1level
        for (int i = 0; i < len; i += block) {
            for (int j = 0; j < block / 2; j++) {
                eliminate[i + j] = 1;
            }
        }

        for (int i = 0; i < len / Block; i++) {
            for (int j = 0; j < block / 2; j++) {
                int swap_bit;
                if (i % 2 == 0)
                    swap_bit = vec_unsort[i * Block + j] < vec_unsort[i * Block + step + j] ? 1 : 0;
                else
                    swap_bit = vec_unsort[i * Block + j] > vec_unsort[i * Block + step + j] ? 1 : 0;
                swap_vector[i * Block + j]        = swap_bit;
                swap_vector[i * Block + step + j] = swap_bit;
                std::cout << "i:" << i * Block + j << ", j" << i * Block + step + j;
                std::cout << ", comp1:" << vec_unsort[i * Block + j] << ", comp2" << vec_unsort[i * Block + step + j]
                          << std::endl;
            }
        }

        swap_martix.push_back(swap_vector);
        std::vector<double> swap_vector_upper(len, 0);
        std::vector<double> swap_vector_below(len, 0);
        std::vector<double> swap_vector_tmp(len, 0);

        // if (block < len) {
        for (int i = 0; i < len / Block; i++)
            for (int j = 0; j < block / 2; j++)
                swap_vector_tmp[i * Block + j] = 1;
        for (int i = 0; i < len; i++)
            swap_vector_below[i] = swap_vector_tmp[i] * (swap_vector_tmp[i] - swap_martix[0][i]);
        rotate_copy(swap_vector_below.begin(), swap_vector_below.begin() + len - step, swap_vector_below.end(),
                    swap_vector_upper.begin());
        swap_martix.push_back(swap_vector_below);
        swap_martix.push_back(swap_vector_upper);
        // }
        // else {
        //     for (int i = 0; i < len; i++)
        //         swap_vector_tmp[i] = 1;
        //     for (int i = 0; i < len; i++)
        //         swap_vector_upper[i] = swap_vector_tmp[i] - swap_martix[0][i];
        //     swap_martix.push_back(swap_vector_upper);
        //     swap_martix.push_back(swap_vector_upper);
        // }
        std::vector<double> vec_sort_temp(len, 0);

        // if (block < len) {
        std::vector<double> vec_rot_upper(len, 0), vec_rot_below(len, 0);
        rotate_copy(vec_unsort.begin(), vec_unsort.begin() + step, vec_unsort.end(), vec_rot_below.begin());
        rotate_copy(vec_unsort.begin(), vec_unsort.begin() + len - step, vec_unsort.end(), vec_rot_upper.begin());

        for (int i = 0; i < len; i++)
            vec_unsort[i] = swap_martix[0][i] * vec_unsort[i] + swap_martix[1][i] * vec_rot_below[i] +
                            swap_martix[2][i] * vec_rot_upper[i];

        plain_sort.push_back(vec_unsort);
        martix.push_back(swap_martix);

        for (int part = 1; part < std::log2(k); part++) {
            // printf("stage: %d part: %d\n", stage, part);
            std::vector<std::vector<double>> swap_martix;  // dimension = 2 or 3, depends on stage and part
            int Block = (1 << int(std::log2(k) + 1)) + (pow(2, i + 1) - 1) * 2 * k;
            int block = 1 << (int(std::log2(k)) + 1);
            int step  = (1 << int(std::log2(k))) + (pow(2, i + 1) - 1) * 2 * k;
            // printf("block: %d\n", block);
            std::vector<double> swap_vector(len, 0);
            std::vector<double> eliminate(len, 0), eliminate_rot(len, 0);
            //1level
            for (int i = 0; i < len; i += block) {
                for (int j = 0; j < block / 2; j++) {
                    eliminate[i + j] = 1;
                }
            }

            for (int i = 0; i < len / Block; i++) {
                for (int j = 0; j < block / 2; j++) {
                    int swap_bit;
                    if ((part == 0 && i % 2 == 0) || (part > 0 && (i * block / Block) % 2 == 0))
                        swap_bit = vec_unsort[i * Block + j] > vec_unsort[i * Block + block / 2 + j] ? 1 : 0;
                    else
                        swap_bit = vec_unsort[i * Block + j] < vec_unsort[i * Block + block / 2 + j] ? 1 : 0;
                    swap_vector[i * Block + j]             = swap_bit;
                    swap_vector[i * Block + block / 2 + j] = swap_bit;
                }
            }

            swap_martix.push_back(swap_vector);
            std::vector<double> swap_vector_upper(len, 0);
            std::vector<double> swap_vector_below(len, 0);
            std::vector<double> swap_vector_tmp(len, 0);

            // if (block < len) {
            for (int i = 0; i < len / Block; i++)
                for (int j = 0; j < block / 2; j++)
                    swap_vector_tmp[i * Block + j] = 1;
            for (int i = 0; i < len; i++)
                swap_vector_below[i] = swap_vector_tmp[i] * (swap_vector_tmp[i] - swap_martix[0][i]);
            rotate_copy(swap_vector_below.begin(), swap_vector_below.begin() + len - step, swap_vector_below.end(),
                        swap_vector_upper.begin());
            swap_martix.push_back(swap_vector_below);
            swap_martix.push_back(swap_vector_upper);
            // }
            std::vector<double> vec_sort_temp(len, 0);

            // if (block < len) {
            std::vector<double> vec_rot_upper(len, 0), vec_rot_below(len, 0);
            rotate_copy(vec_unsort.begin(), vec_unsort.begin() + step, vec_unsort.end(), vec_rot_below.begin());
            rotate_copy(vec_unsort.begin(), vec_unsort.begin() + len - step, vec_unsort.end(), vec_rot_upper.begin());

            for (int i = 0; i < len; i++)
                vec_unsort[i] = swap_martix[0][i] * vec_unsort[i] + swap_martix[1][i] * vec_rot_below[i] +
                                swap_martix[2][i] * vec_rot_upper[i];

            plain_sort.push_back(vec_unsort);
            martix.push_back(swap_martix);
        }
        std::cout << "after topk" << std::endl;
        for (int i = 0; i < len; i++) {
            std::cout << vec_unsort[i] << ",";
        }
        std::cout << std::endl;
    }

    return martix;
}

void topk_test(int plain_bits, int num_slots, int k) {
    std::cout << "Top-" << k << " test " << std::endl;

    // We need to enable Advanced SHE to use the Chebyshev approximation.

    double precision = (1 << (plain_bits - 1)) - 1;
    std::cout << "precision: " << precision << std::endl;

    double lowerBound = -precision - 10;
    double upperBound = precision + 10;
    // double bound      = 3;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, precision);
    std::vector<std::vector<double>> column;
    std::vector<double> input(num_slots, 0);
    for (int i = 0; i < num_slots; i++) {
        input[i] = message(engine);
    }
    input[0]  = 3;
    input[1]  = 2;
    input[2]  = 4;
    input[3]  = 1;
    input[4]  = 5;
    input[5]  = 7;
    input[6]  = 6;
    input[7]  = 8;
    input[8]  = 9;
    input[9]  = 11;
    input[10] = 10;
    input[11] = 12;
    input[12] = 16;
    input[13] = 14;
    input[14] = 15;
    input[15] = 13;
    std::vector<std::vector<double>> plain_sort_vec;
    std::vector<std::vector<std::vector<double>>> matrix = bitonictopk_test_plain(input, plain_sort_vec, num_slots, k);
}
//4096

std::vector<std::vector<std::vector<double>>> bitonictopk_test_plain_1(std::vector<double> vec_unsort,
                                                                       std::vector<std::vector<double>>& plain_sort,
                                                                       uint32_t pack_slots, uint32_t k) {
    // vector<int> vec_unsort = {5, 2, 1, 7, 3, 8, 6, 4};
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<> dis(1, 10000);
    std::vector<std::vector<std::vector<double>>> martix;
    int len = vec_unsort.size();

    for (int stage = 0; int(std::log2(k) + 1) >> stage > 0; stage++) {
        for (int part = 0; stage - part >= 0; part++) {
            // printf("stage: %d part: %d\n", stage, part);
            std::vector<std::vector<double>> swap_martix;  // dimension = 2 or 3, depends on stage and part
            int Block = 1 << (stage + 1);
            int block = 1 << (stage + 1 - part);
            // printf("block: %d\n", block);
            std::vector<double> swap_vector(len, 0);
            std::vector<double> eliminate(len, 0), eliminate_rot(len, 0);
            //1level
            for (int i = 0; i < len; i += block) {
                for (int j = 0; j < block / 2; j++) {
                    eliminate[i + j] = 1;
                }
            }

            for (int i = 0; i < len / block; i++) {
                for (int j = 0; j < block / 2; j++) {
                    int swap_bit;
                    if ((part == 0 && i % 2 == 0) || (part > 0 && (i * block / Block) % 2 == 0))
                        swap_bit = vec_unsort[i * block + j] > vec_unsort[i * block + block / 2 + j] ? 1 : 0;
                    else
                        swap_bit = vec_unsort[i * block + j] < vec_unsort[i * block + block / 2 + j] ? 1 : 0;
                    swap_vector[i * block + j]             = swap_bit;
                    swap_vector[i * block + block / 2 + j] = swap_bit;
                }
            }

            swap_martix.push_back(swap_vector);
            std::vector<double> swap_vector_upper(len, 0);
            std::vector<double> swap_vector_below(len, 0);
            std::vector<double> swap_vector_tmp(len, 0);

            // if (block < len) {
            for (int i = 0; i < len / block; i++)
                for (int j = 0; j < block / 2; j++)
                    swap_vector_tmp[i * block + j] = 1;
            for (int i = 0; i < len; i++)
                swap_vector_below[i] = swap_vector_tmp[i] * (swap_vector_tmp[i] - swap_martix[0][i]);
            rotate_copy(swap_vector_below.begin(), swap_vector_below.begin() + len - block / 2, swap_vector_below.end(),
                        swap_vector_upper.begin());
            swap_martix.push_back(swap_vector_below);
            swap_martix.push_back(swap_vector_upper);
            // }
            // else {
            //     for (int i = 0; i < len; i++)
            //         swap_vector_tmp[i] = 1;
            //     for (int i = 0; i < len; i++)
            //         swap_vector_upper[i] = swap_vector_tmp[i] - swap_martix[0][i];
            //     swap_martix.push_back(swap_vector_upper);
            //     swap_martix.push_back(swap_vector_upper);
            // }
            std::vector<double> vec_sort_temp(len, 0);

            // if (block < len) {
            std::vector<double> vec_rot_upper(len, 0), vec_rot_below(len, 0);
            rotate_copy(vec_unsort.begin(), vec_unsort.begin() + block / 2, vec_unsort.end(), vec_rot_below.begin());
            rotate_copy(vec_unsort.begin(), vec_unsort.begin() + len - block / 2, vec_unsort.end(),
                        vec_rot_upper.begin());

            for (int i = 0; i < len; i++)
                vec_unsort[i] = swap_martix[0][i] * vec_unsort[i] + swap_martix[1][i] * vec_rot_below[i] +
                                swap_martix[2][i] * vec_rot_upper[i];
            // }
            // else {
            //     std::vector<double> vec_rot_below(len, 0);
            //     rotate_copy(vec_unsort.begin(), vec_unsort.begin() + block / 2, vec_unsort.end(),
            //                 vec_rot_below.begin());
            //     for (int i = 0; i < len; i++)
            //         vec_unsort[i] = swap_martix[0][i] * vec_unsort[i] + swap_martix[1][i] * vec_rot_below[i];
            // }
            plain_sort.push_back(vec_unsort);
            martix.push_back(swap_martix);
        }
    }
    std::cout << "partially sort:" << std::endl;
    for (int i = 0; i < len; i++) {
        std::cout << vec_unsort[i] << ",";
    }
    std::cout << std::endl;
    for (int i = 0; i < std::log2(pack_slots) - std::log2(k) - 1; i++) {
        std::cout << "times:" << i << std::endl;
        std::vector<std::vector<double>> swap_martix;  // dimension = 2 or 3, depends on stage and part
        int Block = (1 << int(std::log2(k) + 1)) + (pow(2, i + 1) - 1) * 2 * k;
        int block = 1 << (int(std::log2(k)) + 1);
        int step  = (1 << int(std::log2(k))) + (pow(2, i + 1) - 1) * 2 * k;
        // printf("block: %d\n", block);
        std::vector<double> swap_vector(len, 0);
        std::vector<double> eliminate(len, 0), eliminate_rot(len, 0);
        //1level
        for (int i = 0; i < len; i += block) {
            for (int j = 0; j < block / 2; j++) {
                eliminate[i + j] = 1;
            }
        }

        for (int i = 0; i < len / Block; i++) {
            for (int j = 0; j < block / 2; j++) {
                int swap_bit;
                if (i % 2 == 0)
                    swap_bit = vec_unsort[i * Block + j] > vec_unsort[i * Block + step + j] ? 1 : 0;
                else
                    swap_bit = vec_unsort[i * Block + j] < vec_unsort[i * Block + step + j] ? 1 : 0;
                swap_vector[i * Block + j]        = swap_bit;
                swap_vector[i * Block + step + j] = swap_bit;
                std::cout << "i:" << i * Block + j << ", j" << i * Block + step + j;
                std::cout << ", comp1:" << vec_unsort[i * Block + j] << ", comp2" << vec_unsort[i * Block + step + j]
                          << std::endl;
            }
        }

        swap_martix.push_back(swap_vector);
        std::vector<double> swap_vector_upper(len, 0);
        std::vector<double> swap_vector_below(len, 0);
        std::vector<double> swap_vector_tmp(len, 0);

        // if (block < len) {
        for (int i = 0; i < len / Block; i++)
            for (int j = 0; j < block / 2; j++)
                swap_vector_tmp[i * Block + j] = 1;
        for (int i = 0; i < len; i++)
            swap_vector_below[i] = swap_vector_tmp[i] * (swap_vector_tmp[i] - swap_martix[0][i]);
        rotate_copy(swap_vector_below.begin(), swap_vector_below.begin() + len - step, swap_vector_below.end(),
                    swap_vector_upper.begin());
        swap_martix.push_back(swap_vector_below);
        swap_martix.push_back(swap_vector_upper);
        // }
        // else {
        //     for (int i = 0; i < len; i++)
        //         swap_vector_tmp[i] = 1;
        //     for (int i = 0; i < len; i++)
        //         swap_vector_upper[i] = swap_vector_tmp[i] - swap_martix[0][i];
        //     swap_martix.push_back(swap_vector_upper);
        //     swap_martix.push_back(swap_vector_upper);
        // }
        std::vector<double> vec_sort_temp(len, 0);

        // if (block < len) {
        std::vector<double> vec_rot_upper(len, 0), vec_rot_below(len, 0);
        rotate_copy(vec_unsort.begin(), vec_unsort.begin() + step, vec_unsort.end(), vec_rot_below.begin());
        rotate_copy(vec_unsort.begin(), vec_unsort.begin() + len - step, vec_unsort.end(), vec_rot_upper.begin());

        for (int i = 0; i < len; i++)
            vec_unsort[i] = swap_martix[0][i] * vec_unsort[i] + swap_martix[1][i] * vec_rot_below[i] +
                            swap_martix[2][i] * vec_rot_upper[i];

        plain_sort.push_back(vec_unsort);
        martix.push_back(swap_martix);

        for (int part = 1; part < std::log2(k); part++) {
            // printf("stage: %d part: %d\n", stage, part);
            std::vector<std::vector<double>> swap_martix;  // dimension = 2 or 3, depends on stage and part
            int Block = 1 << int(std::log2(k) + 1);
            int block = 1 << int(std::log2(k) + 1 - part);
            // printf("block: %d\n", block);
            std::vector<double> swap_vector(len, 0);
            std::vector<double> eliminate(len, 0), eliminate_rot(len, 0);
            //1level
            for (int i = 0; i < len; i += block) {
                for (int j = 0; j < block / 2; j++) {
                    eliminate[i + j] = 1;
                }
            }

            for (int i = 0; i < len / Block; i++) {
                for (int j = 0; j < block / 2; j++) {
                    int swap_bit;
                    if ((part == 0 && i % 2 == 0) || (part > 0 && (i * block / Block) % 2 == 0))
                        swap_bit = vec_unsort[i * Block + j] > vec_unsort[i * Block + block / 2 + j] ? 1 : 0;
                    else
                        swap_bit = vec_unsort[i * Block + j] < vec_unsort[i * Block + block / 2 + j] ? 1 : 0;
                    swap_vector[i * Block + j]             = swap_bit;
                    swap_vector[i * Block + block / 2 + j] = swap_bit;
                }
            }

            swap_martix.push_back(swap_vector);
            std::vector<double> swap_vector_upper(len, 0);
            std::vector<double> swap_vector_below(len, 0);
            std::vector<double> swap_vector_tmp(len, 0);

            for (int i = 0; i < len / block; i++)
                for (int j = 0; j < block / 2; j++)
                    swap_vector_tmp[i * block + j] = 1;
            for (int i = 0; i < len; i++)
                swap_vector_below[i] = swap_vector_tmp[i] * (swap_vector_tmp[i] - swap_martix[0][i]);
            rotate_copy(swap_vector_below.begin(), swap_vector_below.begin() + len - block / 2, swap_vector_below.end(),
                        swap_vector_upper.begin());
            swap_martix.push_back(swap_vector_below);
            swap_martix.push_back(swap_vector_upper);
            // }
            std::vector<double> vec_sort_temp(len, 0);

            // if (block < len) {
            std::vector<double> vec_rot_upper(len, 0), vec_rot_below(len, 0);
            rotate_copy(vec_unsort.begin(), vec_unsort.begin() + block / 2, vec_unsort.end(), vec_rot_below.begin());
            rotate_copy(vec_unsort.begin(), vec_unsort.begin() + len - block / 2, vec_unsort.end(),
                        vec_rot_upper.begin());

            for (int i = 0; i < len; i++)
                vec_unsort[i] = swap_martix[0][i] * vec_unsort[i] + swap_martix[1][i] * vec_rot_below[i] +
                                swap_martix[2][i] * vec_rot_upper[i];
            plain_sort.push_back(vec_unsort);
            martix.push_back(swap_martix);
        }
        std::cout << "after topk" << std::endl;
        for (int i = 0; i < len; i++) {
            std::cout << vec_unsort[i] << ",";
        }
        std::cout << std::endl;
    }

    return martix;
}

void bitonic_comp_topk(int stage, int slots, int topk_i, int k, Ciphertext<lbcrypto::DCRTPoly>& ct,
                       Ciphertext<lbcrypto::DCRTPoly>& res, double precision, std::vector<double>& coefficients,
                       lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey) {
    auto cc   = ct->GetCryptoContext();
    int Block = (1 << int(std::log2(k) + 1)) + (pow(2, topk_i + 1) - 1) * 2 * k;
    int block = 1 << (int(std::log2(k)) + 1);
    int step  = (1 << int(std::log2(k))) + (pow(2, topk_i + 1) - 1) * 2 * k;

    std::vector<double> mask(slots, 0), eliminate(slots, 0);
    // if (block < slots)
    for (int i = 0; i < slots; i += Block) {
        for (int j = 0; j < block / 2; j++) {
            mask[i + j] = 1;
        }
        // std::cout << "mask i:" << i << "  block / 2: " << block / 2 << std::endl;
        // std::cout << "mask i:" << i << " i + 2 * block - 1: " << i + 2 * block - 1 << std::endl;
        for (int j = step; j < block / 2; j++) {
            mask[i + j] = 1;
        }
    }
    // else

    for (int i = 0; i < slots; i += Block) {
        for (int j = 0; j < block / 2; j++) {
            eliminate[i + j] = 1;
        }
        for (int j = step; j < block / 2; j++) {
            eliminate[i + j] = 1;
        }
        // std::cout << "eliminate i:" << i << "  i + 2 * block - 1: " << i + 2 * block - 1 << std::endl;
    }
    // for (int i = 0; i < 8; i++) {
    //     std::cout << "i:" << i << " mask: " << mask[i] << ", eliminate: " << eliminate[i] << std::endl;
    // }

    Plaintext plain_mask      = cc->MakeCKKSPackedPlaintext(mask);
    Plaintext plain_eliminate = cc->MakeCKKSPackedPlaintext(eliminate);
    auto ct_mask              = cc->EvalMult(ct, plain_mask);
    // std::cout << "\n EvalRotate: " << block / 2 << std::endl;
    auto ct_rot = cc->EvalRotate(ct_mask, step);

    // std::cout << "\n comp_partial " << std::endl;
    ct_mask = cc->EvalMult(ct_mask, plain_eliminate);

    ct_rot = cc->EvalMult(ct_rot, plain_eliminate);

    // std::cout << " before sign  " << std::endl;
    // output(ct_rot, 8, privateKey);
    // std::cout << std::endl;
    Ciphertext<lbcrypto::DCRTPoly> ploy_res;
    // comp_greater_than(ct_mask, ct_rot, precision, coefficients, ploy_res, privateKey);
    // comp_partial(ct_mask, ct_rot, precision, coefficients, ploy_res, privateKey);
    // std::cout << "\n Eval comp_partial: " << block / 2 << std::endl;
    comp_partial(ct_mask, ct_rot, precision, coefficients, res, privateKey);
    res = cc->EvalMult(res, plain_eliminate);

    // std::cout << "compres  " << std::endl;
    // output(res, 8, privateKey);
    // std::cout << std::endl;

    // auto ct_rot = cc->EvalRotate(res, -Block / 2);
    // res         = cc->EvalAdd(res, ct_rot);
}

void bitonic_swap_topk(int slots, int topk_i, int k, Ciphertext<lbcrypto::DCRTPoly>& compres,
                       Ciphertext<lbcrypto::DCRTPoly>& ct, Ciphertext<lbcrypto::DCRTPoly>& res,
                       lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey) {
    std::chrono::system_clock::time_point start, end;
    double time = 0.;
    start       = std::chrono::system_clock::now();
    auto cc     = ct->GetCryptoContext();
    // int Block = 1 << (stage + 1);
    int Block = (1 << int(std::log2(k) + 1)) + (pow(2, topk_i + 1) - 1) * 2 * k;
    int block = 1 << (int(std::log2(k)) + 1);
    int step  = (1 << int(std::log2(k))) + (pow(2, topk_i + 1) - 1) * 2 * k;
    std::vector<Ciphertext<lbcrypto::DCRTPoly>> res_martrix;
    std::vector<double> eliminate(slots, 0), eliminate_rot(slots, 0);
    //1level
    for (int i = 0; i < slots; i += Block) {
        for (int j = 0; j < block / 2; j++) {
            eliminate[i + j] = 1;
        }
        for (int j = step; j < block / 2; j++) {
            eliminate[i + j] = 1;
        }
    }
    std::rotate_copy(eliminate.begin(), eliminate.begin() + block / 2, eliminate.end(), eliminate_rot.begin());

    Plaintext plaintextDec;
    Plaintext plain_eliminate = cc->MakeCKKSPackedPlaintext(eliminate);
    end                       = std::chrono::system_clock::now();
    time                      = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    // std::cout << "preprocess time: " << time << "ms" << std::endl;
    // Plaintext plain_eliminate_rot = cc->MakeCKKSPackedPlaintext(eliminate_rot);
    // if (block < slots) {
    // std::cout << " EvalRotate: " << std::endl;
    auto ct_upper = cc->EvalRotate(ct, -step);
    auto ct_lower = cc->EvalRotate(ct, step);

    // ct_upper = cc->EvalMult(ct_upper, plain_eliminate_rot);
    // ct_lower = cc->EvalMult(ct_lower, plain_eliminate);

    std::vector<std::complex<double>> enc;
    start            = std::chrono::system_clock::now();
    auto compres_neg = cc->EvalNegate(compres);
    auto compres_rot = cc->EvalRotate(compres, -step);

    auto swap_vec_lower = cc->EvalAdd(compres_neg, plain_eliminate);

    auto swap_vec_upper = cc->EvalRotate(swap_vec_lower, -step);

    auto swap_vec_diagonal = cc->EvalAdd(compres, compres_rot);
    end                    = std::chrono::system_clock::now();
    time                   = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    // std::cout << " EvalMult: " << std::endl;
    start = std::chrono::system_clock::now();
    res_martrix.push_back(cc->EvalMult(ct_upper, swap_vec_upper));
    res_martrix.push_back(cc->EvalMult(ct, swap_vec_diagonal));
    res_martrix.push_back(cc->EvalMult(ct_lower, swap_vec_lower));
    // }
    res = cc->EvalAdd(res_martrix[0], res_martrix[1]);
    res = cc->EvalAdd(res, res_martrix[2]);
    // res = cc->EvalMult(res, plain_eliminate);
    // cc->Decrypt(privateKey, res, &plaintextDec);
    // plaintextDec->SetLength(slots);
    // enc = plaintextDec->GetCKKSPackedValue();
    // for (int i = 0; i < 16; i++) {
    //     std::cout << " swap_res: " << enc[i].real() << std::endl;
    // }
    end  = std::chrono::system_clock::now();
    time = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
}

void topk_sort(int plain_bits, int num_slots, int k) {
    std::cout << "top " << k << " in " << num_slots << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits               = 78;
    usint firstMod               = 89;
#else
    usint dcrtBits = 59;
    usint firstMod = 60;
#endif
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetFirstModSize(firstMod);
    std::vector<uint32_t> levelBudget = {4, 4};
    // SecretKeyDist secretKeyDist       = UNIFORM_TERNARY;
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    // Choosing a higher degree yields better precision, but a longer runtime.
    uint32_t polyDegree = 495;  //

    // The multiplicative depth depends on the polynomial degree.
    // See the FUNCTION_EVALUATION.md file for a table mapping polynomial degrees to multiplicative depths.
    // uint32_t levelsAvailableAfterBootstrap = 35;  //35
    // uint32_t levelsbeforeBootstrap         = 15;
    uint32_t numIterations = 2;
    // usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist) +
    //                   (numIterations - 1) + levelsbeforeBootstrap;
    uint32_t levelsAvailableAfterBootstrap = 25;
    usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    std::cout << "polyDegree " << polyDegree << std::endl;
    std::cout << "CyclotomicOrder " << cc->GetCyclotomicOrder() << std::endl;
    std::cout << "RingDimension " << cc->GetRingDimension() << std::endl;
    // std::cout << "Modulus " << cc->GetModulus() << std::endl;
    std::cout << "PlaintextModulus " << cc->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(FHE);
    cc->Enable(ADVANCEDSHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.

    double precision = (1 << (plain_bits - 1)) - 1;
    std::cout << "precision: " << precision << std::endl;

    double lowerBound = -precision - 10;
    double upperBound = precision + 10;
    // double bound      = 3;
    auto keyPair  = cc->KeyGen();
    usint ringDim = cc->GetRingDimension();
    int length    = ringDim / 2;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, precision);
    std::vector<double> input(length, 0);
    const std::vector<DCRTPoly>& ckks_pk = keyPair.publicKey->GetPublicElements();
    // std::cout << "Moduli chain of pk: " << std::endl;
    // print_moduli_chain(ckks_pk[0]);
    for (int i = 0; i < length; i++) {
        input[i] = message(engine);
    }
    input[0]  = 3;
    input[1]  = 1;
    input[2]  = 4;
    input[3]  = 1;
    input[4]  = 5;
    input[5]  = 7;
    input[6]  = 6;
    input[7]  = 8;
    input[8]  = 9;
    input[9]  = 11;
    input[10] = 10;
    input[11] = 12;
    input[12] = 16;
    input[13] = 14;
    input[14] = 15;
    input[15] = 13;
    // for (int i = 0; i < 8; i++) {
    //     std::cout << input[i] << std::endl;
    // }
    std::vector<std::vector<double>> plain_sort_vec;
    std::vector<std::vector<std::vector<double>>> matrix = bitonic_test_plain(input, plain_sort_vec, num_slots);
    // for (int i = 0; i < plain_sort_vec[0].size(); i++) {
    //     std::cout << plain_sort_vec[0][i] << ",";
    // }
    std::cout << "ringDim: " << ringDim << std::endl;
    cc->EvalMultKeyGen(keyPair.secretKey);
    std::vector<int> rotstep;
    for (int i = 1; i < num_slots; i *= 2) {
        rotstep.push_back(i);
        rotstep.push_back(-i);
    }
    for (int i = 0; i < std::log2(num_slots) - std::log2(k) - 1; i++) {
        rotstep.push_back((1 << int(std::log2(k))) + (pow(2, i + 1) - 1) * 2 * k);
        rotstep.push_back(-((1 << int(std::log2(k))) + (pow(2, i + 1) - 1) * 2 * k));
    }
    cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
    cc->EvalBootstrapSetup(levelBudget);
    cc->EvalBootstrapKeyGen(keyPair.secretKey, length);

    std::vector<double> coefficients = EvalChebyshevCoefficients([](double x) -> double { return ((x > 0) ? 1 : 0); },
                                                                 lowerBound, upperBound, polyDegree);
    Plaintext plain                  = cc->MakeCKKSPackedPlaintext(input);
    auto ciphertext_unsort           = cc->Encrypt(keyPair.publicKey, plain);

    std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
    std::cout << "number of levels fresh: " << multDepth - ciphertext_unsort->GetLevel() << std::endl;
    start             = std::chrono::system_clock::now();
    double time_total = 0;
    for (int stage = 0; int(std::log2(k) + 1) >> stage > 0; stage++) {
        for (int part = 0; stage - part >= 0; part++) {
            std::cout << " stage " << stage << " part " << part << std::endl;
            Ciphertext<lbcrypto::DCRTPoly> comp_res;
            Ciphertext<lbcrypto::DCRTPoly> sort_res;
            std::cout << "bitonic_comp " << std::endl;
            start_comp = std::chrono::system_clock::now();
            bitonic_comp(stage, part, length, ciphertext_unsort, comp_res, precision, coefficients, keyPair.secretKey);
            end_comp = std::chrono::system_clock::now();
            // comp_res = cc->EvalBootstrap(comp_res);
            // Plaintext plaintextDec;
            // cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec);
            // std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
            // for (int i = 0; i < 2 * int(num_slots); i++) {
            //     std::cout << " comp_res: " << enc[i].real() << std::endl;
            // }
            // std::cout << "number of levels remaining before bitonic_swap: " << multDepth - ciphertext_unsort->GetLevel()
            //           << std::endl;
            std::cout << "bitonic_swap " << std::endl;
            start_swap = std::chrono::system_clock::now();
            bitonic_swap(stage, part, length, comp_res, ciphertext_unsort, sort_res, keyPair.secretKey);
            end_swap = std::chrono::system_clock::now();
            // std::cout << "bitonic_swap " << std::endl;
            double err_afterswap = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], sort_res,
                                                  keyPair.secretKey, 2 * num_slots);
            // std::cout << "number of levels remaining after bitonic_swap: " << multDepth - ciphertext_unsort->GetLevel()
            //           << std::endl;
            std::cout << "error after bitonic_swap is " << err_afterswap << " ~ 2^" << std::log2(err_afterswap)
                      << std::endl;

            if ((multDepth - sort_res->GetLevel()) < 23) {  //23
                sort_res = cc->EvalMult(sort_res, 1 / precision);
                std::cout << "boot " << std::endl;
                // int bootprecision = itboot(sort_res, keyPair.secretKey);
                // std::cout << " sort_res bootprecision: " << bootprecision << std::endl;
                start_boot        = std::chrono::system_clock::now();
                ciphertext_unsort = cc->EvalBootstrap(sort_res, numIterations, 8);
                // ciphertext_unsort = cc->EvalBootstrap(sort_res);
                ciphertext_unsort = cc->EvalMult(ciphertext_unsort, precision);
                end_boot          = std::chrono::system_clock::now();
                std::cout << "number of levels remaining  after boot: " << multDepth - ciphertext_unsort->GetLevel()
                          << std::endl;
            }
            else {
                start_boot        = std::chrono::system_clock::now();
                ciphertext_unsort = sort_res;
                end_boot          = std::chrono::system_clock::now();
                std::cout << "number of levels remaining : " << multDepth - ciphertext_unsort->GetLevel() << std::endl;
            }

            double err_afterBoot = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], ciphertext_unsort,
                                                  keyPair.secretKey, 2 * num_slots);

            std::cout << "error afterBoot in stage " << stage << " part " << part << " is " << err_afterBoot << " ~ 2^"
                      << std::log2(err_afterBoot) << std::endl;
            time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
            std::cout << "boot time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count())
                      << "ms" << std::endl;
            std::cout << "swap time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count())
                      << "ms" << std::endl;
            std::cout << "comp time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count())
                      << "ms" << std::endl
                      << std::endl;
        }

        std::cout << "stage " << stage << " sort total time: " << time_total << "ms" << std::endl;
        std::cout << "stage " << stage << " sort amortize time: " << time_total / double(length) * (1 << (stage + 1))
                  << "ms" << std::endl
                  << std::endl
                  << std::endl
                  << std::endl;
    }
    Plaintext plaintextDec_temp;
    for (int i = 0; i < std::log2(num_slots) - std::log2(k) - 1; i++) {
        std::cout << " itration " << i << std::endl;

        cc->Decrypt(keyPair.secretKey, ciphertext_unsort, &plaintextDec_temp);
        std::vector<std::complex<double>> enc = plaintextDec_temp->GetCKKSPackedValue();
        std::cout << " sortres1 " << std::endl;
        for (int i = 0; i < int(num_slots); i++) {
            std::cout << enc[i].real() << ", ";
        }
        std::cout << std::endl;

        int Block = (1 << int(std::log2(k) + 1)) + (pow(2, i + 1) - 1) * 2 * k;
        // std::cout << " Block: " << Block << std::endl;
        int block = 1 << (int(std::log2(k)) + 1);
        // std::cout << " block: " << block << std::endl;
        int step = (1 << int(std::log2(k))) + (pow(2, i + 1) - 1) * 2 * k;
        // std::cout << " step: " << step << std::endl;
        Ciphertext<lbcrypto::DCRTPoly> comp_res;
        Ciphertext<lbcrypto::DCRTPoly> sort_res;
        std::cout << "bitonic_comp_topk " << std::endl;
        bitonic_comp_topk(std::log2(k), length, i, k, ciphertext_unsort, comp_res, precision, coefficients,
                          keyPair.secretKey);
        cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec_temp);
        enc = plaintextDec_temp->GetCKKSPackedValue();
        std::cout << " comp res " << std::endl;
        for (int i = 0; i < int(num_slots); i++) {
            std::cout << enc[i].real() << ", ";
        }
        std::cout << std::endl;
        std::cout << "bitonic_swap_topk " << std::endl;
        bitonic_swap_topk(length, i, k, comp_res, ciphertext_unsort, sort_res, keyPair.secretKey);

        cc->Decrypt(keyPair.secretKey, sort_res, &plaintextDec_temp);
        enc = plaintextDec_temp->GetCKKSPackedValue();
        std::cout << " sortres2 " << std::endl;
        for (int i = 0; i < int(num_slots); i++) {
            std::cout << enc[i].real() << ", ";
        }
        std::cout << std::endl;
        if ((multDepth - sort_res->GetLevel()) < 23) {  //23
            sort_res = cc->EvalMult(sort_res, 1 / precision);
            std::cout << "boot " << std::endl;
            // int bootprecision = itboot(sort_res, keyPair.secretKey);
            // std::cout << " sort_res bootprecision: " << bootprecision << std::endl;
            start_boot        = std::chrono::system_clock::now();
            ciphertext_unsort = cc->EvalBootstrap(sort_res, numIterations, 8);
            // ciphertext_unsort = cc->EvalBootstrap(sort_res);
            ciphertext_unsort = cc->EvalMult(ciphertext_unsort, precision);
            end_boot          = std::chrono::system_clock::now();
            std::cout << "number of levels remaining  after boot: " << multDepth - ciphertext_unsort->GetLevel()
                      << std::endl;
        }
        else {
            start_boot        = std::chrono::system_clock::now();
            ciphertext_unsort = sort_res;
            end_boot          = std::chrono::system_clock::now();
            std::cout << "number of levels remaining : " << multDepth - ciphertext_unsort->GetLevel() << std::endl;
        }
        cc->Decrypt(keyPair.secretKey, ciphertext_unsort, &plaintextDec_temp);
        enc = plaintextDec_temp->GetCKKSPackedValue();
        std::cout << " sortres3 " << std::endl;
        for (int i = 0; i < int(num_slots); i++) {
            std::cout << enc[i].real() << ", ";
        }
        std::cout << std::endl;
        // double err_afterBoot = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], ciphertext_unsort,
        //                                       keyPair.secretKey, 2 * num_slots);

        for (int part = 1; part < std::log2(k); part++) {
            int stage = std::log2(k);
            bitonic_comp(stage, part, length, ciphertext_unsort, comp_res, precision, coefficients, keyPair.secretKey);
            bitonic_swap(stage, part, length, comp_res, ciphertext_unsort, sort_res, keyPair.secretKey);
            end_swap = std::chrono::system_clock::now();
            // double err_afterswap = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], sort_res,
            //                                       keyPair.secretKey, 2 * num_slots);
            // std::cout << "number of levels remaining after bitonic_swap: " << multDepth - ciphertext_unsort->GetLevel()
            //           << std::endl;
            // std::cout << "error after bitonic_swap is " << err_afterswap << " ~ 2^" << std::log2(err_afterswap)
            //           << std::endl;

            if ((multDepth - sort_res->GetLevel()) < 23) {  //23
                sort_res = cc->EvalMult(sort_res, 1 / precision);
                std::cout << "boot " << std::endl;
                // int bootprecision = itboot(sort_res, keyPair.secretKey);
                // std::cout << " sort_res bootprecision: " << bootprecision << std::endl;
                start_boot        = std::chrono::system_clock::now();
                ciphertext_unsort = cc->EvalBootstrap(sort_res, numIterations, 8);
                // ciphertext_unsort = cc->EvalBootstrap(sort_res);
                ciphertext_unsort = cc->EvalMult(ciphertext_unsort, precision);
                end_boot          = std::chrono::system_clock::now();
                std::cout << "number of levels remaining  after boot: " << multDepth - ciphertext_unsort->GetLevel()
                          << std::endl;
            }
            else {
                start_boot        = std::chrono::system_clock::now();
                ciphertext_unsort = sort_res;
                end_boot          = std::chrono::system_clock::now();
                std::cout << "number of levels remaining : " << multDepth - ciphertext_unsort->GetLevel() << std::endl;
            }

            // double err_afterBoot = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], ciphertext_unsort,
            //                                       keyPair.secretKey, 2 * num_slots);

            // std::cout << "error afterBoot in stage " << stage << " part " << part << " is " << err_afterBoot << " ~ 2^"
            //           << std::log2(err_afterBoot) << std::endl;
            time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
            std::cout << "boot time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count())
                      << "ms" << std::endl;
            std::cout << "swap time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count())
                      << "ms" << std::endl;
            std::cout << "comp time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count())
                      << "ms" << std::endl
                      << std::endl;
        }
    }

    end = std::chrono::system_clock::now();
    Plaintext plaintextDec;
    cc->Decrypt(keyPair.secretKey, ciphertext_unsort, &plaintextDec);
    std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
    std::cout << " sortres " << std::endl;
    for (int i = 0; i < int(num_slots); i++) {
        std::cout << enc[i].real() << std::endl;
    }
    std::cout << length << "slots sort total time: " << time_total << "ms" << std::endl;
    std::cout << length << "slots sort amortize time: " << time_total / double(length) * num_slots << "ms" << std::endl
              << std::endl;
}

int main(int argc, char* argv[]) {
    // topk_sort(8, 16, 4);
    // bitonic_sort_small(8, 8);
    // test_gen_rotate(8, 2, 65536);
    // bitonic_sort_modular(8, 2, 4);
    bitonic_sort_modular(8, 2, 65536);
    // topk_test(16, 256, 8);
    // cross_rotate(3, 4, 7);
    // for (int i = 4; i < 65537; i *= 2)
    // sync_test_small(8, i);
    // for (int i = 1; i < 17; i *= 2)
    //     sync_test_big(i, 8, 65536);
    return 0;
}