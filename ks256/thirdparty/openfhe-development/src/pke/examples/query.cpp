#include "openfhe.h"
#include <chrono>
#include <omp.h>
#include "math/chebyshev.h"
using namespace lbcrypto;
template <typename T>
inline uint64_t CeilLog2(T x) {
    return static_cast<uint64_t>(std::ceil(std::log2(x)));
}
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

//ct1-ct2 >0 -> res=1; ct1-ct2 ==0 -> res=0.5; ct1-ct2 <0 -> res=0
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

    ct_mask = cc->EvalMult(ct_mask, plain_eliminate);

    ct_rot = cc->EvalMult(ct_rot, plain_eliminate);

    Ciphertext<lbcrypto::DCRTPoly> ploy_res;

    comp_partial(ct_mask, ct_rot, precision, coefficients, ploy_res, privateKey);
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
    std::cout << "preprocess time: " << time << "ms" << std::endl;

    auto ct_upper = cc->EvalRotate(ct, -block / 2);
    auto ct_lower = cc->EvalRotate(ct, block / 2);

    std::vector<std::complex<double>> enc;
    start            = std::chrono::system_clock::now();
    auto compres_neg = cc->EvalNegate(compres);
    auto compres_rot = cc->EvalRotate(compres, -block / 2);

    auto swap_vec_lower = cc->EvalAdd(compres_neg, plain_eliminate);

    auto swap_vec_upper = cc->EvalRotate(swap_vec_lower, -block / 2);

    auto swap_vec_diagonal = cc->EvalAdd(compres, compres_rot);

    end  = std::chrono::system_clock::now();
    time = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "rotate time: " << time << "ms" << std::endl;

    start = std::chrono::system_clock::now();
    res_martrix.push_back(cc->EvalMult(ct_upper, swap_vec_upper));
    res_martrix.push_back(cc->EvalMult(ct, swap_vec_diagonal));
    res_martrix.push_back(cc->EvalMult(ct_lower, swap_vec_lower));

    res  = cc->EvalAdd(res_martrix[0], res_martrix[1]);
    res  = cc->EvalAdd(res, res_martrix[2]);
    end  = std::chrono::system_clock::now();
    time = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "mul and add time: " << time << "ms" << std::endl;
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

double bitonic_sort(int plain_bits, int num_slots) {
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
    SecretKeyDist secretKeyDist       = SPARSE_TERNARY;
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
    uint32_t levelsAvailableAfterBootstrap = 50;
    usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
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
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, precision);
    std::vector<double> input(length, 0);
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
    // cc->Decrypt(keyPair.secretKey, ciphertext_unsort, &plaintextDec);
    // plaintextDec->SetLength(num_slots);
    // std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
    // for (int i = 0; i < int(num_slots); i++) {
    //     std::cout << " Enc: " << enc[i].real() << std::endl;
    // }

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
            // std::cout << "bitonic_swap " << std::endl;
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
    return time_total / double(length) * num_slots;
}

double Eval_Agg(int plain_bits, int num_slots) {
    // std::cout << "Agg: " << plain_bits << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
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
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetRingDim(65536 * 2);

    usint multDepth = 2;
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    // cc->Enable(FHE);
    // cc->Enable(ADVANCEDSHE);

    // We need to enable Advanced SHE to use the Chebyshev approximation.

    double precision = (1 << (plain_bits - 1)) - 1;

    double lowerBound = -precision - 10;
    double upperBound = precision + 10;
    // double bound      = 3;
    auto keyPair                         = cc->KeyGen();
    const std::vector<DCRTPoly>& ckks_pk = keyPair.publicKey->GetPublicElements();
    // std::cout << "Moduli chain of pk: " << std::endl;
    // print_moduli_chain(ckks_pk[0]);
    usint ringDim = cc->GetRingDimension();
    // std::cout << "ringDim: " << ringDim << std::endl;
    int length = ringDim / 2;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, precision);
    std::vector<double> input_1(num_slots, 0), input_2(num_slots, 0);
    for (int i = 0; i < num_slots; i++) {
        input_1[i] = message(engine);
        input_2[i] = message(engine);
    }

    cc->EvalMultKeyGen(keyPair.secretKey);
    std::vector<int> rotstep;
    for (int i = 1; i < num_slots; i *= 2) {
        rotstep.push_back(i);
        rotstep.push_back(-i);
    }

    cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
    // cc->EvalBootstrapSetup(levelBudget);
    // cc->EvalBootstrapKeyGen(keyPair.secretKey, length);

    Plaintext plain_1 = cc->MakeCKKSPackedPlaintext(input_1);
    auto cipher_1     = cc->Encrypt(keyPair.publicKey, plain_1);
    Plaintext plain_2 = cc->MakeCKKSPackedPlaintext(input_2);
    auto cipher_2     = cc->Encrypt(keyPair.publicKey, plain_2);
    Plaintext plain_3 = cc->MakeCKKSPackedPlaintext(input_1);
    auto cipher_3     = cc->Encrypt(keyPair.publicKey, plain_3);
    Plaintext plain_4 = cc->MakeCKKSPackedPlaintext(input_2);
    auto cipher_4     = cc->Encrypt(keyPair.publicKey, plain_4);
    double time_total = 0;
    std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
    start                                       = std::chrono::system_clock::now();
    auto cipher_mul_1                           = cc->EvalMult(cipher_1, cipher_1);
    auto cipher_mul_2                           = cc->EvalMult(cipher_1, cipher_2);
    auto cipher_mul_3                           = cc->EvalMult(cipher_1, cipher_3);
    auto cipher_mul_4                           = cc->EvalMult(cipher_1, cipher_4);
    int logrow                                  = log2(num_slots);
    Ciphertext<lbcrypto::DCRTPoly> cipher_agg_1 = cipher_mul_1;
    Ciphertext<lbcrypto::DCRTPoly> cipher_agg_2 = cipher_mul_2;
    Ciphertext<lbcrypto::DCRTPoly> cipher_agg_3 = cipher_mul_3;
    Ciphertext<lbcrypto::DCRTPoly> cipher_agg_4 = cipher_mul_4;
    Ciphertext<lbcrypto::DCRTPoly> cipher_rotate_1, cipher_rotate_2, cipher_rotate_3, cipher_rotate_4;
    for (size_t i = 0; i < logrow; i++) {
        int step        = 1 << (logrow - i - 1);
        cipher_rotate_1 = cc->EvalRotate(cipher_mul_1, step);
        cipher_agg_1    = cc->EvalAdd(cipher_rotate_1, cipher_agg_1);

        step            = 1 << (logrow - i - 1);
        auto temp_2     = cipher_agg_2;
        cipher_rotate_2 = cc->EvalRotate(cipher_agg_2, step);
        cipher_agg_2    = cc->EvalAdd(cipher_rotate_2, cipher_agg_2);

        step            = 1 << (logrow - i - 1);
        auto temp_3     = cipher_agg_3;
        cipher_rotate_3 = cc->EvalRotate(cipher_mul_3, step);
        cipher_agg_3    = cc->EvalAdd(cipher_rotate_3, cipher_agg_3);

        step            = 1 << (logrow - i - 1);
        auto temp_4     = cipher_agg_4;
        cipher_rotate_4 = cc->EvalRotate(cipher_mul_4, step);
        cipher_agg_4    = cc->EvalAdd(cipher_rotate_4, cipher_agg_4);
    }
    end        = std::chrono::system_clock::now();
    time_total = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    // time_total = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

    return time_total;
}

double Eval_Agg_1(int plain_bits, int num_slots) {
    // std::cout << "Agg: " << plain_bits << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
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
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetRingDim(65536 * 2);

    usint multDepth = 2;
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    // cc->Enable(FHE);
    // cc->Enable(ADVANCEDSHE);

    // We need to enable Advanced SHE to use the Chebyshev approximation.

    double precision = (1 << (plain_bits - 1)) - 1;

    double lowerBound = -precision - 10;
    double upperBound = precision + 10;
    // double bound      = 3;
    auto keyPair                         = cc->KeyGen();
    const std::vector<DCRTPoly>& ckks_pk = keyPair.publicKey->GetPublicElements();
    // std::cout << "Moduli chain of pk: " << std::endl;
    // print_moduli_chain(ckks_pk[0]);
    usint ringDim = cc->GetRingDimension();
    std::cout << "ringDim: " << ringDim << std::endl;
    int length = ringDim / 2;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, precision);
    std::vector<double> input_1(num_slots, 0), input_2(num_slots, 0);
    for (int i = 0; i < num_slots; i++) {
        input_1[i] = message(engine);
        input_2[i] = message(engine);
    }

    cc->EvalMultKeyGen(keyPair.secretKey);
    std::vector<int> rotstep;
    for (int i = 1; i < num_slots; i *= 2) {
        rotstep.push_back(i);
        rotstep.push_back(-i);
    }

    cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
    // cc->EvalBootstrapSetup(levelBudget);
    // cc->EvalBootstrapKeyGen(keyPair.secretKey, length);

    Plaintext plain_1 = cc->MakeCKKSPackedPlaintext(input_1);
    auto cipher_1     = cc->Encrypt(keyPair.publicKey, plain_1);
    double time_total = 0;
    std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
    start                                       = std::chrono::system_clock::now();
    auto cipher_mul_1                           = cc->EvalMult(cipher_1, cipher_1);
    int logrow                                  = log2(num_slots);
    Ciphertext<lbcrypto::DCRTPoly> cipher_agg_1 = cipher_mul_1;
    Ciphertext<lbcrypto::DCRTPoly> cipher_rotate_1;
    for (size_t i = 0; i < logrow; i++) {
        int step        = 1 << (logrow - i - 1);
        cipher_rotate_1 = cc->EvalRotate(cipher_mul_1, step);
        cipher_agg_1    = cc->EvalAdd(cipher_rotate_1, cipher_agg_1);
    }
    end        = std::chrono::system_clock::now();
    time_total = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    // time_total = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

    return time_total;
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

    parameters.SetScalingModSize(scalingModSize);
    parameters.SetFirstModSize(firstModSize);
    uint32_t polyDegree = 0;
    uint32_t multDepth  = 23;
    // Choosing a higher degree yields better precision, but a longer runtime.
    if (plain_bits >= 8)
        polyDegree = 119;
    else {
        polyDegree = 27;
        multDepth  = 20;
    }
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

double Eval_modular_Strictly_greater_test(uint32_t plain_bits, uint32_t block) {
    // std::cout << plain_bits << " plain_bits " << block << " Strictly modular greater test " << std::endl;
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
    parameters.SetRingDim(65536);
    uint32_t polyDegree = 59;
    uint32_t multDepth  = 20 + block;
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
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    // std::uniform_int_distribution<int> message_part_16_1(-2048, -1);
    std::uniform_int_distribution<int> message(0, precision);
    usint ringDim = cc->GetRingDimension();
    int length    = ringDim / 2;

    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;
    // std::cout << "using polyDegree " << polyDegree << std::endl;
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
                // std::cout << "input_a:" << input_a_plain[i] << ", input_b:" << input_b_plain[i]
                //           << ", res :" << finalResult[i] << std::endl;
            }

            err += std::fabs(expectedOutput_compres[i] - finalResult[i].real());
            totalerr += err / encodedLength;
        }
        // totaltime +=
        //     double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / double(encodedLength));
        totaltime +=
            double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / double(encodedLength));
        if (wrong > 0) {
            // std::cout << "avg error: 2^" << std::log2(err / encodedLength) << std::endl;
            // std::cout << "wrong number:" << wrong << std::endl;
        }
    }
    // std::cout << encodedLength << "slots amortize time: " << double((totaltime / num_test)) << "ms" << std::endl
    //           << std::endl;
    // std::cout << "avg error: 2^" << std::log2(totalerr / num_test) << std::endl;
    return totaltime / num_test;
    // for (int j = 0; j < 8; j++) {
    //     std::cout << "input_a:" << input_a_plain[j] << ", input_b:" << input_b_plain[j] << ", res :" << finalResult[j]
    //               << std::endl;
    // }
}

double Eval_modular_Strictly_equal_test(uint32_t plain_bits, uint32_t block) {
    // std::cout << plain_bits << " plain_bits " << block << " Strictly modular equal test " << std::endl;
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
    parameters.SetRingDim(65536 * 2);
    // uint32_t polyDegree = 0;
    uint32_t polyDegree = 59;
    uint32_t multDepth  = 20 + block;
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

    // std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;
    // std::cout << "using polyDegree " << polyDegree << std::endl;
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

        end = std::chrono::system_clock::now();
        Plaintext plaintextDec, plaintextDec_comp_res, plaintextDec_comp_res_equal;

        cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec_comp_res);

        plaintextDec_comp_res->SetLength(encodedLength);
        std::vector<double> expectedOutput_compres;
        std::vector<std::complex<double>> finalResult = plaintextDec_comp_res->GetCKKSPackedValue();
        double err                                    = 0;
        int wrong                                     = 0;

        for (int i = 0; i < int(encodedLength); i++) {
            expectedOutput_compres.push_back(static_cast<double>(input_a_plain[i] == input_b_plain[i]));

            if (fabs(expectedOutput_compres[i] - finalResult[i].real()) > 0.5) {
                wrong++;
                // std::cout << "input_a:" << input_a_plain[i] << ", input_b:" << input_b_plain[i]
                //           << ", res :" << finalResult[i] << std::endl;
            }

            err += std::fabs(expectedOutput_compres[i] - finalResult[i].real());
            totalerr += err / encodedLength;
        }
        totaltime +=
            double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / double(encodedLength));

        if (wrong > 0) {
            // std::cout << "avg error: 2^" << std::log2(err / encodedLength) << std::endl;
            // std::cout << "wrong number:" << wrong << std::endl;
        }
    }
    // std::cout << encodedLength << "slots amortize time: " << double((totaltime / num_test)) << "ms" << std::endl
    //           << std::endl;
    // std::cout << "avg error: 2^" << std::log2(totalerr / num_test) << std::endl;
    return totaltime / num_test;
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
    uint32_t multDepth = 23;
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

void Eval_query(uint32_t plain_bits, uint32_t num, uint32_t block) {
    //filter

    double filter_time = 0.;
    //8*num comp
    for (int i = 0; i < 4; i++)
        filter_time += Eval_modular_Strictly_greater_test(plain_bits, block);
    for (int i = 0; i < 4; i++)
        filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    filter_time *= num;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;

    double agg_time = Eval_Agg(plain_bits, num);
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    double sort_time = bitonic_sort(plain_bits, num);
    std::cout << num << " sort_time time: " << sort_time << "ms" << std::endl;
    double total_time = filter_time + agg_time + sort_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q1(uint32_t plain_bits, uint32_t num, uint32_t block) {
    //filter

    double filter_time = 0.;
    //8*num comp
#pragma omp parallel for num_threads(4)
    for (int i = 0; i < 4; i++)
        filter_time += Eval_modular_Strictly_greater_test(plain_bits, block);
        // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
#pragma omp parallel for num_threads(4)
    for (int i = 0; i < 4; i++)
        filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    // std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    filter_time *= num;
    double agg_time = Eval_Agg(plain_bits, num);
    // std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;

    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q3(uint32_t plain_bits, uint32_t num, uint32_t block) {
    //filter
    std::cout << "tpc-q3" << std::endl;
    int comp_num = 5;
    int agg_num  = 1;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    //8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++) {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++) {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time          = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}
void Eval_query_simple_q4(uint32_t plain_bits, uint32_t num, uint32_t block) {
    //filter
    //filter
    std::cout << "tpc-q4" << std::endl;
    int comp_num = 5;
    int agg_num  = 1;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    //8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++) {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++) {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time          = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}
// without order by
void Eval_query_simple_q5(uint32_t plain_bits, uint32_t num, uint32_t block) {
    //filter
    //filter
    std::cout << "tpc-q5" << std::endl;
    int comp_num = 9;
    int agg_num  = 1;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    //8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++) {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++) {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time          = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q10(uint32_t plain_bits, uint32_t num, uint32_t block) {
    //filter
    std::cout << "tpc-q10" << std::endl;
    int comp_num = 6;
    int agg_num  = 1;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    //8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++) {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++) {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time          = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q12(uint32_t plain_bits, uint32_t num, uint32_t block) {
    //filter
    std::cout << "tpc-q12" << std::endl;
    int comp_num = 8;
    int agg_num  = 2;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    //8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++) {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++) {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time          = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q17(uint32_t plain_bits, uint32_t num, uint32_t block) {
    //filter
    std::cout << "tpc-q17" << std::endl;
    int comp_num = 5;
    int agg_num  = 3;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    //8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++) {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++) {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time          = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q18(uint32_t plain_bits, uint32_t num, uint32_t block) {
    //filter
    std::cout << "tpc-q18" << std::endl;
    int comp_num = 4;
    int agg_num  = 2;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    //8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++) {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++) {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time          = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q19(uint32_t plain_bits, uint32_t num, uint32_t block) {
    //filter
    std::cout << "tpc-q19" << std::endl;
    int comp_num = 28;
    int agg_num  = 1;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    //8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++) {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++) {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time          = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q21(uint32_t plain_bits, uint32_t num, uint32_t block) {
    //filter
    std::cout << "tpc-q21" << std::endl;
    int comp_num = 11;
    int agg_num  = 1;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    //8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++) {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++) {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time          = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q6(uint32_t plain_bits, uint32_t num, uint32_t block) {
    //filter
    std::cout << "tpc-q6 " << std::endl;
    std::vector<double> filter_times(5, 0);
    //8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(5)
        for (int i = 0; i < 5; i++) {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    double agg_time = 0.;
    for (int test = 0; test < num_test; test++)
        agg_time += Eval_Agg_1(plain_bits, num);
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time          = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

// void Eval_query_simple_q12(uint32_t plain_bits, uint32_t num, uint32_t block) {
//     //filter

//     double filter_time = 0.;
//     //8*num comp
//     for (int i = 0; i < 4; i++)
//         filter_time += Eval_modular_Strictly_greater_test(plain_bits, block);
//     for (int i = 0; i < 4; i++)
//         filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
//     filter_time *= num;
//     std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;

//     double total_time = filter_time;
//     std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
// }

int main(int argc, char* argv[]) {
    Eval_query_simple_q3(8, 32768, 2);
    Eval_query_simple_q4(8, 32768, 2);
    Eval_query_simple_q5(8, 32768, 2);
    Eval_query_simple_q10(8, 32768, 2);
    Eval_query_simple_q12(8, 32768, 2);
    Eval_query_simple_q17(8, 32768, 2);
    Eval_query_simple_q18(8, 32768, 2);
    Eval_query_simple_q19(8, 32768, 2);
    Eval_query_simple_q21(8, 32768, 2);
    // for (int i = 512; i < 32769; i *= 2) {
    //     // Eval_query_simple_q12(8, i, 2);
    //     // double agg_time = Eval_Agg(8, i);
    //     // std::cout << i << " agg_time time: " << agg_time / 1000 << "ms" << std::endl;
    //     Eval_query_simple_q1(8, i, 2);
    // }
    // Eval_query_simple_q1(8, 65536, 2);
    // Eval_query_simple_q6(8, 32768, 2);
    // // Eval_query(8, 64, 2);
    // int count         = 0;
    // int time_values[] = {186953, 190282, 193382, 201848, 215466, 219006, 223919, 223253, 224660, 225923,
    //                      226809, 227132, 226943, 228712, 231627, 232627, 233570, 233402, 233920, 237331,
    //                      240156, 240210, 240908, 241014, 240836, 243348, 247597, 253094, 252720, 254240,
    //                      256164, 257420, 257090, 261306, 262393, 261914, 264352, 264092, 264408, 264926,
    //                      266192, 266980, 267917, 268344, 268302, 268848, 269906, 270678};
    // for (int i = 0; i < 48; i++) {
    //     count += time_values[i];
    // }
    // std::cout << count << std::endl;

    return 0;
}
