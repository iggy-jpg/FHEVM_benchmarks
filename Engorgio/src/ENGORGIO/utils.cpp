#include "comp.h"
#include "quant.h"
#include <chrono>
#include "math/chebyshev.h"
#include <stdexcept>
namespace openfhe
{
    void output(Ciphertext<lbcrypto::DCRTPoly> &res, int slot, lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        Plaintext plaintextDec;
        auto cc = res->GetCryptoContext();
        cc->Decrypt(privateKey, res, &plaintextDec);
        std::vector<std::complex<double>> dec = plaintextDec->GetCKKSPackedValue();
        for (int i = 0; i < slot; i++)
        {
            std::cout << "dec: " << dec[i].real() << std::endl;
        }
    }

    void inverse(Ciphertext<lbcrypto::DCRTPoly> &res)
    {
        auto cc = res->GetCryptoContext();
        auto neg = cc->EvalNegate(res);

        res = cc->EvalAdd(1, res);
    }

    double error_estimate(std::vector<double> plain, Ciphertext<lbcrypto::DCRTPoly> &res,
                          lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey, size_t num_slots)
    {
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
        double err_max = 0.;
        // std::cout << "plain.size:" << plain.size() << std::endl;
        for (int i = 0; i < int(plain.size()); i++)
        {
            err_max = fabs(plain[i] - Result[i].real()) > err_max ? fabs(plain[i] - Result[i].real()) : err_max;
            err += fabs(plain[i] - Result[i].real());
            error_vec.push_back(fabs(plain[i] - plaintextDec->GetCKKSPackedValue()[i].real()));
            if (fabs(plain[i] - Result[i].real()) > 0.5)
                throw std::runtime_error(" err>0.5  ");
        }
        // int maxPosition = max_element(error_vec.begin(), error_vec.end()) - error_vec.begin();
        std::cout << "avg error:  " << err / int(plain.size()) << " ~ 2^" << std::log2(err / int(plain.size())) << std::endl;
        return err_max;
    }

    double error_estimate_modular(std::vector<double> plain, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &res,
                                  lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey, size_t num_slots)
    {
        auto cc = res[0]->GetCryptoContext();
        std::vector<double> error_vec;
        std::vector<double> merge(num_slots, 0);
        Plaintext plaintextDec0, plaintextDec1;
        cc->Decrypt(privateKey, res[0], &plaintextDec0);
        std::vector<std::complex<double>> Result0 = plaintextDec0->GetCKKSPackedValue();
        cc->Decrypt(privateKey, res[1], &plaintextDec1);
        std::vector<std::complex<double>> Result1 = plaintextDec1->GetCKKSPackedValue();
        for (int j = 0; j < int(Result0.size()); j++)
        {
            merge[j] += round(Result0[j].real()) + round(Result1[j].real()) * pow(2, 8);
        }
        double err = 0.;
        double err_max = 0.;
        for (int i = 0; i < int(plain.size()); i++)
        {
            err_max = fabs((int(plain[i]) >> 8) - Result1[i].real()) > err_max ? fabs((int(plain[i]) >> 8) - Result1[i].real()) : err_max;
            err += fabs((int(plain[i]) >> 8) - Result1[i].real());
            error_vec.push_back(fabs(int(plain[i]) >> 8) - Result1[i].real());

            if (error_vec[i] > 0.5)
                std::cout << "error! plain:  " << plain[i] << "  merge: " << merge[i] << std::endl;
        }
        std::cout << "avg error:  " << err / int(plain.size()) << " ~ 2^" << std::log2(err / int(plain.size())) << std::endl;
        return err_max;
    }

    double error_estimate_unbounded(std::vector<double> plain, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &res,
                                    lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        // std::cout << "error_estimate_unbounded" << std::endl;
        auto cc = res[0]->GetCryptoContext();
        std::vector<double> error_vec;
        std::vector<double> dec_res;
        Plaintext plaintextDec;
        for (int i = 0; i < res.size(); i++)
        {
            cc->Decrypt(privateKey, res[i], &plaintextDec);
            std::vector<std::complex<double>> Result = plaintextDec->GetCKKSPackedValue();
            auto real_res = extractRealParts(Result);
            dec_res.insert(dec_res.end(), real_res.begin(), real_res.end());
        }

        double err = 0.;
        double err_max = 0.;
        // std::cout << "plain.size:" << plain.size() << std::endl;
        for (int i = 0; i < int(plain.size()); i++)
        {
            err_max = fabs(plain[i] - dec_res[i]) > err_max ? fabs(plain[i] - dec_res[i]) : err_max;
            err += fabs(plain[i] - dec_res[i]);
            error_vec.push_back(fabs(plain[i] - dec_res[i]));
            if (fabs(plain[i] - dec_res[i]) > 0.5)
                // throw std::runtime_error(" err>0.5  ");
                std::cout << "Result dec: " << dec_res[i] << ", expect res: " << plain[i] << std::endl;
        }
        // int maxPosition = max_element(error_vec.begin(), error_vec.end()) - error_vec.begin();
        std::cout << "avg error:  " << err / int(plain.size()) << " ~ 2^" << std::log2(err / int(plain.size())) << std::endl;
        return err_max;
    }

    double error_estimate_unbounded_modular(std::vector<double> plain, std::vector<std::vector<Ciphertext<lbcrypto::DCRTPoly>>> &res,
                                            lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey, int num_ct, int num_slots, int length)
    {
        // std::cout << "error_estimate_unbounded" << std::endl;
        auto cc = res[0][0]->GetCryptoContext();
        std::vector<double> error_vec;
        std::vector<double> dec_res;
        Plaintext plaintextDec;
        std::vector<double> merge0(num_slots), merge1(num_slots);
        std::vector<Plaintext> plaintextDec0(num_ct), plaintextDec1(num_ct);
        cc->Decrypt(privateKey, res[0][0], &plaintextDec0[0]);
        cc->Decrypt(privateKey, res[1][0], &plaintextDec0[1]);

        cc->Decrypt(privateKey, res[0][1], &plaintextDec1[0]);
        cc->Decrypt(privateKey, res[1][1], &plaintextDec1[1]);
        std::vector<std::complex<double>> Result00 = plaintextDec0[0]->GetCKKSPackedValue();
        std::vector<std::complex<double>> Result10 = plaintextDec0[1]->GetCKKSPackedValue();
        std::vector<std::complex<double>> Result01 = plaintextDec1[0]->GetCKKSPackedValue();
        std::vector<std::complex<double>> Result11 = plaintextDec1[1]->GetCKKSPackedValue();
        for (int j = 0; j < int(Result00.size()); j++)
        {
            merge0[j] += round(Result00[j].real()) + round(Result01[j].real()) * pow(2, 8);
            merge1[j] += round(Result10[j].real()) + round(Result11[j].real()) * pow(2, 8);
        }
        merge0.insert(merge0.end(), merge1.begin(), merge1.end());
        double err = 0.;
        double err_max = 0.;

        // std::cout << "plain.size:" << plain.size() << std::endl;
        for (int i = 0; i < int(plain.size()); i++)
        {
            err_max = fabs(plain[i] - merge0[i]) > err_max ? fabs(plain[i] - merge0[i]) : err_max;
            err += fabs(plain[i] - merge0[i]);
            error_vec.push_back(fabs(plain[i] - merge0[i]));
            if (fabs(plain[i] - merge0[i]) > 0.5)
                // throw std::runtime_error(" err>0.5  ");
                std::cout << "Result dec: " << merge0[i] << ", expect res: " << plain[i] << std::endl;
        }
        // int maxPosition = max_element(error_vec.begin(), error_vec.end()) - error_vec.begin();
        std::cout << "avg error:  " << err / int(plain.size()) << " ~ 2^" << std::log2(err / int(plain.size())) << std::endl;
        return err_max;
    }

    std::vector<int> dot_product(std::vector<int> &A, std::vector<int> &B)
    {
        std::vector<int> result(B.size());
        for (int i = 0; i < A.size(); i++)
        {
            result[i] += A[i] * B[i];
        }
        return result;
    }

    std::vector<int> vadd(std::vector<int> &A, std::vector<int> &B)
    {
        std::vector<int> result(B.size());
        for (int i = 0; i < A.size(); i++)
        {
            result[i] += A[i] + B[i];
        }
        return result;
    }

    void cross_rotate(int N, int l, int T)
    {
        std::vector<std::vector<int>> big_vec;
        std::vector<std::vector<int>> big_vec_rot;
        std::vector<int> test;
        // int N     = 8;
        int baby = T % N;
        int giant = floor(T / N);
        for (int i = 0; i < l; i++)
        {
            std::vector<int> small_vec;
            for (int j = 0; j < N; j++)
            {
                small_vec.push_back(i * N + j);
                test.push_back(i * N + j);
            }
            big_vec.push_back(small_vec);
        }
        for (int i = 0; i < test.size(); i++)
            std::cout << test[i] << ", " << std::endl;
        std::vector<std::vector<int>> big_vec_prerot;
        for (int i = 0; i < l; i++)
        {
            big_vec_prerot.push_back(big_vec[(i + giant) % l]);
        }
        std::vector<int> m0(N, 0);
        std::vector<int> m1(N, 1);
        for (int i = 0; i < N - baby; i++)
        {
            m0[i] = 1;
            m1[i] = 0;
        }
        for (int i = 0; i < l - 1; i++)
        {
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
        for (int i = 0; i < l; i++)
        {
            for (int j = 0; j < N; j++)
            {
                // if (big_vec_rot[i][j] != test[i * l + j])
                std::cout << "i: " << i << ", j: " << j << ", big_vec_rot: " << big_vec_rot[i][j]
                          << ", test:" << test[i * N + j] << std::endl;
            }
        }
    }
    // N slots length, l vector number, T rotate step
    std::vector<Ciphertext<lbcrypto::DCRTPoly>> gen_rotate(int N, int l, int T,
                                                           std::vector<Ciphertext<lbcrypto::DCRTPoly>> &big_vec,
                                                           std::vector<double> &test,
                                                           lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        // std::cout << "N:" << N << ", l: " << l << ", T:" << T << std::endl;
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> big_vec_rot;
        // std::cout << "gen_rotate start" << std::endl;
        if (test.size() < l * N)
        {
            std::cerr << "Error: test vector size is insufficient." << std::endl;
            return {};
        }
        // int N     = 8;
        int baby = T % N;
        auto cc = big_vec[0]->GetCryptoContext();

        int giant = floor(std::abs(T) / N);
        for (int i = 0; i < l; i++)
        {
            Plaintext plaintextDec;
            cc->Decrypt(privateKey, big_vec[i], &plaintextDec);
            std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
            // for (int j = 0; j < 20; j++)
            // {
            //     std::cout << enc[j].real() << std::endl;
            // }
            for (int j = 0; j < N; j++)
            {
                if (round(enc[j].real()) - test[i * N + j] > 0.5)
                    std::cout << "i: " << i << ", j: " << j << ", big_vec_rot: " << enc[j].real()
                              << ", test:" << test[i * N + j] << std::endl;
            }
        }

        std::vector<double> plain_rot(test.size());
        std::rotate_copy(test.begin(), test.begin() + T, test.end(), plain_rot.begin());
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> big_vec_prerot;
        for (int i = 0; i < l; i++)
        {
            // std::cout << "i: " << i << ", (i + giant) % l: " << (i + giant) % l << std::endl;
            big_vec_prerot.push_back(big_vec[(i + giant) % l]);
        }

        std::vector<double> mask0(N, 0);
        std::vector<double> mask1(N, 1);
        // std::cout << "assign mask" << std::endl;
        for (int i = 0; i < N - baby; i++)
        {
            mask0[i] = 1;
            mask1[i] = 0;
        }
        Plaintext m0 = cc->MakeCKKSPackedPlaintext(mask0);
        Plaintext m1 = cc->MakeCKKSPackedPlaintext(mask1);
        if (baby != 0)
        {
            for (int i = 0; i < l - 1; i++)
            {
                auto v0 = cc->EvalRotate(big_vec_prerot[i], baby);
                auto v1 = cc->EvalRotate(big_vec_prerot[i + 1], baby);
                auto rm0 = cc->EvalMult(v0, m0);
                auto rm1 = cc->EvalMult(v1, m1);
                big_vec_rot.push_back(cc->EvalAdd(rm0, rm1));
            }

            auto v0 = cc->EvalRotate(big_vec_prerot[l - 1], baby);
            auto v1 = cc->EvalRotate(big_vec_prerot[0], baby);
            auto rm0 = cc->EvalMult(v0, m0);
            auto rm1 = cc->EvalMult(v1, m1);
            big_vec_rot.push_back(cc->EvalAdd(rm0, rm1));
        }
        else
        {
            big_vec_rot = big_vec_prerot;
        }
        Plaintext plaintextDec;
        for (int i = 0; i < l; i++)
        {
            cc->Decrypt(privateKey, big_vec_rot[i], &plaintextDec);
            std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
            for (int j = 0; j < N; j++)
            {
                if (round(enc[j].real()) - plain_rot[i * N + j] > 0.5)
                    std::cout << "i: " << i << ", j: " << j << ", big_vec_rot: " << enc[i].real()
                              << ", plain_rot:" << plain_rot[i * N + j] << std::endl;
            }
        }
        // std::cout << " finish " << std::endl;
        return big_vec_rot;
    }

    void test_gen_rotate(int plain_bits, int l, int num_slots)
    {
        std::cout << "\ntest_gen_rotate: " << plain_bits << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;

        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
        ScalingTechnique rescaleTech = FIXEDAUTO;
        usint dcrtBits = 78;
        usint firstMod = 89;
#else
        usint dcrtBits = 59;
        usint firstMod = 60;
#endif
        parameters.SetScalingModSize(dcrtBits);
        parameters.SetFirstModSize(firstMod);
        std::vector<std::uint32_t> levelBudget = {4, 4};
        // SecretKeyDist secretKeyDist       = UNIFORM_TERNARY;
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetSecretKeyDist(secretKeyDist);
        // Choosing a higher degree yields better precision, but a longer runtime.
        std::uint32_t polyDegree = 495; //
        std::uint32_t numIterations = 2;
        std::uint32_t levelsAvailableAfterBootstrap = 10;
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
        auto keyPair = cc->KeyGen();
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        std::vector<std::vector<double>> input_big;
        std::vector<double> test;
        std::vector<Plaintext> plain_big(l);
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> cipher_big(l);

        for (int i = 0; i < l; i++)
        {
            std::vector<double> input_temp;
            for (int j = 0; j < length; j++)
            {
                input_temp.push_back(j);
                test.push_back(input_temp[j]);
            }
            input_big.push_back(input_temp);
            plain_big[i] = cc->MakeCKKSPackedPlaintext(input_temp);
            cipher_big[i] = cc->Encrypt(keyPair.publicKey, plain_big[i]);
        }
        std::vector<std::vector<double>> plain_sort_vec;
        std::cout << "ringDim: " << ringDim << std::endl;
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::vector<int> rotstep;
        for (int i = 1; i < num_slots; i *= 2)
        {
            rotstep.push_back(i);
            rotstep.push_back(-i);
        }

        cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_rotate =
            gen_rotate(length, l, 16, cipher_big, test, keyPair.secretKey);
    }

    double Eval_Agg_4(int plain_bits, int num_slots)
    {
        CCParams<CryptoContextCKKSRNS> parameters;
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
        ScalingTechnique rescaleTech = FIXEDAUTO;
        usint dcrtBits = 78;
        usint firstMod = 89;
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

        double precision = (1 << (plain_bits - 1)) - 1;

        double lowerBound = -precision - 10;
        double upperBound = precision + 10;
        // double bound      = 3;
        auto keyPair = cc->KeyGen();
        const std::vector<DCRTPoly> &ckks_pk = keyPair.publicKey->GetPublicElements();
        // std::cout << "Moduli chain of pk: " << std::endl;
        // print_moduli_chain(ckks_pk[0]);
        usint ringDim = cc->GetRingDimension();
        // std::cout << "ringDim: " << ringDim << std::endl;
        int length = ringDim / 2;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        std::vector<double> input_1(num_slots, 0), input_2(num_slots, 0);
        for (int i = 0; i < num_slots; i++)
        {
            input_1[i] = message(engine);
            input_2[i] = message(engine);
        }

        cc->EvalMultKeyGen(keyPair.secretKey);
        std::vector<int> rotstep;
        for (int i = 1; i < num_slots; i *= 2)
        {
            rotstep.push_back(i);
            rotstep.push_back(-i);
        }

        cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
        // cc->EvalBootstrapSetup(levelBudget);
        // cc->EvalBootstrapKeyGen(keyPair.secretKey, length);

        Plaintext plain_1 = cc->MakeCKKSPackedPlaintext(input_1);
        auto cipher_1 = cc->Encrypt(keyPair.publicKey, plain_1);
        Plaintext plain_2 = cc->MakeCKKSPackedPlaintext(input_2);
        auto cipher_2 = cc->Encrypt(keyPair.publicKey, plain_2);
        Plaintext plain_3 = cc->MakeCKKSPackedPlaintext(input_1);
        auto cipher_3 = cc->Encrypt(keyPair.publicKey, plain_3);
        Plaintext plain_4 = cc->MakeCKKSPackedPlaintext(input_2);
        auto cipher_4 = cc->Encrypt(keyPair.publicKey, plain_4);
        double time_total = 0;
        std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
        start = std::chrono::system_clock::now();
        auto cipher_mul_1 = cc->EvalMult(cipher_1, cipher_1);
        auto cipher_mul_2 = cc->EvalMult(cipher_1, cipher_2);
        auto cipher_mul_3 = cc->EvalMult(cipher_1, cipher_3);
        auto cipher_mul_4 = cc->EvalMult(cipher_1, cipher_4);
        int logrow = log2(num_slots);
        Ciphertext<lbcrypto::DCRTPoly> cipher_agg_1 = cipher_mul_1;
        Ciphertext<lbcrypto::DCRTPoly> cipher_agg_2 = cipher_mul_2;
        Ciphertext<lbcrypto::DCRTPoly> cipher_agg_3 = cipher_mul_3;
        Ciphertext<lbcrypto::DCRTPoly> cipher_agg_4 = cipher_mul_4;
        Ciphertext<lbcrypto::DCRTPoly> cipher_rotate_1, cipher_rotate_2, cipher_rotate_3, cipher_rotate_4;
        for (size_t i = 0; i < logrow; i++)
        {
            int step = 1 << (logrow - i - 1);
            cipher_rotate_1 = cc->EvalRotate(cipher_mul_1, step);
            cipher_agg_1 = cc->EvalAdd(cipher_rotate_1, cipher_agg_1);

            step = 1 << (logrow - i - 1);
            auto temp_2 = cipher_agg_2;
            cipher_rotate_2 = cc->EvalRotate(cipher_agg_2, step);
            cipher_agg_2 = cc->EvalAdd(cipher_rotate_2, cipher_agg_2);

            step = 1 << (logrow - i - 1);
            auto temp_3 = cipher_agg_3;
            cipher_rotate_3 = cc->EvalRotate(cipher_mul_3, step);
            cipher_agg_3 = cc->EvalAdd(cipher_rotate_3, cipher_agg_3);

            step = 1 << (logrow - i - 1);
            auto temp_4 = cipher_agg_4;
            cipher_rotate_4 = cc->EvalRotate(cipher_mul_4, step);
            cipher_agg_4 = cc->EvalAdd(cipher_rotate_4, cipher_agg_4);
        }
        end = std::chrono::system_clock::now();
        time_total = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
        // time_total = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

        return time_total;
    }

    double Eval_Agg_1(int plain_bits, int num_slots)
    {
        // std::cout << "Agg: " << plain_bits << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
        ScalingTechnique rescaleTech = FIXEDAUTO;
        usint dcrtBits = 78;
        usint firstMod = 89;
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

        double precision = (1 << (plain_bits - 1)) - 1;

        double lowerBound = -precision - 10;
        double upperBound = precision + 10;
        // double bound      = 3;
        auto keyPair = cc->KeyGen();
        const std::vector<DCRTPoly> &ckks_pk = keyPair.publicKey->GetPublicElements();
        usint ringDim = cc->GetRingDimension();
        cc->EvalMultKeyGen(keyPair.secretKey);
        int length = ringDim / 2;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        std::vector<double> input_1(num_slots, 0), input_2(num_slots, 0);
        for (int i = 0; i < num_slots; i++)
        {
            input_1[i] = message(engine);
            input_2[i] = message(engine);
        }

        cc->EvalMultKeyGen(keyPair.secretKey);
        std::vector<int> rotstep;
        for (int i = 1; i < num_slots; i *= 2)
        {
            rotstep.push_back(i);
            rotstep.push_back(-i);
        }

        cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);

        Plaintext plain_1 = cc->MakeCKKSPackedPlaintext(input_1);
        auto cipher_1 = cc->Encrypt(keyPair.publicKey, plain_1);
        double time_total = 0;
        std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
        start = std::chrono::system_clock::now();
        auto cipher_mul_1 = cc->EvalMult(cipher_1, cipher_1);
        int logrow = log2(num_slots);
        Ciphertext<lbcrypto::DCRTPoly> cipher_agg_1 = cipher_mul_1;
        Ciphertext<lbcrypto::DCRTPoly> cipher_rotate_1;
        for (size_t i = 0; i < logrow; i++)
        {
            int step = 1 << (logrow - i - 1);
            cipher_rotate_1 = cc->EvalRotate(cipher_mul_1, step);
            cipher_agg_1 = cc->EvalAdd(cipher_rotate_1, cipher_agg_1);
        }
        end = std::chrono::system_clock::now();
        time_total = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
        // cc->ClearEvalSumKeys();
        // cc->ClearEvalMultKeys();
        // cc->ClearEvalAutomorphismKeys();
        // CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
        return time_total;
    }

    double Eval_SUM(int plain_bits, int num_slots)
    {
        // std::cout << "Agg: " << plain_bits << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
        ScalingTechnique rescaleTech = FIXEDAUTO;
        usint dcrtBits = 78;
        usint firstMod = 89;
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
        auto keyPair = cc->KeyGen();
        const std::vector<DCRTPoly> &ckks_pk = keyPair.publicKey->GetPublicElements();
        // std::cout << "Moduli chain of pk: " << std::endl;
        // print_moduli_chain(ckks_pk[0]);
        usint ringDim = cc->GetRingDimension();
        std::cout << "ringDim: " << ringDim << std::endl;
        int length = ringDim / 2;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        std::vector<double> input_1(num_slots, 0), input_2(num_slots, 0);
        for (int i = 0; i < num_slots; i++)
        {
            input_1[i] = message(engine);
            input_2[i] = message(engine);
        }

        cc->EvalMultKeyGen(keyPair.secretKey);
        std::vector<int> rotstep;
        for (int i = 1; i < num_slots; i *= 2)
        {
            rotstep.push_back(i);
            rotstep.push_back(-i);
        }

        cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
        // cc->EvalBootstrapSetup(levelBudget);
        // cc->EvalBootstrapKeyGen(keyPair.secretKey, length);

        Plaintext plain_1 = cc->MakeCKKSPackedPlaintext(input_1);
        auto cipher_1 = cc->Encrypt(keyPair.publicKey, plain_1);
        double time_total = 0;
        std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
        start = std::chrono::system_clock::now();
        auto cipher_mul_1 = cc->EvalMult(cipher_1, cipher_1);
        int logrow = log2(num_slots);
        Ciphertext<lbcrypto::DCRTPoly> cipher_agg_1 = cipher_mul_1;
        Ciphertext<lbcrypto::DCRTPoly> cipher_rotate_1;
        for (size_t i = 0; i < logrow; i++)
        {
            int step = 1 << (logrow - i - 1);
            cipher_rotate_1 = cc->EvalRotate(cipher_mul_1, step);
            cipher_agg_1 = cc->EvalAdd(cipher_rotate_1, cipher_agg_1);
        }
        end = std::chrono::system_clock::now();
        time_total = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
        // time_total = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

        return time_total;
    }

    std::vector<std::vector<double>> read_vectors_from_csv(const std::string &filename, int N)
    {
        std::vector<std::vector<double>> sift;
        std::ifstream file(filename);

        if (!file.is_open())
        {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            return sift;
        }

        std::string line;
        int count = 0;

        while (std::getline(file, line) && count < N)
        {
            std::vector<double> vector;
            std::stringstream ss(line);
            std::string value;

            while (std::getline(ss, value, ','))
            {
                vector.push_back(std::stod(value));
            }

            sift.push_back(vector);
            count++;
        }

        file.close();
        return sift;
    }

    std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> &matrix)
    {
        if (matrix.empty() || matrix[0].empty())
        {
            return {};
        }

        int rows = matrix.size();
        int cols = matrix[0].size();

        std::vector<std::vector<double>> transposed(cols, std::vector<double>(rows));

        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                transposed[j][i] = matrix[i][j];
            }
        }

        return transposed;
    }

    Ciphertext<lbcrypto::DCRTPoly> homdis(
        Ciphertext<lbcrypto::DCRTPoly> &ct0, Ciphertext<lbcrypto::DCRTPoly> &ct1)
    {
        auto cc = ct0->GetCryptoContext();
        auto sub_res = cc->EvalSub(ct0, ct1);
        auto dis_res = cc->EvalSquare(sub_res);
        return dis_res;
    }

    double Euclid_distance(int plain_bits)
    {
        CCParams<CryptoContextCKKSRNS> parameters;
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
        ScalingTechnique rescaleTech = FIXEDAUTO;
        usint dcrtBits = 78;
        usint firstMod = 89;
#else
        usint dcrtBits = 40;
        usint firstMod = 60;
#endif
        parameters.SetScalingModSize(dcrtBits);
        parameters.SetFirstModSize(firstMod);
        parameters.SetSecretKeyDist(secretKeyDist);
        parameters.SetRingDim(65536 * 2);

        usint multDepth = 20;
        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);

        double precision = (1 << (plain_bits - 1)) - 1;

        double lowerBound = -precision - 10;
        double upperBound = precision + 10;
        // double bound      = 3;
        auto keyPair = cc->KeyGen();
        cc->EvalMultKeyGen(keyPair.secretKey);
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;
        int vec_size = 128;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        std::vector<std::vector<double>> input_vec(vec_size, std::vector<double>(length)), read_sift, db_feature;
        std::string filename = "../sift_base_100k.csv";
        read_sift = read_vectors_from_csv(filename, length);
        for (int i = 0; i < vec_size; i++)
        {
            input_vec[i][0] = message(engine);
            for (int j = 1; j < length; j++)
            {
                input_vec[i][j] = input_vec[i][0];
            }
        }
        db_feature = transpose(read_sift);
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> edb_feature, enc_input;
        for (int i = 0; i < vec_size; i++)
        {
            Plaintext plain_feat = cc->MakeCKKSPackedPlaintext(db_feature[i]);
            edb_feature.push_back(cc->Encrypt(keyPair.publicKey, plain_feat));
            Plaintext plain_vec = cc->MakeCKKSPackedPlaintext(input_vec[i]);
            enc_input.push_back(cc->Encrypt(keyPair.publicKey, plain_vec));
        }
        double time_total = 0;
        std::chrono::system_clock::time_point start, end;
        start = std::chrono::system_clock::now();
        auto enc_dis = homdis(enc_input[0], edb_feature[0]);
        for (int i = 1; i < vec_size; i++)
        {
            enc_dis = cc->EvalAdd(enc_dis, homdis(enc_input[i], edb_feature[i]));
        }
        end = std::chrono::system_clock::now();
        time_total = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
        cc->ClearEvalSumKeys();
        cc->ClearEvalMultKeys();
        cc->ClearEvalAutomorphismKeys();
        CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
        return time_total;
    }
}