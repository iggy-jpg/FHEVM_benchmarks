#include "utils.h"
#include "comp.h"
#include <chrono>
#include "math/chebyshev.h"
#include <stdexcept>
namespace openfhe
{

    void bitonic_swap_unbounded(int stage, int part, int slots, int length, Ciphertext<lbcrypto::DCRTPoly> &compres,
                                std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ct, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &res, lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        int num_ct = length / slots;
        auto cc = ct[0]->GetCryptoContext();
        // int Block = 1 << (stage + 1);
        int block = 1 << (stage + 1 - part);
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> res_matrix_diag(2);
        std::vector<double> eliminate(length, 0), eliminate_rot(length, 0), test_ct, test_compres, test_vec_lower;

        auto compres_neg = cc->EvalNegate(compres);
        auto swap_vec_lower = cc->EvalAdd(compres_neg, 1);
        auto res_matrix_upper = cc->EvalMult(ct[0], swap_vec_lower);
        auto res_matrix_lower = cc->EvalMult(ct[1], swap_vec_lower);
        for (int i = 0; i < num_ct; i++)
        {
            res_matrix_diag[i] = cc->EvalMult(ct[i], compres);
        }

        res[0] = cc->EvalAdd(res_matrix_diag[0], res_matrix_lower);
        res[1] = cc->EvalAdd(res_matrix_diag[1], res_matrix_upper);
    }

    void bitonic_swap_unbounded_modular(int stage, int part, int slots, int length, Ciphertext<lbcrypto::DCRTPoly> &compres,
                                        std::vector<std::vector<Ciphertext<lbcrypto::DCRTPoly>>> &ct,
                                        std::vector<std::vector<Ciphertext<lbcrypto::DCRTPoly>>> &res,
                                        lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey, int blocks)
    {
        int num_ct = length / slots;
        auto cc = ct[0][0]->GetCryptoContext();
        // int Block = 1 << (stage + 1);
        int block = 1 << (stage + 1 - part);
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> res_matrix_diag(2);
        std::vector<double> eliminate(length, 0), eliminate_rot(length, 0), test_ct, test_compres, test_vec_lower;

        auto compres_neg = cc->EvalNegate(compres);
        auto swap_vec_lower = cc->EvalAdd(compres_neg, 1);
        for (int i = 0; i < blocks; i++)
        {

            auto res_matrix_upper = cc->EvalMult(ct[0][i], swap_vec_lower);
            auto res_matrix_lower = cc->EvalMult(ct[1][i], swap_vec_lower);
            for (int j = 0; j < num_ct; j++)
            {
                res_matrix_diag[j] = cc->EvalMult(ct[j][i], compres);
            }

            res[0][i] = cc->EvalAdd(res_matrix_diag[0], res_matrix_lower);
            res[1][i] = cc->EvalAdd(res_matrix_diag[1], res_matrix_upper);
        }
    }

    void bitonic_swap_topk(int slots, int topk_i, int k, Ciphertext<lbcrypto::DCRTPoly> &compres,
                           Ciphertext<lbcrypto::DCRTPoly> &ct, Ciphertext<lbcrypto::DCRTPoly> &res,
                           lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        std::chrono::system_clock::time_point start, end;
        double time = 0.;
        start = std::chrono::system_clock::now();
        auto cc = ct->GetCryptoContext();
        int Block = (1 << int(std::log2(k) + 1)) + (pow(2, topk_i + 1) - 1) * 2 * k;
        int block = 1 << (int(std::log2(k)) + 1);
        int step = (1 << int(std::log2(k))) + (pow(2, topk_i + 1) - 1) * 2 * k;
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> res_martrix;
        std::vector<double> eliminate(slots, 0), eliminate_rot(slots, 0);
        // 1level
        for (int i = 0; i < slots; i += Block)
        {
            for (int j = 0; j < block / 2; j++)
            {
                eliminate[i + j] = 1;
            }
        }

        Plaintext plaintextDec;
        Plaintext plain_eliminate = cc->MakeCKKSPackedPlaintext(eliminate);
        Plaintext plain_eliminate_rot = cc->MakeCKKSPackedPlaintext(eliminate_rot);
        auto ct_upper = cc->EvalRotate(ct, -step);
        auto ct_lower = cc->EvalRotate(ct, step);

        std::vector<std::complex<double>> enc;
        auto compres_neg = cc->EvalNegate(compres);
        auto compres_rot = cc->EvalRotate(compres, -step);

        auto swap_vec_lower = cc->EvalAdd(compres_neg, plain_eliminate);

        auto swap_vec_upper = cc->EvalRotate(swap_vec_lower, -step);

        auto swap_vec_diagonal = cc->EvalAdd(compres, compres_rot);
        start = std::chrono::system_clock::now();
        res_martrix.push_back(cc->EvalMult(ct_upper, swap_vec_upper));
        res_martrix.push_back(cc->EvalMult(ct, swap_vec_diagonal));
        res_martrix.push_back(cc->EvalMult(ct_lower, swap_vec_lower));
        res = cc->EvalAdd(res_martrix[0], res_martrix[1]);
        res = cc->EvalAdd(res, res_martrix[2]);
    }

    void bitonic_swap(int stage, int part, int slots, Ciphertext<lbcrypto::DCRTPoly> &compres,
                      Ciphertext<lbcrypto::DCRTPoly> &ct, Ciphertext<lbcrypto::DCRTPoly> &res,
                      lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        std::chrono::system_clock::time_point start, end;
        double time = 0.;
        start = std::chrono::system_clock::now();
        auto cc = ct->GetCryptoContext();
        // int Block = 1 << (stage + 1);
        int block = 1 << (stage + 1 - part);
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> res_martrix;
        std::vector<double> eliminate(slots, 0), eliminate_rot(slots, 0);
        // 1level
        for (int i = 0; i < slots; i += block)
        {
            for (int j = 0; j < block / 2; j++)
            {
                eliminate[i + j] = 1;
            }
        }
        std::rotate_copy(eliminate.begin(), eliminate.begin() + block / 2, eliminate.end(), eliminate_rot.begin());

        Plaintext plaintextDec;
        Plaintext plain_eliminate = cc->MakeCKKSPackedPlaintext(eliminate);
        end = std::chrono::system_clock::now();
        time = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

        auto ct_upper = cc->EvalRotate(ct, -block / 2);
        auto ct_lower = cc->EvalRotate(ct, block / 2);

        std::vector<std::complex<double>> enc;
        start = std::chrono::system_clock::now();
        auto compres_neg = cc->EvalNegate(compres);
        auto compres_rot = cc->EvalRotate(compres, -block / 2);

        auto swap_vec_lower = cc->EvalAdd(compres_neg, plain_eliminate);

        auto swap_vec_upper = cc->EvalRotate(swap_vec_lower, -block / 2);

        auto swap_vec_diagonal = cc->EvalAdd(compres, compres_rot);
        end = std::chrono::system_clock::now();
        time = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

        start = std::chrono::system_clock::now();
        res_martrix.push_back(cc->EvalMult(ct_upper, swap_vec_upper));
        res_martrix.push_back(cc->EvalMult(ct, swap_vec_diagonal));
        res_martrix.push_back(cc->EvalMult(ct_lower, swap_vec_lower));
        res = cc->EvalAdd(res_martrix[0], res_martrix[1]);
        res = cc->EvalAdd(res, res_martrix[2]);
        end = std::chrono::system_clock::now();
        time = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    }
}