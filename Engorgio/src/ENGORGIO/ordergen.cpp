#include "utils.h"
#include "comp.h"
#include <chrono>
#include "math/chebyshev.h"
#include <stdexcept>
namespace openfhe
{

    void bitonic_comp_unbounded(int stage, int part, int slots, int length, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ct,
                                Ciphertext<lbcrypto::DCRTPoly> &res, double precision, std::vector<double> &coefficients,
                                lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        auto cc = ct[0]->GetCryptoContext();
        int num_ct = length / slots;
        int Block = 1 << (stage + 1);
        int block = 1 << (stage + 1 - part);
        auto ct_mask = cc->EvalNegate(ct[0]);
        auto ct_rot = cc->EvalNegate(ct[1]);
        comp_partial(ct_mask, ct_rot, precision, coefficients, res, privateKey);
    }

    void bitonic_comp_unbounded_modular(int stage, int part, int slots, int length, std::vector<std::vector<Ciphertext<lbcrypto::DCRTPoly>>> &ct,
                                        Ciphertext<lbcrypto::DCRTPoly> &res, double precision, std::vector<double> &coefficients,
                                        lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey, int blocks)
    {
        auto cc = ct[0][0]->GetCryptoContext();
        int num_ct = length / slots;
        int Block = 1 << (stage + 1);
        int block = 1 << (stage + 1 - part);
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ct_mask(num_ct), ct_rot(num_ct);
        for (int i = 0; i < num_ct; i++)
        {
            ct_mask[i] = cc->EvalNegate(ct[0][i]);
            ct_rot[i] = cc->EvalNegate(ct[1][i]);
        }
        comp_partial_modular(ct_mask, ct_rot, precision, coefficients, res, privateKey);
    }

    void bitonic_comp_topk(int stage, int slots, int topk_i, int k, Ciphertext<lbcrypto::DCRTPoly> &ct,
                           Ciphertext<lbcrypto::DCRTPoly> &res, double precision, std::vector<double> &coefficients,
                           lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        auto cc = ct->GetCryptoContext();
        int Block = (1 << int(std::log2(k) + 1)) + (pow(2, topk_i + 1) - 1) * 2 * k; // 16
        int block = 1 << (int(std::log2(k)) + 1);                                    // 8
        int step = (1 << int(std::log2(k))) + (pow(2, topk_i + 1) - 1) * 2 * k;      // 12
        std::vector<double> mask(slots, 1), eliminate(slots, 0);
        // if (block < slots)

        for (int i = 0; i < slots; i += Block)
        {
            for (int j = 0; j < block / 2; j++)
            {
                eliminate[i + j] = 1;
            }
        }
        for (int i = 0; i < mask.size(); i++)
        {
            mask[i] = mask[i] * eliminate[i];
        }
        for (int i = 0; i < slots; i += 2 * Block)
        {
            for (int j = 0; j < block / 2; j++)
            {
                mask[i + j] = -1;
                mask[i + j + step] = -1;
                mask[i + Block + j] = 1;
                mask[i + Block + j + step] = 1;
            }
        }
        Plaintext plain_mask = cc->MakeCKKSPackedPlaintext(mask);
        Plaintext plain_eliminate = cc->MakeCKKSPackedPlaintext(eliminate);
        auto ct_mask = cc->EvalMult(ct, plain_mask);
        auto ct_rot = cc->EvalRotate(ct_mask, step);

        ct_mask = cc->EvalMult(ct_mask, plain_eliminate);

        ct_rot = cc->EvalMult(ct_rot, plain_eliminate);

        Ciphertext<lbcrypto::DCRTPoly> ploy_res;
        comp_partial(ct_mask, ct_rot, precision, coefficients, res, privateKey);
    }
    void bitonic_comp(int stage, int part, int slots, Ciphertext<lbcrypto::DCRTPoly> &ct,
                      Ciphertext<lbcrypto::DCRTPoly> &res, double precision, std::vector<double> &coefficients,
                      lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        auto cc = ct->GetCryptoContext();
        int Block = 1 << (stage + 1);
        int block = 1 << (stage + 1 - part);
        std::vector<double> mask(slots, 1), eliminate(slots, 0);
        // if (block < slots)
        for (int i = 0; i < slots; i += 2 * Block)
        {
            for (int j = 0; j < Block; j++)
            {
                mask[i + j] = -1;
            }
        }
        // else

        for (int i = 0; i < slots; i += block)
        {
            for (int j = 0; j < block / 2; j++)
            {
                eliminate[i + j] = 1;
            }
        }
        Plaintext plain_mask = cc->MakeCKKSPackedPlaintext(mask);
        Plaintext plain_eliminate = cc->MakeCKKSPackedPlaintext(eliminate);
        auto ct_mask = cc->EvalMult(ct, plain_mask);
        // std::cout << "\n EvalRotate: " << block / 2 << std::endl;
        auto ct_rot = cc->EvalRotate(ct_mask, block / 2);

        // std::cout << "\n comp_partial " << std::endl;
        ct_mask = cc->EvalMult(ct_mask, plain_eliminate);

        ct_rot = cc->EvalMult(ct_rot, plain_eliminate);
        Ciphertext<lbcrypto::DCRTPoly> ploy_res;
        comp_partial(ct_mask, ct_rot, precision, coefficients, res, privateKey);
    }

    void bitonic_comp_modular(int stage, int part, int slots, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ct,
                              Ciphertext<lbcrypto::DCRTPoly> &res, double precision, std::vector<double> &coefficients,
                              lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        auto cc = ct[0]->GetCryptoContext();
        int Block = 1 << (stage + 1);
        int block = 1 << (stage + 1 - part);
        int num_ct = ct.size();
        std::vector<double> mask(slots, 1), eliminate(slots, 0);
        for (int i = 0; i < slots; i += 2 * Block)
        {
            for (int j = 0; j < Block; j++)
            {
                mask[i + j] = -1;
            }
        }
        // else

        for (int i = 0; i < slots; i += block)
        {
            for (int j = 0; j < block / 2; j++)
            {
                eliminate[i + j] = 1;
            }
        }
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ct_mask(num_ct), ct_rot(num_ct);
        Plaintext plain_mask = cc->MakeCKKSPackedPlaintext(mask);
        Plaintext plain_eliminate = cc->MakeCKKSPackedPlaintext(eliminate);
        for (int i = 0; i < num_ct; i++)
        {
            ct_mask[i] = cc->EvalMult(ct[i], plain_mask);
            ct_rot[i] = cc->EvalRotate(ct_mask[i], block / 2);
            ct_mask[i] = cc->EvalMult(ct_mask[i], plain_eliminate);
            ct_rot[i] = cc->EvalMult(ct_rot[i], plain_eliminate);
        }
        Ciphertext<lbcrypto::DCRTPoly> ploy_res;
        comp_partial_modular(ct_mask, ct_rot, precision, coefficients, ploy_res, privateKey);
        res = cc->EvalMult(ploy_res, plain_eliminate);
    }
}