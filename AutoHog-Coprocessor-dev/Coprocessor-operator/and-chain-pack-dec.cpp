// evaluate a compound Homgate (PBS logic inlined)

#include "Autohog_Cop.hpp"
#include "timer.h"
#include "HomPack.h"
#include "SecureRandom.h"
#include <mulfft.hpp>
#include <tfhe++.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>
#include <chrono>  
#include "params/128bit.hpp"
#include "drift_ms_tfhepp.hpp"   
#include "drift_ms_params.hpp"   
#include "HomPKE.h"        
#include "HomDecompress.h" // [MOD] lvl1->lvl2 Decompress

using namespace std;
using namespace TFHEpp;

using Clock       = std::chrono::high_resolution_clock;
using MilliDouble = std::chrono::duration<double, std::milli>;

using DriftKeyLvl0 = driftms::DriftMsKey<lvl0param>;
using DriftKeyLvl2 = driftms::DriftMsKey<lvl2param>;

using P   = lvl2param;   
using Pin = lvl1param;   

template <class P2>
std::vector<TLWE<P2>> TRLWE_Decompress_wrapper(const TRLWE<P2>& trlwe, uint32_t k)
{
    return TRLWE_Decompress<P2>(trlwe, k);
}


template <class Q>
inline double tlweSymDecrypt_noise_standard(const TLWE<Q>& c, const SecretKey& sk)
{
    using T = typename Q::T;
    const auto phase_u = tlweSymPhase<Q>(c, sk.key.template get<Q>());
    const auto phase   = static_cast<std::make_signed_t<T>>(phase_u);
    const auto mu_s    = static_cast<std::make_signed_t<T>>(Q::μ);
    const auto m_ctr   = (phase > 0) ? mu_s : -mu_s;
    const auto noise   = static_cast<std::make_signed_t<T>>(phase - m_ctr);
    const int bits     = std::numeric_limits<T>::digits;
    const double scale = std::ldexp(1.0, bits);
    return static_cast<double>(noise) / scale;
}


template <class Q>
inline double tlweSymDecrypt_noise_linear(const TLWE<Q>& c,
                                          const SecretKey& sk,
                                          int /*Kscale*/)
{
    using T = typename Q::T;
    const auto phase_u  = tlweSymPhase<Q>(c, sk.key.template get<Q>());
    const auto phase    = static_cast<std::make_signed_t<T>>(phase_u);
    const auto mu_s     = static_cast<std::make_signed_t<T>>(Q::μ);
    const auto noise    = static_cast<std::make_signed_t<T>>(phase +  mu_s);
    (void)noise;
    const int bits      = std::numeric_limits<T>::digits;
    const double scale  = std::ldexp(1.0, bits);
    return static_cast<double>(phase_u - mu_s) / scale - 0.5;
}


template <class Q>
int tlweSymPhase_rho_mod(const TLWE<Q>& c, const Key<Q>& key)
{
    using T = typename Q::T;
    using BaseP = lvl2param;
    constexpr int N   = static_cast<int>(BaseP::n);
    constexpr int MOD = 2 * N;

    const int bits = std::numeric_limits<T>::digits;
    const unsigned __int128 DEN  = (unsigned __int128)1u << bits;
    const unsigned __int128 HALF = DEN >> 1;

    const unsigned __int128 b_u = (unsigned __int128)c[Q::k * Q::n];
    const unsigned long long floor_b =
        (unsigned long long)((b_u * (unsigned __int128)MOD) / DEN);

    unsigned __int128 sum = 0;
    for (int k = 0; k < Q::k; ++k) {
        const int base = k * Q::n;
        for (int i = 0; i < Q::n; ++i) {
            const unsigned __int128 a_u = (unsigned __int128)c[base + i];
            const unsigned long long round_ai =
                (unsigned long long)((a_u * (unsigned __int128)MOD + HALF) / DEN);
            const unsigned long long sbit = (unsigned long long)(key[base + i]) & 1ULL;
            sum += (unsigned __int128)round_ai * (unsigned __int128)sbit;
        }
    }

    const unsigned long long sum_mod = (unsigned long long)(sum % (unsigned __int128)MOD);
    int rho = (int)((floor_b % MOD + MOD) % MOD);
    rho -= (int)sum_mod;
    rho %= MOD;
    if (rho < 0) rho += MOD;
    return rho;
}

template <class Q>
int tlweSymPhase_rho_mod(const TLWE<Q>& c, const SecretKey& sk)
{
    return tlweSymPhase_rho_mod<Q>(c, sk.key.template get<Q>());
}

int main()
{
    autohog::SecureRng engine;
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    const int input_num  = 16;
    const int output_num = 2;
    const int REPS       = 1000;   

    std::vector<uint8_t> plain(input_num);
    std::vector<uint8_t> pres(3);

    for (int i = 0; i < input_num; i++) plain[i] = static_cast<uint8_t>(binary(engine) > 0);
    for (int i = 0; i < input_num; i++) plain[i] = 1;

    std::vector<std::vector<int>> lut = {
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    };
    std::vector<int> weight(input_num, 1);

    auto* sk = new SecretKey();

    EvalKey ek;
    ek.emplacebkfft<lvl02param>(*sk);
    ek.emplaceiksk<lvl20param>(*sk);
    ek.emplaceiksk<lvl12param>(*sk);                                 
    PKE::SeededRlwePkPair<P> spk;
    spk = PKE::generate_rlwe_pk_pair_seeded<P>(*sk, engine);

    std::array<bool, P::n> bits_array{};
    for (size_t i = 0; i < plain.size() && i < P::n; ++i)
        bits_array[i] = (plain[i] != 0);

    
    DriftKeyLvl0 drift_key_lvl0;
    {
        using TestParam = DriftMsTestParamLvl0;
        using P0        = lvl0param;
        drift_key_lvl0.r_sigma        = TestParam::r_sigma;
        drift_key_lvl0.bound          = TestParam::bound();
        drift_key_lvl0.input_variance = TestParam::input_variance();
        drift_key_lvl0.log_mod        = TestParam::log_mod;
        drift_key_lvl0.zeros.resize(TestParam::zeros_count);

        for (std::size_t i = 0; i < drift_key_lvl0.zeros.size(); ++i) {
            drift_key_lvl0.zeros[i] = tlweSymEncrypt<P0>(
                0,              
                P0::α,          
                sk->key.get<P0>()
            );
        }
    }
    
    DriftKeyLvl2 drift_key_lvl2;
    {
        using TestParam = DriftMsTestParamLvl2; 
        using P2        = lvl2param;

        drift_key_lvl2.r_sigma        = TestParam::r_sigma;
        drift_key_lvl2.bound          = TestParam::bound();          
        drift_key_lvl2.input_variance = TestParam::input_variance(); 
        drift_key_lvl2.log_mod        = TestParam::log_mod;
        drift_key_lvl2.zeros.resize(TestParam::zeros_count);

        for (std::size_t i = 0; i < drift_key_lvl2.zeros.size(); ++i) {
            drift_key_lvl2.zeros[i] = tlweSymEncrypt<P2>(
                0,              
                P2::α,          
                sk->key.get<P2>()
            );
        }
    }

    // ================== [DRIFT END] ==================
    auto  ksk12_sp = ek.get<KeySwitchingKey<lvl12param>>();  // [MOD]
    auto &ksk12    = *ksk12_sp;                              // [MOD]

    
    AnnihilateKey<Pin> ahk_l1;
    AnnihilateKey<P>   ahk_l2;
    annihilatekeygen<Pin>(ahk_l1, *sk);
    annihilatekeygen<P>(ahk_l2,   *sk);

    
    auto pk_lwe = PKE::generate_public_key<Pin>(*sk);

    std::vector<uint8_t>   m_bits(input_num);   
    std::vector<TLWE<Pin>> lwe_vec(input_num);  

    
    double sum_r1_ms = 0.0;  // PBS single
    double sum_r2_ms = 0.0;  // PBS multi
    double sum_r3_ms = 0.0;  // pack+unpack+PBS single

    
    for (int rep = 0; rep < REPS; ++rep)
    {
        // --- fresh encrypt on lvl1 (LWE, public-key) + packing ---
        
        for (int i = 0; i < input_num; ++i) {
            m_bits[i]  = plain[i] & 1u;
            lwe_vec[i] = PKE::encrypt_public_key_bit<Pin>(pk_lwe, m_bits[i], engine);
        }

        
        TRLWE<Pin> packed_lvl1;
        HomPack::TLWE2TRLWEPacking<Pin>(packed_lvl1, lwe_vec, ahk_l1);

        
        TRLWE<P> ct_trlwe;  
        HomDecompress::Decompress_DirectIKS_lvl1_to_lvl2(
            ct_trlwe,            // out: lvl2 TRLWE
            packed_lvl1,         // in : lvl1 TRLWE
            ksk12,               // IKS: lvl1 -> lvl2
            ahk_l2,              
            input_num,          
            &drift_key_lvl2);          

        std::vector<TLWE<P>> cipher_in(input_num), cipher_out(output_num);
        TLWE<P> linear;

        cipher_in = TRLWE_Decompress_wrapper<P>(ct_trlwe, input_num);
        double nz_decomp0 = tlweSymDecrypt_noise_standard<P>(cipher_in[0], *sk);

        
        const int t_size = static_cast<int>(lut[0].size());
        const int sum_w  = std::accumulate(weight.begin(), weight.end(), 0);
        const int Kscale = 2 * (16 - t_size) + sum_w + 1;

        // ======== [INLINE PBS SINGLE / round1] ========
        double nz_iks0, nz0;
        int    rho_ks0;
        auto t_r1_start = Clock::now();
        {
            TLWE<lvl0param> tlwe0;
            auto &ek_nc = ek;
            auto  ksk_sp = ek_nc.get<KeySwitchingKey<lvl20param>>();
            auto &ksk    = *ksk_sp;

            Linear<P>(linear, cipher_in, weight, t_size);
            IdentityKeySwitch<lvl20param>(tlwe0, linear, ksk);
            driftms::improve_lwe_ciphertext_modulus_switch_noise_for_binary_key<lvl0param>(tlwe0, drift_key_lvl0);
            //centered::apply_centered_correction_to_body<lvl0param>(tlwe0, drift_key_lvl0.log_mod);
            

            rho_ks0 = tlweSymPhase_rho_mod<lvl0param>(tlwe0, *sk);
            nz_iks0 = tlweSymDecrypt_noise_linear<lvl0param>(tlwe0, *sk, Kscale);

            Polynomial<lvl2param> poly;
            for (auto &p : poly) p = 0;
            constexpr int padding = 128;
            for (int j = 0; j < padding; j++) poly[j] = lvl2param::μ;

            std::array<typename lvl2param::T, 16> t{};
            for (int i = 0; i < t_size; i++) t[(i + 16 - t_size)] = int(lut[0][i]);
            for (int i = 0; i < 16; i++)
                for (int j = 0; j < padding; j++)
                    poly[i * padding + j] = (t[i] * 2 - 1) * lvl2param::μ;

            TRLWE<lvl02param::targetP> acc;
            auto  bkfft_sp = ek_nc.get<BootstrappingKeyFFT<lvl02param>>();
            auto &bkfft    = *bkfft_sp;
            BlindRotate<lvl02param>(acc, tlwe0, bkfft, poly);

            SampleExtractIndex<lvl02param::targetP>(cipher_out[0], acc, 0);
        }
        auto t_r1_end = Clock::now();
        sum_r1_ms += MilliDouble(t_r1_end - t_r1_start).count();

        pres[0] = tlweSymDecrypt<P>(cipher_out[0], *sk);
        nz0 = tlweSymDecrypt_noise_standard<P>(cipher_out[0], *sk);

        // ======== [INLINE PBS MULTI / round2] ========
        double nz_iks1, nz1;
        int    rho_ks1;
        auto t_r2_start = Clock::now();
        {
            std::fill(cipher_in.begin(), cipher_in.end(), cipher_out[0]);

            TLWE<lvl0param> tlwe0;
            auto &ek_nc = ek;
            auto  ksk_sp = ek_nc.get<KeySwitchingKey<lvl20param>>();
            auto &ksk    = *ksk_sp;

            Linear<P>(linear, cipher_in, weight, t_size);
            IdentityKeySwitch<lvl20param>(tlwe0, linear, ksk);

            
            driftms::improve_lwe_ciphertext_modulus_switch_noise_for_binary_key<lvl0param>(tlwe0, drift_key_lvl0);
            //centered::apply_centered_correction_to_body<lvl0param>(tlwe0, drift_key_lvl0.log_mod);

            rho_ks1 = tlweSymPhase_rho_mod<lvl0param>(tlwe0, *sk);
            nz_iks1 = tlweSymDecrypt_noise_linear<lvl0param>(tlwe0, *sk, Kscale);

            Polynomial<lvl2param> poly_base;
            for (auto &p : poly_base) p = 0;
            constexpr int padding = 128;
            for (int j = 0; j < padding; j++) poly_base[j] = lvl2param::μ;

            TRLWE<lvl02param::targetP> acc, cres;
            auto  bkfft_sp = ek_nc.get<BootstrappingKeyFFT<lvl02param>>();
            auto &bkfft    = *bkfft_sp;
            BlindRotate<lvl02param>(acc, tlwe0, bkfft, poly_base);

            const int output_num_local = (int)lut.size();
            cipher_out.resize(output_num_local);
            for (int k = 0; k < output_num_local; k++)
            {
                std::array<typename lvl2param::T, lvl2param::n> t{};
                const int t_size_k = (int)lut[k].size();
                for (int i = 0; i < t_size_k; i++)
                    t[(i + 16 - t_size_k) * padding] = (lut[k][i]) * 2 - 1;

                for (int j = 0; j < lvl2param::k + 1; j++)
                    PolyMul<lvl2param>(cres[j], acc[j], t);

                SampleExtractIndex<lvl02param::targetP>(cipher_out[k], cres, 0);
            }
        }
        auto t_r2_end = Clock::now();
        sum_r2_ms += MilliDouble(t_r2_end - t_r2_start).count();

        pres[1] = tlweSymDecrypt<P>(cipher_out[0], *sk);
        nz1 = tlweSymDecrypt_noise_standard<P>(cipher_out[0], *sk);

        // ======== [PACK -> UNPACK -> PBS SINGLE / round3] ========
        double nz_unpacked0, nz_iks2, nz2;
        int    rho_ks2;
        auto t_r3_start = Clock::now();
        {
            constexpr uint32_t K = 256;
            static_assert(K <= P::n, "K must be <= poly degree");
            std::vector<TLWE<P>> clones(K, cipher_out[0]);
            AnnihilateKey<P> ahk;
            annihilatekeygen<P>(ahk, *sk);
            TRLWE<P> packed_trlwe;
            HomPack::TLWE2TRLWEPacking<P>(packed_trlwe, clones, ahk);
            std::vector<TLWE<P>> unpacked = TRLWE_Decompress_wrapper<P>(packed_trlwe, K);

            nz_unpacked0 = tlweSymDecrypt_noise_standard<P>(unpacked[0], *sk);

            std::fill(cipher_in.begin(), cipher_in.end(), unpacked[0]);

            TLWE<lvl0param> tlwe0;
            auto &ek_nc = ek;
            auto  ksk_sp = ek_nc.get<KeySwitchingKey<lvl20param>>();
            auto &ksk    = *ksk_sp;

            Linear<P>(linear, cipher_in, weight, t_size);
            IdentityKeySwitch<lvl20param>(tlwe0, linear, ksk);

            
            driftms::improve_lwe_ciphertext_modulus_switch_noise_for_binary_key<lvl0param>(tlwe0, drift_key_lvl0);
            //centered::apply_centered_correction_to_body<lvl0param>(tlwe0, drift_key_lvl0.log_mod);

            rho_ks2 = tlweSymPhase_rho_mod<lvl0param>(tlwe0, *sk);
            nz_iks2 = tlweSymDecrypt_noise_linear<lvl0param>(tlwe0, *sk, Kscale);

            Polynomial<lvl2param> poly;
            for (auto &p : poly) p = 0;
            constexpr int padding = 128;
            for (int j = 0; j < padding; j++) poly[j] = lvl2param::μ;

            std::array<typename lvl2param::T, 16> t{};
            for (int i = 0; i < t_size; i++) t[(i + 16 - t_size)] = int(lut[0][i]);
            for (int i = 0; i < 16; i++)
                for (int j = 0; j < padding; j++)
                    poly[i * padding + j] = (t[i] * 2 - 1) * lvl2param::μ;

            TRLWE<lvl02param::targetP> acc;
            auto  bkfft_sp = ek_nc.get<BootstrappingKeyFFT<lvl02param>>();
            auto &bkfft    = *bkfft_sp;
            BlindRotate<lvl02param>(acc, tlwe0, bkfft, poly);

            SampleExtractIndex<lvl02param::targetP>(cipher_out[0], acc, 0);
        }
        auto t_r3_end = Clock::now();
        sum_r3_ms += MilliDouble(t_r3_end - t_r3_start).count();

        pres[2] = tlweSymDecrypt<P>(cipher_out[0], *sk);
        nz2 = tlweSymDecrypt_noise_standard<P>(cipher_out[0], *sk);

        
        uint8_t gt_and16 = 1;
        for (int i = 0; i < input_num; ++i) gt_and16 &= (plain[i] & 1u);
        uint8_t gt_chain = gt_and16;

        if (pres[2] != gt_chain) {
            std::printf("Plain AND16 = %u, Chain(AND self x16) = %u, Decrypted pres[0] = %u, pres[1] = %u, pres[2] = %u\n",
                        (unsigned)gt_and16, (unsigned)gt_chain,
                        (unsigned)pres[0], (unsigned)pres[1], (unsigned)pres[2]);
            puts("Result: WRONG ❌");
        }

        
        std::printf("%.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f\n",
            nz_decomp0,
            nz_iks0,  (2048.0 + 64.0 - rho_ks0) / 4096.0,  nz0,
            nz_iks1,  (2048.0 + 64.0 - rho_ks1) / 4096.0,  nz1,
            nz_unpacked0,
            nz_iks2,  (2048.0 + 64.0 - rho_ks2) / 4096.0,  nz2);
    }

    
    std::printf(
        "# AVG time over %d reps (ms):\n"
        "#   round1_PBS_single      = %.6f ms\n"
        "#   round2_PBS_multi       = %.6f ms\n"
        "#   round3_pack_PBS_single = %.6f ms\n",
        REPS,
        sum_r1_ms / REPS,
        sum_r2_ms / REPS,
        sum_r3_ms / REPS);

    delete sk;
    return 0;
}
