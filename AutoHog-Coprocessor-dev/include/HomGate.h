#pragma once
#include <mulfft.hpp>
#include <tfhe++.hpp>

#include "drift_ms_tfhepp.hpp"
#include "HomTrace.h"

using DriftKeyLvl0 = driftms::DriftMsKey<TFHEpp::lvl0param>;

namespace TFHEpp
{
    template <class P>
    void Linear(TLWE<P> &cres,
                const std::vector<TLWE<P>> &cipher,
                const std::vector<int> &weight,
                const int t_size)
    {
        if (hom_trace_schedule_only()) {
            int level = 0;
            for (const auto &input : cipher)
                level = std::max(level,
                    hom_trace_decode_level(static_cast<std::uint64_t>(input[0])));
            cres = {};
            cres[0] = static_cast<typename P::T>(hom_trace_encode_level(level));
            return;
        }
        const int input_num = static_cast<int>(weight.size());
        for (int j = 0; j <= P::k * P::n; j++)
            cres[j] = 0;
        int sum = 2 * (16 - t_size);
        for (int i = 0; i < input_num; i++)
            sum += weight[i];
        for (int i = 0; i < input_num; i++)
            for (int j = 0; j <= P::k * P::n; j++)
                cres[j] += static_cast<typename P::T>(
                    static_cast<long long>(cipher[i][j]) * static_cast<long long>(weight[i]));
        cres[P::k * P::n] += static_cast<typename P::T>(
            static_cast<long long>(sum + 1) * static_cast<long long>(P::μ));
    }

    // ============================ PBS_Single ============================

    inline void PBS_Single(TLWE<lvl2param> &res,
                       const TLWE<lvl2param> &tlwe,
                       const EvalKey &ek,
                       const std::vector<std::vector<int>> &lut,
                        const DriftKeyLvl0 *drift_key_lvl0 = nullptr)
    {
        if (hom_trace_schedule_only()) {
            res = {};
            const int level = hom_trace_decode_level(
                                  static_cast<std::uint64_t>(tlwe[0])) + 1;
            hom_trace_scheduled_at(level, "PBS_Single");
            res[0] = static_cast<lvl2param::T>(hom_trace_encode_level(level));
            return;
        }
        hom_trace_scheduled("PBS_Single");
        const int t_size = std::size(lut[0]);

        TLWE<lvl0param> tlwe0;
        auto &ek_nc = const_cast<EvalKey &>(ek);

        // 1) lvl2 -> lvl0：KeySwitch
        {
            auto ksk_sp = ek_nc.get<KeySwitchingKey<lvl20param>>();
            auto &ksk   = *ksk_sp;
            IdentityKeySwitch<lvl20param>(tlwe0, tlwe, ksk);
            ksk_sp.reset();
        }

        
        if (drift_key_lvl0 != nullptr) {
            driftms::improve_lwe_ciphertext_modulus_switch_noise_for_binary_key<lvl0param>(
                tlwe0, *drift_key_lvl0);
        }

        
        {
            
            Polynomial<lvl2param> poly;
            constexpr int padding = 128;

            std::array<typename lvl2param::T, 16> t{};
            for (int i = 0; i < t_size; i++) t[(i + 16 - t_size)] = int(lut[0][i]);
            for (int i = 0; i < 16; i++)
                for (int j = 0; j < padding; j++)
                    poly[i * padding + j] = (t[i] * 2 - 1) * lvl2param::μ;

            TRLWE<lvl02param::targetP> acc;
            {
                auto bkfft_sp = ek_nc.get<BootstrappingKeyFFT<lvl02param>>();
                auto &bkfft   = *bkfft_sp;
                BlindRotate<lvl02param>(acc, tlwe0, bkfft, poly);
                hom_trace_actual("PBS_Single");
                bkfft_sp.reset();
            }

            SampleExtractIndex<lvl02param::targetP>(res, acc, 0);
        }
    }

    // ============================ PBS_Multi ============================

    
    
    inline void PBS_Multi_body(std::vector<TLWE<lvl2param>> &res,
                               const TLWE<lvl0param>        &tlwe0,
                               const EvalKey                &ek,
                               const std::vector<std::vector<int>> &lut);

    
    
    
    
    //
    
    
    
    
    template <int CAT>
    inline int BatchKS20_chunks(std::vector<TLWE<lvl0param>> &outs,
                                const TLWE<lvl2param> *ins,
                                const KeySwitchingKey<lvl20param> &ksk,
                                int begin, int n, int num_threads_hint)
    {
        const int nchunk = (n - begin) / CAT;
        if (nchunk <= 0) return begin;
        const int useT = std::max(1, std::min(nchunk, num_threads_hint));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (useT > 1) num_threads(useT)
#endif
        for (int ck = 0; ck < nchunk; ++ck) {
            const int base = begin + ck * CAT;
            
            
            static thread_local std::array<TLWE<lvl2param>, CAT> in;
            static thread_local std::array<TLWE<lvl0param>, CAT> out;
            for (int c = 0; c < CAT; ++c) in[c] = ins[base + c];
            CatIdentityKeySwitch<lvl20param, CAT>(out, in, ksk);
            for (int c = 0; c < CAT; ++c) outs[base + c] = out[c];
        }
        return begin + nchunk * CAT;
    }

    inline void BatchKeySwitch20(std::vector<TLWE<lvl0param>> &outs,
                                 const std::vector<TLWE<lvl2param>> &ins_all,
                                 const EvalKey &ek,
                                 int num_threads_hint = 1,
                                 int off = 0,
                                 int len = -1)
    {
        const int n = (len < 0) ? (static_cast<int>(ins_all.size()) - off) : len;
        const TLWE<lvl2param> *ins = ins_all.data() + off;
        outs.resize(n);
        auto &ek_nc = const_cast<EvalKey &>(ek);
        auto ksk_sp = ek_nc.get<KeySwitchingKey<lvl20param>>();
        auto &ksk   = *ksk_sp;

        int done = 0;
        done = BatchKS20_chunks<16>(outs, ins, ksk, done, n, num_threads_hint);
        done = BatchKS20_chunks<8>(outs, ins, ksk, done, n, num_threads_hint);
        done = BatchKS20_chunks<4>(outs, ins, ksk, done, n, num_threads_hint);
        done = BatchKS20_chunks<2>(outs, ins, ksk, done, n, num_threads_hint);
        for (int i = done; i < n; ++i)
            IdentityKeySwitch<lvl20param>(outs[i], ins[i], ksk);
        ksk_sp.reset();
    }

    
    inline void PBS_Multi_batch(
        std::vector<std::vector<TLWE<lvl2param>>> &res,
        const std::vector<TLWE<lvl2param>> &linears,
        const std::vector<const std::vector<std::vector<int>>*> &luts,
        const EvalKey &ek,
        bool parallel,
        int num_threads_hint)
    {
        const int n = static_cast<int>(linears.size());
        res.resize(n);
        if (hom_trace_schedule_only()) {
            for (int i = 0; i < n; ++i) {
                const int level = hom_trace_decode_level(
                                      static_cast<std::uint64_t>(linears[i][0])) + 1;
                hom_trace_scheduled_at(level, "PBS_Multi_batch");
                res[i].resize(luts[i] != nullptr ? luts[i]->size() : 0);
                for (auto &output : res[i]) {
                    output = {};
                    output[0] = static_cast<lvl2param::T>(hom_trace_encode_level(level));
                }
            }
            return;
        }
        hom_trace_scheduled("PBS_Multi_batch", static_cast<std::uint64_t>(n));
        const int useT = parallel ? std::max(1, num_threads_hint) : 1;

        
        
        
        const int slice = std::max(64, std::min(512, useT * 16));
        std::vector<TLWE<lvl0param>> tlwe0s;
        for (int base = 0; base < n; base += slice) {
            const int len = std::min(slice, n - base);
            BatchKeySwitch20(tlwe0s, linears, ek, useT, base, len);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (parallel) num_threads(useT)
#endif
            for (int t = 0; t < len; ++t) {
                PBS_Multi_body(res[base + t], tlwe0s[t], ek, *luts[base + t]);
            }
        }
        hom_trace_actual("PBS_Multi_batch", static_cast<std::uint64_t>(n));
    }

    inline void PBS_Multi(std::vector<TLWE<lvl2param>> &res,
                      const TLWE<lvl2param>        &tlwe,
                      const EvalKey                &ek,
                      const std::vector<std::vector<int>> &lut,
                      const DriftKeyLvl0          *drift_key_lvl0 = nullptr)
    {
        if (hom_trace_schedule_only()) {
            res.resize(lut.size());
            const int level = hom_trace_decode_level(
                                  static_cast<std::uint64_t>(tlwe[0])) + 1;
            hom_trace_scheduled_at(level, "PBS_Multi");
            for (auto &output : res) {
                output = {};
                output[0] = static_cast<lvl2param::T>(hom_trace_encode_level(level));
            }
            return;
        }
        hom_trace_scheduled("PBS_Multi");
        TLWE<lvl0param> tlwe0;
        auto &ek_nc = const_cast<EvalKey &>(ek);

        // 1) lvl2 -> lvl0: KeySwitch
        {
            auto ksk_sp = ek_nc.get<KeySwitchingKey<lvl20param>>();
            auto &ksk   = *ksk_sp;
            IdentityKeySwitch<lvl20param>(tlwe0, tlwe, ksk);
            ksk_sp.reset();
        }

        
        if (drift_key_lvl0 != nullptr) {
            driftms::improve_lwe_ciphertext_modulus_switch_noise_for_binary_key<lvl0param>(
                tlwe0, *drift_key_lvl0);
        }

        PBS_Multi_body(res, tlwe0, ek, lut);
        hom_trace_actual("PBS_Multi");
    }

    inline void PBS_Multi_body(std::vector<TLWE<lvl2param>> &res,
                               const TLWE<lvl0param>        &tlwe0,
                               const EvalKey                &ek,
                               const std::vector<std::vector<int>> &lut)
    {
        const int output_num = std::size(lut);
        const int t_size     = std::size(lut[0]);
        auto &ek_nc = const_cast<EvalKey &>(ek);

        res.resize(output_num);

        
        {
            Polynomial<lvl2param> poly;
            for (auto &p : poly) p = 0;
            constexpr int padding = 128;
            for (int j = 0; j < padding; j++)
                poly[j] = lvl2param::μ;

            TRLWE<lvl02param::targetP> acc;

            {
                auto bkfft_sp = ek_nc.get<BootstrappingKeyFFT<lvl02param>>();
                auto &bkfft   = *bkfft_sp;
                BlindRotate<lvl02param>(acc, tlwe0, bkfft, poly);
                bkfft_sp.reset();
            }

            
            
            
            
            {
                using T = typename lvl2param::T;
                constexpr int N = lvl2param::n;
                alignas(64) static thread_local std::array<Polynomial<lvl2param>,
                                                           lvl2param::k> spa;

                for (int k = 0; k < output_num; k++) {
                    T b0 = 0;
                    bool first = true;
                    for (int i = 0; i < t_size; i++) {
                        const int p = (i + 16 - t_size) * padding;
                        const bool pos = lut[k][i] != 0;  

                        
                        for (int r = 0; r < lvl2param::k; r++) {
                            const auto &a = acc[r];
                            auto &d = spa[r];
                            if (first) {
                                if (pos) {
                                    for (int j = 0; j < p; j++) d[j] = -a[N - p + j];
                                    for (int j = p; j < N; j++) d[j] = a[j - p];
                                } else {
                                    for (int j = 0; j < p; j++) d[j] = a[N - p + j];
                                    for (int j = p; j < N; j++) d[j] = -a[j - p];
                                }
                            } else {
                                if (pos) {
                                    for (int j = 0; j < p; j++) d[j] -= a[N - p + j];
                                    for (int j = p; j < N; j++) d[j] += a[j - p];
                                } else {
                                    for (int j = 0; j < p; j++) d[j] += a[N - p + j];
                                    for (int j = p; j < N; j++) d[j] -= a[j - p];
                                }
                            }
                        }

                        
                        const auto &b = acc[lvl2param::k];
                        const T contrib = (p == 0) ? b[0] : static_cast<T>(-b[N - p]);
                        b0 += pos ? contrib : static_cast<T>(-contrib);
                        first = false;
                    }

                    
                    for (int r = 0; r < lvl2param::k; r++) {
                        res[k][r * N] = spa[r][0];
                        for (int i = 1; i < N; i++)
                            res[k][r * N + i] = -spa[r][N - i];
                    }
                    res[k][lvl2param::k * N] = b0;
                }
            }
            
        }
    }

    template <class P, int casign, int cbsign, typename P::T offset>
    inline void HomGate(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
                        const EvalKey &ek)
    {
        for (int i = 0; i <= P::k * P::n; i++)
            res[i] = 16 * (casign * ca[i] + cbsign * cb[i]);
        res[P::k * P::n] += 16 * offset;
        GateBootstrapping(res, res, ek);
    }
}
