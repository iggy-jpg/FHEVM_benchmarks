#include "HomOP.h"
#include "HomUtil.h"


#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{

    using namespace TFHEpp;

    static void EQ2_PBS(TLWE<lvl2param> &out,
                        const std::vector<TLWE<lvl2param>> &bits,
                        const EvalKey &ek)
    {
        static const std::vector<int> weight = {1, 4, 1, 4};
        static const std::vector<std::vector<int>> lut = {{1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1}};
        TLWE<lvl2param> linear;
        Linear<lvl2param>(linear, bits, weight, 11);
        PBS_Single(out, linear, ek, lut);
    }

    static const std::vector<int> &eq3_weight()
    { static const std::vector<int> w = {4, 2, 1, 4, 2, 1}; return w; }
    static const std::vector<std::vector<int>> &eq3_lut()
    {
        static const std::vector<std::vector<int>> lut = {
            {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0}};
        return lut;
    }
    static const std::vector<int> &and16_weight()
    { static const std::vector<int> w(16, 1); return w; }
    static const std::vector<std::vector<int>> &and16_lut()
    {   
        static const std::vector<std::vector<int>> lut = {std::vector<int>(16, 0)};
        return lut;
    }

    
    static void EQ3_PBS(TLWE<lvl2param> &out,
                        const std::vector<TLWE<lvl2param>> &a3_msb,
                        const std::vector<TLWE<lvl2param>> &b3_msb,
                        const EvalKey &ek)
    {
        std::vector<TLWE<lvl2param>> in(6);
        in[0] = a3_msb[0];
        in[1] = a3_msb[1];
        in[2] = a3_msb[2];
        HomNOT<lvl2param>(in[3], b3_msb[0]);
        HomNOT<lvl2param>(in[4], b3_msb[1]);
        HomNOT<lvl2param>(in[5], b3_msb[2]);
        TLWE<lvl2param> linear;
        Linear<lvl2param>(linear, in, eq3_weight(), 15);
        PBS_Single(out, linear, ek, eq3_lut());
    }

    static void AND16_PBS(
        TFHEpp::TLWE<TFHEpp::lvl2param> &out,
        const std::vector<TFHEpp::TLWE<TFHEpp::lvl2param>> &bits,
        const TFHEpp::EvalKey &ek)
    {
        using namespace TFHEpp;
        if (bits.size() != 16)
            throw std::runtime_error("AND16_PBS: bits.size() must be 16");
        TLWE<lvl2param> linear;
        Linear<lvl2param>(linear, bits, and16_weight(), 16);
        PBS_Single(out, linear, ek, and16_lut());
    }

    
    static void AND_chain_threads(TLWE<lvl2param> &out,
                                  const std::vector<TLWE<lvl2param>> &bits,
                                  const EvalKey &ek,
                                  int /*threads_hint*/ = -1)
    {
        using TL = TLWE<lvl2param>;
        constexpr int FANIN = 16;
        if (bits.empty())
            throw std::runtime_error("AND_chain_threads: bits.size() must be > 0");

        if (bits.size() <= (size_t)FANIN)
        {
            std::vector<TL> slice(FANIN, CT_ONE);
            for (size_t j = 0; j < bits.size(); ++j)
                slice[j] = bits[j];
            AND16_PBS(out, slice, ek);
            return;
        }

        std::vector<TL> cur(bits.begin(), bits.end()), next;
        while (cur.size() > 1)
        {
            const int ngrp = (int)((cur.size() + FANIN - 1) / FANIN);
            next.assign(ngrp, CT_ONE);
            
            std::vector<TLWE<lvl2param>> linears(ngrp);
#pragma omp parallel for schedule(static) if (ngrp > 1)
            for (int g = 0; g < ngrp; ++g)
            {
                const size_t lo = (size_t)g * FANIN;
                const size_t hi = std::min(cur.size(), lo + FANIN);
                std::vector<TL> slice(FANIN, CT_ONE);
                for (size_t j = lo; j < hi; ++j) slice[j - lo] = cur[j];
                Linear<lvl2param>(linears[g], slice, and16_weight(), 16);
            }
            std::vector<const std::vector<std::vector<int>>*> luts(ngrp, &and16_lut());
            std::vector<std::vector<TLWE<lvl2param>>> outs;
            PBS_Multi_batch(outs, linears, luts, ek, ngrp > 1, std::max(1, ngrp));
            for (int g = 0; g < ngrp; ++g) next[g] = outs[g][0];
            cur.swap(next);
        }
        out = cur[0];
    }

    static void AND_chain(
        TFHEpp::TLWE<TFHEpp::lvl2param> &out,
        const std::vector<TFHEpp::TLWE<TFHEpp::lvl2param>> &bits,
        const TFHEpp::EvalKey &ek)
    {
        using namespace TFHEpp;
        using TL = TLWE<lvl2param>;
        if (bits.empty())
            throw std::runtime_error("AND_chain: bits.size() must be > 0");

        if (bits.size() <= 16)
        {
            std::vector<TL> slice(16, CT_ONE);
            for (size_t j = 0; j < bits.size(); ++j)
                slice[j] = bits[j];
            AND16_PBS(out, slice, ek);
            return;
        }

        std::vector<TL> cur = bits;
        std::vector<TL> next;
        next.reserve((cur.size() + 15) / 16);
        std::vector<TL> slice(16);

        while (cur.size() > 16)
        {
            const size_t csz = cur.size();
            const size_t n_groups = (csz + 15) / 16;
            next.clear();
            next.resize(n_groups);
            for (ptrdiff_t g = 0; g < (ptrdiff_t)n_groups; ++g)
            {
                const size_t base = (size_t)g * 16;
                const size_t upto = std::min(base + 16, csz);
                for (int j = 0; j < 16; ++j)
                    slice[j] = CT_ONE;
                for (size_t idx = base, j = 0; idx < upto; ++idx, ++j)
                    slice[j] = cur[idx];
                AND16_PBS(next[g], slice, ek);
            }
            cur.swap(next);
        }

        for (int j = 0; j < 16; ++j)
            slice[j] = CT_ONE;
        for (size_t j = 0; j < cur.size(); ++j)
            slice[j] = cur[j];
        AND16_PBS(out, slice, ek);
    }

} // namespace

namespace TFHEpp
{

    TLWE<lvl2param> HomEQ(
        const std::vector<TLWE<lvl2param>> &A_in,
        const std::vector<TLWE<lvl2param>> &B_in,
        int W,
        const EvalKey &ek)
    {

    if ((int)A_in.size() != W || (int)B_in.size() != W)
        throw std::runtime_error("HomEQ: input size must equal W");
        TLWE<lvl2param> eq_ct;
        
        const int Wpad = ((W + 2) / 3) * 3;
        const int pad = Wpad - W;
        const int blocks = Wpad / 3;

        std::vector<TLWE<lvl2param>> A(Wpad, CT_ZERO), B(Wpad, CT_ZERO);
        for (int i = 0; i < W; ++i)
        {
            A[pad + i] = A_in[i];
            B[pad + i] = B_in[i];
        }

        std::vector<TLWE<lvl2param>> eq_flags(blocks);
        for (int blk = 0; blk < blocks; ++blk)
        {
            const int base = blk * 3;
            EQ3_PBS(eq_flags[blk], {A[base], A[base + 1], A[base + 2]},
                    {B[base], B[base + 1], B[base + 2]}, ek);
        }
        AND_chain(eq_ct, eq_flags, ek);
        return eq_ct;
    }

    TLWE<lvl2param> HomEQ_pal(const std::vector<TLWE<lvl2param>> &A_in,
                              const std::vector<TLWE<lvl2param>> &B_in,
                              int W,
                              const EvalKey &ek,
                              int /*threads_hint*/)
    {

    if ((int)A_in.size() != W || (int)B_in.size() != W)
        throw std::runtime_error("HomEQ_pal: input size must equal W");
        
        TLWE<lvl2param> eq_ct;
        const int Wpad = ((W + 2) / 3) * 3;
        const int pad = Wpad - W;
        const int blocks = Wpad / 3;

        std::vector<TLWE<lvl2param>> A(Wpad, CT_ZERO), B(Wpad, CT_ZERO);
        for (int i = 0; i < W; ++i)
        {
            A[pad + i] = A_in[i];
            B[pad + i] = B_in[i];
        }

        std::vector<TLWE<lvl2param>> eq_flags(blocks);

        
        {
            std::vector<TLWE<lvl2param>> linears(blocks);
#pragma omp parallel for if (blocks > 1)
            for (int blk = 0; blk < blocks; ++blk)
            {
                const int base = blk * 3;
                std::vector<TLWE<lvl2param>> in(6);
                in[0] = A[base]; in[1] = A[base + 1]; in[2] = A[base + 2];
                HomNOT<lvl2param>(in[3], B[base]);
                HomNOT<lvl2param>(in[4], B[base + 1]);
                HomNOT<lvl2param>(in[5], B[base + 2]);
                Linear<lvl2param>(linears[blk], in, eq3_weight(), 15);
            }
            std::vector<const std::vector<std::vector<int>>*> luts(blocks, &eq3_lut());
            std::vector<std::vector<TLWE<lvl2param>>> outs;
            PBS_Multi_batch(outs, linears, luts, ek, blocks > 1, std::max(1, blocks));
            for (int blk = 0; blk < blocks; ++blk) eq_flags[blk] = outs[blk][0];
        }

        
        AND_chain_threads(eq_ct, eq_flags, ek /*, threads_hint*/);
        return eq_ct;
    }

    TLWE<lvl2param> HomNE(
        const std::vector<TLWE<lvl2param>> &A_in,
        const std::vector<TLWE<lvl2param>> &B_in,
        int W,
        const EvalKey &ek)
    {
        TLWE<lvl2param> eq_ct = HomEQ(A_in, B_in, W, ek);
        TLWE<lvl2param> neq_ct;
        HomNOT<lvl2param>(neq_ct, eq_ct);
        return neq_ct;
    }

    TLWE<lvl2param> HomNE_pal(const std::vector<TLWE<lvl2param>> &A_in,
                              const std::vector<TLWE<lvl2param>> &B_in,
                              int W,
                              const EvalKey &ek,
                              int threads_hint)
    {
        TLWE<lvl2param> eq_ct = HomEQ_pal(A_in, B_in, W, ek, threads_hint);
        TLWE<lvl2param> neq_ct;
        HomNOT<lvl2param>(neq_ct, eq_ct);
        return neq_ct;
    }

} // namespace TFHEpp