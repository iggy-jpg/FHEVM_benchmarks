#include "HomOP.h"
#include "HomUtil.h"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using namespace TFHEpp;



static void SuffixOr2_PBS(
    std::vector<TLWE<lvl2param>> &out,
    const TLWE<lvl2param> &s,
    const TLWE<lvl2param> &bhi,
    const TLWE<lvl2param> &blo,
    const EvalKey &ek)
{
    static const std::vector<int> weight = {1, 2, 4};
    static const std::vector<std::vector<int>> lut = [] {
        std::vector<int> hi(8), lo(8);
        for (int idx = 0; idx < 8; ++idx) {
            const int sv = idx & 1, bh = (idx >> 1) & 1, bl = (idx >> 2) & 1;
            hi[idx] = sv | bh;
            lo[idx] = sv | bh | bl;
        }
        return std::vector<std::vector<int>>{hi, lo};
    }();
    TLWE<lvl2param> linear;
    Linear<lvl2param>(linear, {s, bhi, blo}, weight, 8);
    PBS_Multi(out, linear, ek, lut);
}




static void QBit3_PBS(
    std::vector<TLWE<lvl2param>> &out,
    const TLWE<lvl2param> &borrow,
    const TLWE<lvl2param> &sfx_w,
    const TLWE<lvl2param> &sfx_0,
    const EvalKey &ek)
{
    static const std::vector<int> weight = {1, 2, 4};
    static const std::vector<std::vector<int>> lut = [] {
        std::vector<int> ge(8), q(8);
        for (int idx = 0; idx < 8; ++idx) {
            const int bo = idx & 1, sw = (idx >> 1) & 1, s0 = (idx >> 2) & 1;
            ge[idx] = (!bo && !sw) ? 1 : 0;
            q[idx] = (!bo && !sw && s0) ? 1 : 0;
        }
        return std::vector<std::vector<int>>{ge, q};
    }();
    TLWE<lvl2param> linear;
    Linear<lvl2param>(linear, {borrow, sfx_w, sfx_0}, weight, 8);
    PBS_Multi(out, linear, ek, lut);
}


static void QBit2_PBS(
    std::vector<TLWE<lvl2param>> &out,
    const TLWE<lvl2param> &borrow,
    const TLWE<lvl2param> &sfx_0,
    const EvalKey &ek)
{
    static const std::vector<int> weight = {1, 2};
    static const std::vector<std::vector<int>> lut = [] {
        std::vector<int> ge(4), q(4);
        for (int idx = 0; idx < 4; ++idx) {
            const int bo = idx & 1, s0 = (idx >> 1) & 1;
            ge[idx] = bo ? 0 : 1;
            q[idx] = (!bo && s0) ? 1 : 0;
        }
        return std::vector<std::vector<int>>{ge, q};
    }();
    TLWE<lvl2param> linear;
    Linear<lvl2param>(linear, {borrow, sfx_0}, weight, 4);
    PBS_Multi(out, linear, ek, lut);
}





static const std::vector<int> &masknot_weight(int nb)
{
    static const std::vector<int> w1 = {1, 2};
    static const std::vector<int> w2 = {1, 2, 4};
    static const std::vector<int> w3 = {1, 2, 4, 8};
    return (nb == 1) ? w1 : (nb == 2 ? w2 : w3);
}
static int masknot_tsize(int nb) { return (nb == 1) ? 4 : (nb == 2 ? 8 : 16); }
static const std::vector<std::vector<int>> &masknot_lut(int nb)
{
    static const auto make = [](int n) {
        const int t = (n == 1) ? 4 : (n == 2 ? 8 : 16);
        std::vector<std::vector<int>> lut(n, std::vector<int>(t));
        for (int k = 0; k < n; ++k)
            for (int idx = 0; idx < t; ++idx) {
                const int gv = idx & 1, bv = (idx >> (k + 1)) & 1;
                lut[k][idx] = (bv && gv) ? 0 : 1;
            }
        return lut;
    };
    static const std::vector<std::vector<int>> l1 = make(1);
    static const std::vector<std::vector<int>> l2 = make(2);
    static const std::vector<std::vector<int>> l3 = make(3);
    return (nb == 1) ? l1 : (nb == 2 ? l2 : l3);
}

static void MaskNot_PBS(
    std::vector<TLWE<lvl2param>> &out,
    const TLWE<lvl2param> &ge,
    const std::vector<TLWE<lvl2param>> &bbits,
    const EvalKey &ek)
{
    const int nb = static_cast<int>(bbits.size());
    if (nb < 1 || nb > 3) throw std::runtime_error("MaskNot_PBS: need 1..3 bits");

    std::vector<TLWE<lvl2param>> in;
    in.push_back(ge);
    for (const auto &b : bbits) in.push_back(b);

    TLWE<lvl2param> linear;
    Linear<lvl2param>(linear, in, masknot_weight(nb), masknot_tsize(nb));
    PBS_Multi(out, linear, ek, masknot_lut(nb));
    out.resize(nb);
}


static void HomFA3_PBS(
    std::vector<TLWE<lvl2param>> &cipher_out,
    const std::vector<TLWE<lvl2param>> &cipher_in,
    const EvalKey &ek)
{
    cipher_out.resize(4);
    static const std::vector<std::vector<int>> lut = {
        {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1},
        {0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1},
        {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}};
    static const std::vector<int> weight = {4, 2, 1, 4, 2, 1, 1};
    TLWE<lvl2param> linear;
    Linear<lvl2param>(linear, cipher_in, weight, 16);
    PBS_Multi(cipher_out, linear, ek, lut);
}



static void add_prepared_lsb(
    std::vector<TLWE<lvl2param>> &rem,
    const std::vector<TLWE<lvl2param>> &addend,
    int w,
    const EvalKey &ek)
{
    const int Wpad = ((w + 2) / 3) * 3;
    std::vector<TLWE<lvl2param>> A(Wpad, CT_ZERO), Bv(Wpad, CT_ONE);
    for (int k = 0; k < w; ++k) {
        A[k] = rem[k];
        Bv[k] = addend[k];
    }

    TLWE<lvl2param> carry = CT_ONE;
    for (int blk = 0; blk * 3 < Wpad; ++blk) {
        const int lo = blk * 3;
        std::vector<TLWE<lvl2param>> in7 = {A[lo + 2], A[lo + 1], A[lo],
                                            Bv[lo + 2], Bv[lo + 1], Bv[lo], carry};
        std::vector<TLWE<lvl2param>> out4;
        HomFA3_PBS(out4, in7, ek);
        for (int k = 0; k < 3; ++k)
            if (lo + k < w) rem[lo + k] = out4[k];
        carry = out4[3];
    }
}

static HomDivResult div_impl(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<TLWE<lvl2param>> &B_in,
    int W,
    const EvalKey &ek,
    bool parallel,
    int threads_hint)
{
    if (W <= 0) throw std::runtime_error("HomDiv: W must be positive");
    if (static_cast<int>(A_in.size()) != W || static_cast<int>(B_in.size()) != W) {
        throw std::runtime_error("HomDiv: input size must equal W");
    }

    auto lsb = [&](const std::vector<TLWE<lvl2param>> &V, int k) -> const TLWE<lvl2param> & {
        return V[W - 1 - k];
    };

    
    
    std::vector<TLWE<lvl2param>> sfx(W + 1, CT_ZERO);
    sfx[W - 1] = lsb(B_in, W - 1);
    {
        int k = W - 3;
        for (; k >= 0; k -= 2) {
            hom_trace_next_layer("div_suffix_or");
            std::vector<TLWE<lvl2param>> out2;
            SuffixOr2_PBS(out2, sfx[k + 2], lsb(B_in, k + 1), lsb(B_in, k), ek);
            sfx[k + 1] = out2[0];
            sfx[k] = out2[1];
        }
        if (k == -1) {
            
            hom_trace_next_layer("div_suffix_or");
            std::vector<TLWE<lvl2param>> out2;
            SuffixOr2_PBS(out2, sfx[1], lsb(B_in, 0), CT_ZERO, ek);
            sfx[0] = out2[0];
        }
    }
    const int suffix_depth = hom_trace_current_level();
    hom_trace_set_level(0);

    
    std::vector<TLWE<lvl2param>> rem;
    rem.reserve(W);
    std::vector<TLWE<lvl2param>> quotient(W, CT_ZERO);

    for (int i = 0; i < W; ++i) {
        ScopedHomTracePrefix trace_prefix("div_round_" + std::to_string(i) + "_");
        const int w = i + 1;
        rem.insert(rem.begin(), A_in[i]);  // rem = 2*rem + a_i

        
        std::vector<TLWE<lvl2param>> r_msb(w), b_msb(w);
        for (int k = 0; k < w; ++k) {
            r_msb[w - 1 - k] = rem[k];
            b_msb[w - 1 - k] = lsb(B_in, k);
        }

        
        HomAddResult sub = (parallel && w > 3)
            ? HomSub_unsigned_pal(r_msb, b_msb, w, ek, true, threads_hint)
            : HomSub_unsigned(r_msb, b_msb, w, ek, true);
        if (!sub.ovf_ct.has_value()) {
            throw std::runtime_error("HomDiv: subtraction did not return borrow");
        }

        std::vector<TLWE<lvl2param>> geq;
        hom_trace_ensure_level(suffix_depth);
        hom_trace_next_layer("qbit");
        if (w < W) {
            QBit3_PBS(geq, *sub.ovf_ct, sfx[w], sfx[0], ek);
        } else {
            QBit2_PBS(geq, *sub.ovf_ct, sfx[0], ek);
        }
        const TLWE<lvl2param> ge = geq[0];
        quotient[i] = geq[1];

        if (parallel) {
            
            std::vector<TLWE<lvl2param>> sel = HomSelect(
                r_msb, sub.SUM, ge, w, ek, threads_hint);
            for (int k = 0; k < w; ++k) rem[k] = sel[w - 1 - k];
        } else {
            
            
            std::vector<TLWE<lvl2param>> notmask(w);
            {
                const int ngrp = (w + 2) / 3;
                std::vector<TLWE<lvl2param>> linears(ngrp);
                std::vector<const std::vector<std::vector<int>>*> luts(ngrp);
                std::vector<int> nbs(ngrp);
                for (int g = 0; g < ngrp; ++g) {
                    const int k = g * 3;
                    const int nb = std::min(3, w - k);
                    nbs[g] = nb;
                    std::vector<TLWE<lvl2param>> in;
                    in.push_back(ge);
                    for (int j = 0; j < nb; ++j) in.push_back(lsb(B_in, k + j));
                    Linear<lvl2param>(linears[g], in, masknot_weight(nb),
                                      masknot_tsize(nb));
                    luts[g] = &masknot_lut(nb);
                }
                std::vector<std::vector<TLWE<lvl2param>>> outs;
                PBS_Multi_batch(outs, linears, luts, ek, false, 1);
                for (int g = 0; g < ngrp; ++g) {
                    const int k = g * 3;
                    for (int j = 0; j < nbs[g]; ++j) notmask[k + j] = outs[g][j];
                }
            }
            add_prepared_lsb(rem, notmask, w, ek);
        }
    }

    
    
    HomDivResult result;
    result.quotient = std::move(quotient);
    result.remainder.resize(W);
    for (int k = 0; k < W; ++k) result.remainder[W - 1 - k] = rem[k];
    if (hom_trace_schedule_only()) {
        for (const auto &signal : result.quotient)
            hom_trace_mark_output_level(hom_trace_decode_level(signal[0]));
        for (const auto &signal : result.remainder)
            hom_trace_mark_output_level(hom_trace_decode_level(signal[0]));
    } else {
        hom_trace_mark_outputs();
    }
    return result;
}

} // namespace

namespace TFHEpp {

HomDivResult HomDiv(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<TLWE<lvl2param>> &B_in,
    int W,
    const EvalKey &ek)
{
    return div_impl(A_in, B_in, W, ek, false, 1);
}

HomDivResult HomDiv_pal(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<TLWE<lvl2param>> &B_in,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    return div_impl(A_in, B_in, W, ek, true, threads_hint);
}

} // namespace TFHEpp
