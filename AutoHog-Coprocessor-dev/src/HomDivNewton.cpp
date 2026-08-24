


//







//        E = 2^(f+1) - t
//        X = (X * E) >> f



#include "HomOP.h"
#include "HomUtil.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using namespace TFHEpp;

constexpr int SEED_BITS = 4;
constexpr int SEED_FRAC = 2 * SEED_BITS;
constexpr int GUARD = 3;      
constexpr int WIN_GUARD = 8;  
constexpr int CORR_ROUNDS = 3;

using CtVec = std::vector<TLWE<lvl2param>>;



static const TLWE<lvl2param> &bit_at_lsb(const CtVec &V, int k)
{
    return V[static_cast<int>(V.size()) - 1 - k];
}


static CtVec slice_bits(const CtVec &V, int lsb_from, int count)
{
    const int n = static_cast<int>(V.size());
    CtVec out(count, CT_ZERO);
    for (int k = 0; k < count; ++k) {
        const int src = lsb_from + k;
        if (src >= 0 && src < n) out[count - 1 - k] = V[n - 1 - src];
    }
    return out;
}


static CtVec append_lsb_zeros(const CtVec &V, int count)
{
    CtVec out = V;
    out.insert(out.end(), count, CT_ZERO);
    return out;
}


static CtVec const_bits(unsigned long long hi, unsigned long long lo, int width)
{
    CtVec out(width, CT_ZERO);
    for (int k = 0; k < width && k < 64; ++k)
        if ((lo >> k) & 1ULL) out[width - 1 - k] = CT_ONE;
    for (int k = 64; k < width && k - 64 < 64; ++k)
        if ((hi >> (k - 64)) & 1ULL) out[width - 1 - k] = CT_ONE;
    return out;
}


static CtVec const_pow2(int bit, int width)
{
    CtVec out(width, CT_ZERO);
    if (bit >= 0 && bit < width) out[width - 1 - bit] = CT_ONE;
    return out;
}




static void SuffixOr2_PBS(
    CtVec &out,
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



static void OneHot2_PBS(
    CtVec &out,
    const TLWE<lvl2param> &s0,
    const TLWE<lvl2param> &s1,
    const TLWE<lvl2param> &s2,
    const EvalKey &ek)
{
    static const std::vector<int> weight = {1, 2, 4};
    static const std::vector<std::vector<int>> lut = [] {
        std::vector<int> h0(8), h1(8);
        for (int idx = 0; idx < 8; ++idx) {
            const int a = idx & 1, b = (idx >> 1) & 1, c = (idx >> 2) & 1;
            h0[idx] = (a && !b) ? 1 : 0;
            h1[idx] = (b && !c) ? 1 : 0;
        }
        return std::vector<std::vector<int>>{h0, h1};
    }();
    TLWE<lvl2param> linear;
    Linear<lvl2param>(linear, {s0, s1, s2}, weight, 8);
    PBS_Multi(out, linear, ek, lut);
}


static void OrG_PBS(
    TLWE<lvl2param> &out,
    const CtVec &bits,
    const EvalKey &ek)
{
    const int g = static_cast<int>(bits.size());
    if (g < 1 || g > 15) throw std::runtime_error("OrG_PBS: need 1..15 bits");
    if (g == 1) {
        out = bits[0];
        return;
    }
    const int t_size = (g >= 8) ? 16 : (g >= 4 ? 8 : 4);
    std::vector<std::vector<int>> lut(1, std::vector<int>(t_size));
    for (int idx = 0; idx < t_size; ++idx) lut[0][idx] = (idx > 0) ? 1 : 0;
    const std::vector<int> weight(g, 1);
    TLWE<lvl2param> linear;
    Linear<lvl2param>(linear, bits, weight, t_size);
    std::vector<TLWE<lvl2param>> res;
    PBS_Multi(res, linear, ek, lut);
    out = res[0];
}


static TLWE<lvl2param> wide_or(const CtVec &bits, const EvalKey &ek)
{
    CtVec cur = bits;
    while (cur.size() > 1) {
        CtVec next;
        for (size_t i = 0; i < cur.size(); i += 15) {
            CtVec grp(cur.begin() + i,
                      cur.begin() + std::min(cur.size(), i + 15));
            TLWE<lvl2param> o;
            OrG_PBS(o, grp, ek);
            next.push_back(o);
        }
        cur.swap(next);
    }
    return cur.empty() ? CT_ZERO : cur[0];
}



static void Seed_PBS(
    CtVec &out,
    const CtVec &idx_bits,  
    const EvalKey &ek)
{
    static const std::vector<int> weight = {1, 2, 4, 8};
    static const std::vector<std::vector<int>> lut = [] {
        const int nb = SEED_FRAC + 2;
        std::vector<std::vector<int>> rows(nb, std::vector<int>(16));
        for (int idx = 0; idx < 16; ++idx) {
            const double b_mid = (1.0 + (idx + 0.5) / 16.0) / 2.0;
            const long long x =
                std::llround((1.0 / b_mid) * static_cast<double>(1 << SEED_FRAC));
            for (int j = 0; j < nb; ++j) rows[j][idx] = (x >> j) & 1;
        }
        return rows;
    }();
    TLWE<lvl2param> linear;
    Linear<lvl2param>(linear, idx_bits, weight, 16);
    PBS_Multi(out, linear, ek, lut);
    
    CtVec msb(out.rbegin(), out.rend());
    out = msb;
}


static void MaskAnd_PBS(
    CtVec &out,
    const TLWE<lvl2param> &ge,
    const CtVec &bbits,
    const EvalKey &ek)
{
    const int nb = static_cast<int>(bbits.size());
    if (nb < 1 || nb > 3) throw std::runtime_error("MaskAnd_PBS: need 1..3 bits");
    const int t_size = (nb == 1) ? 4 : (nb == 2 ? 8 : 16);
    static const std::vector<int> full_weight = {1, 2, 4, 8};
    const std::vector<int> weight(full_weight.begin(), full_weight.begin() + nb + 1);

    std::vector<std::vector<int>> lut(nb, std::vector<int>(t_size));
    for (int k = 0; k < nb; ++k)
        for (int idx = 0; idx < t_size; ++idx) {
            const int gv = idx & 1;
            const int bv = (idx >> (k + 1)) & 1;
            lut[k][idx] = (bv && gv) ? 1 : 0;
        }

    CtVec in;
    in.push_back(ge);
    for (const auto &b : bbits) in.push_back(b);
    TLWE<lvl2param> linear;
    Linear<lvl2param>(linear, in, weight, t_size);
    PBS_Multi(out, linear, ek, lut);
    out.resize(nb);
}


static void And2_PBS(
    TLWE<lvl2param> &out,
    const TLWE<lvl2param> &a,
    const TLWE<lvl2param> &b,
    const EvalKey &ek)
{
    static const std::vector<int> weight = {1, 2};
    static const std::vector<std::vector<int>> lut = {{0, 0, 0, 1}};
    TLWE<lvl2param> linear;
    Linear<lvl2param>(linear, {a, b}, weight, 4);
    std::vector<TLWE<lvl2param>> res;
    PBS_Multi(res, linear, ek, lut);
    out = res[0];
}


static CtVec mask_and_extend(
    const CtVec &B,
    int W,
    const TLWE<lvl2param> &ge,
    int width,
    const EvalKey &ek,
    int useT)
{
    const int ngrp = (W + 2) / 3;
    std::vector<CtVec> grp_out(ngrp);
#pragma omp parallel for schedule(static) num_threads(useT)
    for (int g = 0; g < ngrp; ++g) {
        const int k = g * 3;
        const int nb = std::min(3, W - k);
        CtVec bb;
        for (int j = 0; j < nb; ++j) bb.push_back(bit_at_lsb(B, k + j));
        MaskAnd_PBS(grp_out[g], ge, bb, ek);
    }
    CtVec out(width, CT_ZERO);
    for (int g = 0; g < ngrp; ++g) {
        const int k = g * 3;
        for (size_t j = 0; j < grp_out[g].size(); ++j)
            out[width - 1 - (k + static_cast<int>(j))] = grp_out[g][j];
    }
    return out;
}


static CtVec barrel_shift_left(
    const CtVec &V,
    const CtVec &lzb,  
    const EvalKey &ek,
    int threads_hint)
{
    CtVec cur = V;
    const int n = static_cast<int>(V.size());
    for (int j = static_cast<int>(lzb.size()) - 1; j >= 0; --j) {
        const int s = 1 << j;
        if (s >= n) continue;
        
        CtVec shifted(n, CT_ZERO);
        for (int k = s; k < n; ++k) shifted[n - 1 - k] = cur[n - 1 - (k - s)];
        cur = HomSelect(cur, shifted, lzb[j], n, ek, threads_hint);
    }
    return cur;
}

} // namespace

namespace TFHEpp {

HomDivResult HomDivNewton(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<TLWE<lvl2param>> &B_in,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    if (W <= 0) throw std::runtime_error("HomDivNewton: W must be positive");
    if (static_cast<int>(A_in.size()) != W || static_cast<int>(B_in.size()) != W) {
        throw std::runtime_error("HomDivNewton: input size must equal W");
    }
    if (W < 8) {
        
        return HomDiv_pal(A_in, B_in, W, ek, threads_hint);
    }

    const int useT = std::max(1, pick_thread_count_hint(threads_hint));

    
    CtVec sfx(W + 1, CT_ZERO);
    sfx[W - 1] = bit_at_lsb(B_in, W - 1);
    {
        int k = W - 3;
        for (; k >= 0; k -= 2) {
            CtVec out2;
            SuffixOr2_PBS(out2, sfx[k + 2], bit_at_lsb(B_in, k + 1),
                          bit_at_lsb(B_in, k), ek);
            sfx[k + 1] = out2[0];
            sfx[k] = out2[1];
        }
        if (k == -1) {
            CtVec out2;
            SuffixOr2_PBS(out2, sfx[1], bit_at_lsb(B_in, 0), CT_ZERO, ek);
            sfx[0] = out2[0];
        }
    }
    TLWE<lvl2param> b_is_zero;
    HomNOT<lvl2param>(b_is_zero, sfx[0]);

    
    CtVec h(W, CT_ZERO);
    {
        const int npair = (W + 1) / 2;
#pragma omp parallel for schedule(static) num_threads(useT)
        for (int p = 0; p < npair; ++p) {
            const int k = 2 * p;
            CtVec out2;
            OneHot2_PBS(out2, sfx[k], sfx[k + 1],
                        (k + 2 <= W) ? sfx[k + 2] : CT_ZERO, ek);
            h[k] = out2[0];
            if (k + 1 < W) h[k + 1] = out2[1];
        }
    }

    
    int lz_bits = 0;
    while ((1 << lz_bits) < W) ++lz_bits;
    CtVec lzb(lz_bits, CT_ZERO);
#pragma omp parallel for schedule(dynamic) num_threads(std::max(1, std::min(useT, lz_bits)))
    for (int j = 0; j < lz_bits; ++j) {
        CtVec sel;
        for (int k = 0; k < W; ++k)
            if (((W - 1 - k) >> j) & 1) sel.push_back(h[k]);
        lzb[j] = sel.empty() ? CT_ZERO : wide_or(sel, ek);
    }

    
    const CtVec Bn = barrel_shift_left(B_in, lzb, ek, threads_hint);
    CtVec A_ext2W(2 * W, CT_ZERO);
    for (int k = 0; k < W; ++k) A_ext2W[2 * W - 1 - k] = bit_at_lsb(A_in, k);
    const CtVec Ash = barrel_shift_left(A_ext2W, lzb, ek, threads_hint);

    
    
    CtVec idx_bits = {bit_at_lsb(Bn, W - 5), bit_at_lsb(Bn, W - 4),
                      bit_at_lsb(Bn, W - 3), bit_at_lsb(Bn, W - 2)};
    CtVec X;  
    Seed_PBS(X, idx_bits, ek);
    int frac = SEED_FRAC;

    
    const int F = W + GUARD;
    std::vector<int> precs;
    {
        int acc = SEED_BITS + 1;
        while (acc < F - 2) {
            acc = 2 * acc - 1;
            precs.push_back(acc < F - 2 ? std::min(2 * acc + 3, F) : F);
        }
    }

    for (int f : precs) {
        
        X = append_lsb_zeros(X, f - frac);  
        
        CtVec b_fx = (f <= W) ? slice_bits(Bn, W - f, f)
                              : append_lsb_zeros(Bn, f - W);
        const int mc = std::max(0, f - WIN_GUARD);

        
        CtVec T = HomMulRectWindow(b_fx, f, X, f + 2, mc, f + (f + 2), ek, threads_hint);
        CtVec t = slice_bits(T, f - mc, f + 2);

        
        CtVec E = HomSub_unsigned_pal(const_pow2(f + 1, f + 2), t, f + 2, ek, false).SUM;
        CtVec E1 = slice_bits(E, 0, f + 1);

        
        CtVec XK = HomMulRectWindow(X, f + 2, E1, f + 1, mc, (f + 2) + (f + 1), ek, threads_hint);
        X = slice_bits(XK, f - mc, f + 2);
        frac = f;
    }

    
    const int mc2 = std::max(0, W + F - WIN_GUARD);
    CtVec QK = HomMulRectWindow(Ash, 2 * W, X, F + 2, mc2, 2 * W + F + 2, ek, threads_hint);
    CtVec Qhat = slice_bits(QK, (W + F) - mc2, W + 2);

    const int RW = W + 3;
    CtVec P = HomMulRectWindow(Qhat, W + 2, B_in, W, 0, RW, ek, threads_hint);
    CtVec A_ext(RW, CT_ZERO);
    for (int k = 0; k < W; ++k) A_ext[RW - 1 - k] = bit_at_lsb(A_in, k);
    CtVec R = HomSub_unsigned_pal(A_ext, P, RW, ek, false).SUM;  

    
    CtVec Bext(RW, CT_ZERO);
    for (int k = 0; k < W; ++k) Bext[RW - 1 - k] = bit_at_lsb(B_in, k);

    for (int round = 0; round < CORR_ROUNDS; ++round) {
        // R < 0 -> R += B, Q -= 1
        const TLWE<lvl2param> s_neg = R[0];
        CtVec maskB = mask_and_extend(B_in, W, s_neg, RW, ek, useT);
        R = HomAdd_unsigned_pal(R, maskB, RW, ek, false).SUM;
        {
            CtVec dec(W + 2, CT_ZERO);
            dec[W + 1] = s_neg;
            Qhat = HomSub_unsigned_pal(Qhat, dec, W + 2, ek, false).SUM;
        }

        
        HomAddResult cmp = HomSub_unsigned_pal(R, Bext, RW, ek, true);
        TLWE<lvl2param> not_borrow;
        HomNOT<lvl2param>(not_borrow, *cmp.ovf_ct);
        TLWE<lvl2param> not_sign;
        HomNOT<lvl2param>(not_sign, R[0]);
        TLWE<lvl2param> ge;
        And2_PBS(ge, not_borrow, not_sign, ek);

        CtVec maskB2 = mask_and_extend(B_in, W, ge, RW, ek, useT);
        R = HomSub_unsigned_pal(R, maskB2, RW, ek, false).SUM;
        {
            CtVec inc(W + 2, CT_ZERO);
            inc[W + 1] = ge;
            Qhat = HomAdd_unsigned_pal(Qhat, inc, W + 2, ek, false).SUM;
        }
    }

    
    CtVec zero_w(W, CT_ZERO);
    CtVec Qfin = HomSelect(slice_bits(Qhat, 0, W), zero_w, b_is_zero, W, ek, threads_hint);
    CtVec Rfin = HomSelect(slice_bits(R, 0, W), A_in, b_is_zero, W, ek, threads_hint);

    HomDivResult result;
    result.quotient = std::move(Qfin);
    result.remainder = std::move(Rfin);
    return result;
}

} // namespace TFHEpp
