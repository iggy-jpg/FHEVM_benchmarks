
//


//


//     floor(A/d) = (A*m) >> (W+L)




#include "HomOP.h"
#include "HomUtil.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <stdexcept>
#include <vector>

namespace {

using namespace TFHEpp;
using boost::multiprecision::cpp_int;

static cpp_int bits_to_int(const std::vector<uint8_t> &bits)
{
    cpp_int v = 0;
    for (uint8_t b : bits) v = (v << 1) | (b & 1u);
    return v;
}


static void push_constant(
    std::vector<std::vector<TLWE<lvl2param>>> &columns,
    const cpp_int &c,
    int outW)
{
    for (int i = 0; i < outW; ++i)
        if (boost::multiprecision::bit_test(c, i)) columns[i].push_back(CT_ONE);
}





static std::vector<TLWE<lvl2param>> mul_plain_impl(
    const std::vector<TLWE<lvl2param>> &A_in,
    int W,
    const std::vector<uint8_t> &k_bits,
    int outW,
    const EvalKey &ek,
    bool parallel,
    int threads_hint)
{
    const int Wk = static_cast<int>(k_bits.size());
    const cpp_int k = bits_to_int(k_bits);

    
    std::vector<int> naf;  
    {
        cpp_int v = k;
        while (v > 0) {
            if (boost::multiprecision::bit_test(v, 0)) {
                const int di = 2 - static_cast<int>(v % 4);
                naf.push_back(di);
                v -= di;
            } else {
                naf.push_back(0);
            }
            v >>= 1;
        }
    }
    int pop = 0, csd = 0;
    for (int i = 0; i < Wk; ++i) pop += k_bits[i] & 1;
    for (int d : naf) csd += (d != 0);
    const bool use_csd = csd < pop;

    std::vector<std::vector<TLWE<lvl2param>>> columns(outW);
    cpp_int constant = 0;

    if (!use_csd) {
        for (int i = 0; i < Wk; ++i) {
            if (!(k_bits[Wk - 1 - i] & 1)) continue;  
            for (int j = 0; j < W; ++j) {
                const int c = i + j;
                if (c < outW) columns[c].push_back(A_in[W - 1 - j]);
            }
        }
    } else {
        for (int i = 0; i < static_cast<int>(naf.size()); ++i) {
            if (naf[i] == 0 || i >= outW) continue;
            if (naf[i] == 1) {
                for (int j = 0; j < W; ++j) {
                    const int c = i + j;
                    if (c < outW) columns[c].push_back(A_in[W - 1 - j]);
                }
            } else {
                
                
                for (int j = 0; j < W; ++j) {
                    const int c = i + j;
                    if (c < outW) {
                        TLWE<lvl2param> nb;
                        HomNOT<lvl2param>(nb, A_in[W - 1 - j]);
                        columns[c].push_back(nb);
                    }
                }
                for (int c = 0; c < i && c < outW; ++c)
                    constant += cpp_int(1) << c;
                for (int c = std::min(i + W, outW); c < outW; ++c)
                    constant += cpp_int(1) << c;
                constant += 1;
            }
        }
        constant &= (cpp_int(1) << outW) - 1;
        push_constant(columns, constant, outW);
    }
    return HomColumnsCompressAdd(columns, outW, ek, parallel, threads_hint);
}





static const std::vector<std::vector<int>> &andg_lut(int g)
{
    static std::vector<std::vector<std::vector<int>>> cache = [] {
        std::vector<std::vector<std::vector<int>>> c(16);
        for (int g = 2; g <= 15; ++g) {
            const int t = (g >= 8) ? 16 : (g >= 4 ? 8 : 4);
            std::vector<int> row(t, 0);
            row[g] = 1;
            c[g] = {row};
        }
        return c;
    }();
    return cache[g];
}
static int andg_tsize(int g) { return (g >= 8) ? 16 : (g >= 4 ? 8 : 4); }


static const std::vector<std::vector<int>> &leafcmp_lut(int nb, int kv)
{
    static std::vector<std::vector<std::vector<std::vector<int>>>> cache = [] {
        std::vector<std::vector<std::vector<std::vector<int>>>> c(5);
        for (int nb = 1; nb <= 4; ++nb) {
            const int t = 1 << nb;
            c[nb].resize(t);
            for (int kv = 0; kv < t; ++kv) {
                std::vector<int> G(t), GE(t);
                for (int v = 0; v < t; ++v) { G[v] = (v > kv); GE[v] = (v >= kv); }
                c[nb][kv] = {G, GE};
            }
        }
        return c;
    }();
    return cache[nb][kv];
}
static const std::vector<int> &leafcmp_weight(int nb)
{
    static const std::vector<int> w1 = {1}, w2 = {2, 1}, w3 = {4, 2, 1}, w4 = {8, 4, 2, 1};
    return (nb == 1) ? w1 : (nb == 2 ? w2 : (nb == 3 ? w3 : w4));
}


static void AndG_PBS(
    TLWE<lvl2param> &out,
    const std::vector<TLWE<lvl2param>> &bits,
    const EvalKey &ek)
{
    const int g = static_cast<int>(bits.size());
    if (g < 1 || g > 15) throw std::runtime_error("AndG_PBS: need 1..15 bits");
    if (g == 1) { out = bits[0]; return; }
    const int t_size = (g >= 8) ? 16 : (g >= 4 ? 8 : 4);
    std::vector<std::vector<int>> lut(1, std::vector<int>(t_size, 0));
    lut[0][g] = 1;
    const std::vector<int> weight(g, 1);
    TLWE<lvl2param> linear;
    Linear<lvl2param>(linear, bits, weight, t_size);
    std::vector<TLWE<lvl2param>> res;
    PBS_Multi(res, linear, ek, lut);
    out = res[0];
}



static void LeafCmpPlain_PBS(
    std::vector<TLWE<lvl2param>> &out,
    const std::vector<TLWE<lvl2param>> &abits_msb,
    int kv,
    const EvalKey &ek)
{
    const int nb = static_cast<int>(abits_msb.size());
    if (nb < 1 || nb > 4) throw std::runtime_error("LeafCmpPlain_PBS: need 1..4 bits");
    static const int wtab[4][4] = {{1}, {2, 1}, {4, 2, 1}, {8, 4, 2, 1}};
    const int t_size = 1 << nb;
    std::vector<int> weight(wtab[nb - 1], wtab[nb - 1] + nb);
    std::vector<std::vector<int>> lut(2, std::vector<int>(t_size));
    for (int v = 0; v < t_size; ++v) {
        lut[0][v] = (v > kv) ? 1 : 0;
        lut[1][v] = (v >= kv) ? 1 : 0;
    }
    TLWE<lvl2param> linear;
    Linear<lvl2param>(linear, abits_msb, weight, t_size);
    PBS_Multi(out, linear, ek, lut);
    out.resize(2);
}


static void CmpC2_PBS(
    std::vector<TLWE<lvl2param>> &outs,
    const TLWE<lvl2param> &Ghi, const TLWE<lvl2param> &GPhi,
    const TLWE<lvl2param> &Glo, const TLWE<lvl2param> &GPlo,
    const EvalKey &ek)
{
    outs.resize(2);
    static const std::vector<int> weight = {2, 2, 1, 1};
    static const std::vector<std::vector<int>> lut = {
        {0, 0, 0, 0, 1, 1, 1},
        {0, 0, 0, 1, 1, 1, 1}};
    std::vector<TLWE<lvl2param>> in = {Ghi, GPhi, Glo, GPlo};
    TLWE<lvl2param> lin;
    Linear<lvl2param>(lin, in, weight, 7);
    PBS_Multi(outs, lin, ek, lut);
}

static void CmpC3_PBS(
    std::vector<TLWE<lvl2param>> &outs,
    const TLWE<lvl2param> &Ghi, const TLWE<lvl2param> &GPhi,
    const TLWE<lvl2param> &Gmi, const TLWE<lvl2param> &GPmi,
    const TLWE<lvl2param> &Glo, const TLWE<lvl2param> &GPlo,
    const EvalKey &ek)
{
    outs.resize(2);
    static const std::vector<int> weight = {4, 4, 2, 2, 1, 1};
    static const std::vector<std::vector<int>> lut = {
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1},
        {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}};
    std::vector<TLWE<lvl2param>> in = {Ghi, GPhi, Gmi, GPmi, Glo, GPlo};
    TLWE<lvl2param> lin;
    Linear<lvl2param>(lin, in, weight, 15);
    PBS_Multi(outs, lin, ek, lut);
}


static void cmp_plain_tree(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek,
    int threads_hint,
    TLWE<lvl2param> &G_root,
    TLWE<lvl2param> &GE_root)
{
    
    const int nleaf = (W + 3) / 4;
    std::vector<std::array<TLWE<lvl2param>, 2>> nodes(nleaf);
    const int useT = std::max(1, std::min(nleaf, pick_thread_count_hint(threads_hint)));
    
    {
        std::vector<TLWE<lvl2param>> linears(nleaf);
        std::vector<const std::vector<std::vector<int>>*> luts(nleaf);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(useT)
#endif
        for (int i = 0; i < nleaf; ++i) {
            const int base = 4 * i;
            const int nb = std::min(4, W - base);
            std::vector<TLWE<lvl2param>> ab;
            int kv = 0;
            for (int j = 0; j < nb; ++j) {
                ab.push_back(A_in[base + j]);
                kv = (kv << 1) | (k_bits[base + j] & 1);
            }
            Linear<lvl2param>(linears[i], ab, leafcmp_weight(nb), 1 << nb);
            luts[i] = &leafcmp_lut(nb, kv);
        }
        std::vector<std::vector<TLWE<lvl2param>>> outs;
        PBS_Multi_batch(outs, linears, luts, ek, useT > 1, useT);
        for (int i = 0; i < nleaf; ++i) nodes[i] = {outs[i][0], outs[i][1]};
    }

    
    while (nodes.size() > 1) {
        const int n = static_cast<int>(nodes.size());
        std::vector<std::array<TLWE<lvl2param>, 2>> next;
        std::vector<std::array<int, 3>> plan;  // hi, mi(-1), lo(-1)
        int i = 0;
        while (i + 3 <= n) { plan.push_back({i, i + 1, i + 2}); i += 3; }
        if (n - i == 2) { plan.push_back({i, i + 1, -1}); i += 2; }
        else if (n - i == 1) { plan.push_back({i, -1, -1}); i += 1; }
        next.resize(plan.size());
        const int pT = std::max(1, std::min<int>(plan.size(), pick_thread_count_hint(threads_hint)));
        
        {
            static const std::vector<int> w2 = {2, 2, 1, 1};
            static const std::vector<int> w3 = {4, 4, 2, 2, 1, 1};
            static const std::vector<std::vector<int>> lut2 = {
                {0, 0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1, 1}};
            static const std::vector<std::vector<int>> lut3 = {
                {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1},
                {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}};
            const int np = static_cast<int>(plan.size());
            std::vector<TLWE<lvl2param>> linears;
            std::vector<const std::vector<std::vector<int>>*> luts;
            std::vector<int> slot;
            linears.reserve(np); luts.reserve(np); slot.reserve(np);
            for (int g = 0; g < np; ++g) {
                const auto &d = plan[g];
                if (d[1] < 0) { next[g] = nodes[d[0]]; continue; }
                TLWE<lvl2param> lin;
                if (d[2] >= 0) {
                    std::vector<TLWE<lvl2param>> in = {
                        nodes[d[0]][0], nodes[d[0]][1], nodes[d[1]][0],
                        nodes[d[1]][1], nodes[d[2]][0], nodes[d[2]][1]};
                    Linear<lvl2param>(lin, in, w3, 15);
                    luts.push_back(&lut3);
                } else {
                    std::vector<TLWE<lvl2param>> in = {
                        nodes[d[0]][0], nodes[d[0]][1], nodes[d[1]][0], nodes[d[1]][1]};
                    Linear<lvl2param>(lin, in, w2, 7);
                    luts.push_back(&lut2);
                }
                linears.push_back(lin);
                slot.push_back(g);
            }
            if (!linears.empty()) {
                std::vector<std::vector<TLWE<lvl2param>>> outs;
                PBS_Multi_batch(outs, linears, luts, ek, pT > 1, pT);
                for (size_t t = 0; t < slot.size(); ++t)
                    next[slot[t]] = {outs[t][0], outs[t][1]};
            }
        }
        nodes.swap(next);
    }
    G_root = nodes[0][0];
    GE_root = nodes[0][1];
}

} // namespace

namespace TFHEpp {

std::vector<TLWE<lvl2param>> HomMulPlain(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek)
{
    if (W <= 0 || static_cast<int>(A_in.size()) != W)
        throw std::runtime_error("HomMulPlain: input size must equal W");
    if (k_bits.empty())
        throw std::runtime_error("HomMulPlain: k_bits must be non-empty");
    const int outW = W + static_cast<int>(k_bits.size());
    return mul_plain_impl(A_in, W, k_bits, outW, ek, false, 1);
}

std::vector<TLWE<lvl2param>> HomMulPlain_pal(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    if (W <= 0 || static_cast<int>(A_in.size()) != W)
        throw std::runtime_error("HomMulPlain_pal: input size must equal W");
    if (k_bits.empty())
        throw std::runtime_error("HomMulPlain_pal: k_bits must be non-empty");
    const int outW = W + static_cast<int>(k_bits.size());
    return mul_plain_impl(A_in, W, k_bits, outW, ek, true, threads_hint);
}



static HomDivResult div_plain_impl(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &d_bits,
    int W,
    const EvalKey &ek,
    bool parallel,
    int threads_hint,
    bool want_remainder = true)
{
    if (W <= 0 || static_cast<int>(A_in.size()) != W)
        throw std::runtime_error("HomDivPlain: input size must equal W");
    if (static_cast<int>(d_bits.size()) != W)
        throw std::runtime_error("HomDivPlain: d_bits size must equal W");

    const cpp_int d = bits_to_int(d_bits);
    HomDivResult res;

    
    if (d == 0) {
        res.quotient.assign(W, CT_ZERO);
        if (want_remainder) res.remainder = A_in;
        return res;
    }
    // d == 1：q = A, r = 0
    if (d == 1) {
        res.quotient = A_in;
        if (want_remainder) res.remainder.assign(W, CT_ZERO);
        return res;
    }

    int L = 0;
    while ((cpp_int(1) << L) < d) ++L;  // L = ceil(log2 d)

    if ((cpp_int(1) << L) == d) {
        
        res.quotient.assign(W, CT_ZERO);
        for (int j = L; j < W; ++j)
            res.quotient[W - 1 - (j - L)] = A_in[W - 1 - j];
        if (want_remainder) {
            res.remainder.assign(W, CT_ZERO);
            for (int j = 0; j < L; ++j)
                res.remainder[W - 1 - j] = A_in[W - 1 - j];
        }
        return res;
    }

    
    const int K = W + L;
    const cpp_int m = ((cpp_int(1) << K) + d - 1) / d;
    int Wm = 1;
    while ((cpp_int(1) << Wm) <= m) ++Wm;  

    std::vector<uint8_t> m_bits(Wm);
    for (int i = 0; i < Wm; ++i)
        m_bits[Wm - 1 - i] =
            static_cast<uint8_t>(boost::multiprecision::bit_test(m, i));

    
    const int prodW = W + Wm;
    std::vector<TLWE<lvl2param>> prod =
        mul_plain_impl(A_in, W, m_bits, prodW, ek, parallel, threads_hint);

    
    res.quotient.assign(W, CT_ZERO);
    for (int j = 0; j < W; ++j) {
        const int src = K + j;
        if (src < prodW) res.quotient[W - 1 - j] = prod[prodW - 1 - src];
    }

    if (!want_remainder) return res;   

    
    std::vector<uint8_t> d_trim;  
    {
        int msb = W - 1;
        while (msb > 0 && !(d_bits[W - 1 - msb] & 1)) --msb;
        for (int i = msb; i >= 0; --i) d_trim.push_back(d_bits[W - 1 - i]);
    }
    std::vector<TLWE<lvl2param>> qd =
        mul_plain_impl(res.quotient, W, d_trim, W, ek, parallel, threads_hint);
    res.remainder = parallel
        ? HomSub_unsigned_pal(A_in, qd, W, ek, false).SUM
        : HomSub_unsigned(A_in, qd, W, ek, false).SUM;

    return res;
}

HomDivResult HomDivPlain(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &d_bits,
    int W,
    const EvalKey &ek)
{
    return div_plain_impl(A_in, d_bits, W, ek, false, 1);
}

HomDivResult HomDivPlain_pal(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &d_bits,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    return div_plain_impl(A_in, d_bits, W, ek, true, threads_hint);
}

std::vector<TLWE<lvl2param>> HomDivPlainQuotient(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &d_bits,
    int W,
    const EvalKey &ek)
{
    return div_plain_impl(A_in, d_bits, W, ek, false, 1, false).quotient;
}

std::vector<TLWE<lvl2param>> HomDivPlainQuotient_pal(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &d_bits,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    return div_plain_impl(A_in, d_bits, W, ek, true, threads_hint, false).quotient;
}

std::vector<std::vector<TLWE<lvl2param>>> HomMulPlainBatch(
    const std::vector<std::vector<TLWE<lvl2param>>> &As,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    const int N = static_cast<int>(As.size());
    std::vector<std::vector<TLWE<lvl2param>>> out(N);
    if (N == 0) return out;
    const int useT = std::max(1, pick_thread_count_hint(threads_hint));

    if (N >= useT) {
        
        
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(useT)
#endif
        for (int i = 0; i < N; ++i) {
            out[i] = HomMulPlain(As[i], k_bits, W, ek);
        }
    } else {
        for (int i = 0; i < N; ++i) {
            out[i] = HomMulPlain_pal(As[i], k_bits, W, ek, threads_hint);
        }
    }
    return out;
}

std::vector<std::vector<TLWE<lvl2param>>> HomDivPlainQuotientBatch(
    const std::vector<std::vector<TLWE<lvl2param>>> &As,
    const std::vector<uint8_t> &d_bits,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    const int N = static_cast<int>(As.size());
    std::vector<std::vector<TLWE<lvl2param>>> out(N);
    if (N == 0) return out;
    const int useT = std::max(1, pick_thread_count_hint(threads_hint));

    if (N >= useT) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(useT)
#endif
        for (int i = 0; i < N; ++i)
            out[i] = HomDivPlainQuotient(As[i], d_bits, W, ek);
    } else {
        for (int i = 0; i < N; ++i)
            out[i] = HomDivPlainQuotient_pal(As[i], d_bits, W, ek, threads_hint);
    }
    return out;
}

std::vector<HomDivResult> HomDivPlainBatch(
    const std::vector<std::vector<TLWE<lvl2param>>> &As,
    const std::vector<uint8_t> &d_bits,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    const int N = static_cast<int>(As.size());
    std::vector<HomDivResult> out(N);
    if (N == 0) return out;
    const int useT = std::max(1, pick_thread_count_hint(threads_hint));

    if (N >= useT) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(useT)
#endif
        for (int i = 0; i < N; ++i) {
            out[i] = HomDivPlain(As[i], d_bits, W, ek);
        }
    } else {
        for (int i = 0; i < N; ++i) {
            out[i] = HomDivPlain_pal(As[i], d_bits, W, ek, threads_hint);
        }
    }
    return out;
}



std::vector<TLWE<lvl2param>> HomANDPlain(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W)
{
    if ((int)A_in.size() != W || (int)k_bits.size() != W)
        throw std::runtime_error("HomANDPlain: size mismatch");
    std::vector<TLWE<lvl2param>> out(W);
    for (int i = 0; i < W; ++i) out[i] = (k_bits[i] & 1) ? A_in[i] : CT_ZERO;
    return out;
}

std::vector<TLWE<lvl2param>> HomORPlain(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W)
{
    if ((int)A_in.size() != W || (int)k_bits.size() != W)
        throw std::runtime_error("HomORPlain: size mismatch");
    std::vector<TLWE<lvl2param>> out(W);
    for (int i = 0; i < W; ++i) out[i] = (k_bits[i] & 1) ? CT_ONE : A_in[i];
    return out;
}

std::vector<TLWE<lvl2param>> HomXORPlain(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W)
{
    if ((int)A_in.size() != W || (int)k_bits.size() != W)
        throw std::runtime_error("HomXORPlain: size mismatch");
    std::vector<TLWE<lvl2param>> out(W);
    for (int i = 0; i < W; ++i) {
        if (k_bits[i] & 1) {
            HomNOT<lvl2param>(out[i], A_in[i]);
        } else {
            out[i] = A_in[i];
        }
    }
    return out;
}

std::vector<TLWE<lvl2param>> HomSelectPlain(
    const std::vector<uint8_t> &a_bits,
    const std::vector<uint8_t> &b_bits,
    const TLWE<lvl2param> &sel,
    int W)
{
    if ((int)a_bits.size() != W || (int)b_bits.size() != W)
        throw std::runtime_error("HomSelectPlain: size mismatch");
    std::vector<TLWE<lvl2param>> out(W);
    for (int i = 0; i < W; ++i) {
        const int a = a_bits[i] & 1, b = b_bits[i] & 1;
        if (a == b) {
            out[i] = a ? CT_ONE : CT_ZERO;
        } else if (b == 1) {
            out[i] = sel;                      // sel ? 1 : 0
        } else {
            HomNOT<lvl2param>(out[i], sel);    // sel ? 0 : 1
        }
    }
    return out;
}

TLWE<lvl2param> HomEQPlain_pal(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    if ((int)A_in.size() != W || (int)k_bits.size() != W)
        throw std::runtime_error("HomEQPlain: size mismatch");
    
    std::vector<TLWE<lvl2param>> matched(W);
    for (int i = 0; i < W; ++i) {
        if (k_bits[i] & 1) matched[i] = A_in[i];
        else HomNOT<lvl2param>(matched[i], A_in[i]);
    }
    
    while (matched.size() > 1) {
        std::vector<TLWE<lvl2param>> next;
        const int ngrp = (static_cast<int>(matched.size()) + 14) / 15;
        next.resize(ngrp);
        const int useT = std::max(1, std::min(ngrp, pick_thread_count_hint(threads_hint)));
        {
            std::vector<TLWE<lvl2param>> linears(ngrp);
            std::vector<const std::vector<std::vector<int>>*> luts(ngrp);
            std::vector<int> gsz(ngrp);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(useT)
#endif
            for (int g = 0; g < ngrp; ++g) {
                const size_t lo = static_cast<size_t>(g) * 15;
                const size_t hi = std::min(matched.size(), lo + 15);
                std::vector<TLWE<lvl2param>> grp(matched.begin() + lo, matched.begin() + hi);
                gsz[g] = static_cast<int>(grp.size());
                if (gsz[g] == 1) continue;   
                const std::vector<int> weight(gsz[g], 1);
                Linear<lvl2param>(linears[g], grp, weight, andg_tsize(gsz[g]));
                luts[g] = &andg_lut(gsz[g]);
            }
            
            std::vector<TLWE<lvl2param>> bl;
            std::vector<const std::vector<std::vector<int>>*> bluts;
            std::vector<int> slot;
            for (int g = 0; g < ngrp; ++g) {
                if (gsz[g] == 1) { next[g] = matched[static_cast<size_t>(g) * 15]; continue; }
                bl.push_back(linears[g]); bluts.push_back(luts[g]); slot.push_back(g);
            }
            if (!bl.empty()) {
                std::vector<std::vector<TLWE<lvl2param>>> outs;
                PBS_Multi_batch(outs, bl, bluts, ek, useT > 1, useT);
                for (size_t t = 0; t < slot.size(); ++t) next[slot[t]] = outs[t][0];
            }
        }
        matched.swap(next);
    }
    return matched[0];
}

TLWE<lvl2param> HomEQPlain(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek)
{
    return HomEQPlain_pal(A_in, k_bits, W, ek, 1);
}

TLWE<lvl2param> HomGTPlain_pal(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    if ((int)A_in.size() != W || (int)k_bits.size() != W)
        throw std::runtime_error("HomGTPlain: size mismatch");
    TLWE<lvl2param> G, GE;
    cmp_plain_tree(A_in, k_bits, W, ek, threads_hint, G, GE);
    return G;
}

TLWE<lvl2param> HomGTPlain(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek)
{
    return HomGTPlain_pal(A_in, k_bits, W, ek, 1);
}

TLWE<lvl2param> HomGEPlain_pal(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    if ((int)A_in.size() != W || (int)k_bits.size() != W)
        throw std::runtime_error("HomGEPlain: size mismatch");
    TLWE<lvl2param> G, GE;
    cmp_plain_tree(A_in, k_bits, W, ek, threads_hint, G, GE);
    return GE;
}

TLWE<lvl2param> HomGEPlain(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek)
{
    return HomGEPlain_pal(A_in, k_bits, W, ek, 1);
}

TLWE<lvl2param> HomLEPlain_pal(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    TLWE<lvl2param> gt = HomGTPlain_pal(A_in, k_bits, W, ek, threads_hint);
    TLWE<lvl2param> res;
    HomNOT<lvl2param>(res, gt);
    return res;
}

TLWE<lvl2param> HomLEPlain(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek)
{
    return HomLEPlain_pal(A_in, k_bits, W, ek, 1);
}

TLWE<lvl2param> HomLTPlain_pal(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    TLWE<lvl2param> ge = HomGEPlain_pal(A_in, k_bits, W, ek, threads_hint);
    TLWE<lvl2param> res;
    HomNOT<lvl2param>(res, ge);
    return res;
}

TLWE<lvl2param> HomLTPlain(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek)
{
    return HomLTPlain_pal(A_in, k_bits, W, ek, 1);
}

std::vector<TLWE<lvl2param>> HomMulPlainLow(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek)
{
    if (W <= 0 || static_cast<int>(A_in.size()) != W)
        throw std::runtime_error("HomMulPlainLow: input size must equal W");
    if (k_bits.empty())
        throw std::runtime_error("HomMulPlainLow: k_bits must be non-empty");
    return mul_plain_impl(A_in, W, k_bits, W, ek, false, 1);
}

std::vector<TLWE<lvl2param>> HomMulPlainLow_pal(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    if (W <= 0 || static_cast<int>(A_in.size()) != W)
        throw std::runtime_error("HomMulPlainLow_pal: input size must equal W");
    if (k_bits.empty())
        throw std::runtime_error("HomMulPlainLow_pal: k_bits must be non-empty");
    return mul_plain_impl(A_in, W, k_bits, W, ek, true, threads_hint);
}

std::vector<std::vector<TLWE<lvl2param>>> HomMulPlainLowBatch(
    const std::vector<std::vector<TLWE<lvl2param>>> &As,
    const std::vector<uint8_t> &k_bits,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    const int N = static_cast<int>(As.size());
    std::vector<std::vector<TLWE<lvl2param>>> out(N);
    if (N == 0) return out;
    const int useT = std::max(1, pick_thread_count_hint(threads_hint));

    if (N >= useT) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(useT)
#endif
        for (int i = 0; i < N; ++i) {
            out[i] = HomMulPlainLow(As[i], k_bits, W, ek);
        }
    } else {
        for (int i = 0; i < N; ++i) {
            out[i] = HomMulPlainLow_pal(As[i], k_bits, W, ek, threads_hint);
        }
    }
    return out;
}

} // namespace TFHEpp
