#include "HomOP.h"
#include "HomUtil.h"

#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
    using namespace TFHEpp;

    

    
    static const std::vector<std::vector<int>> &fa3_lut()
    {
        static const std::vector<std::vector<int>> lut = {
            {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1},
            {0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1},
            {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1},
            {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}};
        return lut;
    }
    static const std::vector<int> &fa3_weight()
    {
        static const std::vector<int> w = {4, 2, 1, 4, 2, 1, 1};
        return w;
    }
    static const std::vector<std::vector<int>> &pg_lut()
    {
        static const std::vector<std::vector<int>> lut = {
            {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1},
            {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}};
        return lut;
    }
    static const std::vector<int> &pg_weight()
    {
        static const std::vector<int> w = {4, 2, 1, 4, 2, 1};
        return w;
    }

    
    static void HomFA3_unsigned_PBS(std::vector<TLWE<lvl2param>> &cipher_out,
                                    const std::vector<TLWE<lvl2param>> &cipher_in,
                                    const EvalKey &ek)
    {
        cipher_out.resize(4);
        TLWE<lvl2param> linear;
        Linear<lvl2param>(linear, cipher_in, fa3_weight(), 16);
        PBS_Multi(cipher_out, linear, ek, fa3_lut());
    }

    
    
    
    static void pgblock_PBS(std::vector<TLWE<lvl2param>> &outs,
                            const std::vector<TLWE<lvl2param>> &cipher_in,
                            const EvalKey &ek)
    {
        outs.resize(3);
        TLWE<lvl2param> linear;
        Linear<lvl2param>(linear, cipher_in, pg_weight(), 15);
        PBS_Multi(outs, linear, ek, pg_lut());
    }

    
    
    // C2: S = 2·s_hi + s_lo，G' = [S>=4]，GP' = [S>=3]（t7, Σw²=10）
    static void C2_PBS(std::vector<TLWE<lvl2param>> &outs,
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

    // C3: S = 4·s_hi + 2·s_mid + s_lo，G' = [S>=8]，GP' = [S>=7]（t15, Σw²=42）
    static void C3_PBS(std::vector<TLWE<lvl2param>> &outs,
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

    struct GGP {
        TLWE<lvl2param> G, GP;
    };

    
    
    static std::vector<GGP> bk_r3_prefix(
        const std::vector<GGP> &vals, const EvalKey &ek, int recursion_depth = 0)
    {
        const int n = (int)vals.size();
        if (n == 1) return vals;

        struct GDef { int hi, mi, lo; };  
        std::vector<GDef> gdef;
        int i = 0;
        while (i + 3 <= n) { gdef.push_back({i + 2, i + 1, i}); i += 3; }
        if (n - i == 2) gdef.push_back({-1, i + 1, i});
        else if (n - i == 1) gdef.push_back({-1, -1, i});

        std::vector<GGP> groups(gdef.size());
        hom_trace_next_layer("prefix_up_" + std::to_string(recursion_depth));
#pragma omp parallel for if (gdef.size() > 1)
        for (int g = 0; g < (int)gdef.size(); ++g) {
            const GDef &d = gdef[g];
            std::vector<TLWE<lvl2param>> o(2);
            if (d.hi >= 0) {
                C3_PBS(o, vals[d.hi].G, vals[d.hi].GP, vals[d.mi].G, vals[d.mi].GP,
                       vals[d.lo].G, vals[d.lo].GP, ek);
                groups[g] = {o[0], o[1]};
            } else if (d.mi >= 0) {
                C2_PBS(o, vals[d.mi].G, vals[d.mi].GP, vals[d.lo].G, vals[d.lo].GP, ek);
                groups[g] = {o[0], o[1]};
            } else {
                groups[g] = vals[d.lo];
            }
        }

        const std::vector<GGP> pf = bk_r3_prefix(groups, ek, recursion_depth + 1);

        std::vector<GGP> out(n);
        int down_operations = 0;
        for (int j = 0; j < n; ++j) {
            const int g = j / 3, r = j % 3;
            const bool last = (j == n - 1);
            if (r == 2 || (last && r == 1 && n % 3 == 2) ||
                (last && r == 0 && n % 3 == 1) || (g == 0 && r == 0))
                continue;
            ++down_operations;
        }
        if (down_operations > 0)
            hom_trace_next_layer("prefix_down_" + std::to_string(recursion_depth));
#pragma omp parallel for if (n > 1)
        for (int j = 0; j < n; ++j) {
            const int g = j / 3, r = j % 3;
            const bool last = (j == n - 1);
            std::vector<TLWE<lvl2param>> o(2);
            if (r == 2 || (last && r == 1 && n % 3 == 2) ||
                (last && r == 0 && n % 3 == 1)) {
                out[j] = pf[g];
            } else if (g == 0) {
                if (r == 0) {
                    out[j] = vals[0];
                } else {
                    C2_PBS(o, vals[1].G, vals[1].GP, vals[0].G, vals[0].GP, ek);
                    out[j] = {o[0], o[1]};
                }
            } else if (r == 0) {
                C2_PBS(o, vals[j].G, vals[j].GP, pf[g - 1].G, pf[g - 1].GP, ek);
                out[j] = {o[0], o[1]};
            } else {
                C3_PBS(o, vals[j].G, vals[j].GP, vals[j - 1].G, vals[j - 1].GP,
                       pf[g - 1].G, pf[g - 1].GP, ek);
                out[j] = {o[0], o[1]};
            }
        }
        return out;
    }

    
    
    
    static std::vector<TLWE<lvl2param>> carry_from_block_pg_prefix(
        const std::vector<TLWE<lvl2param>> &G,
        const std::vector<TLWE<lvl2param>> &GP,
        bool carry_one,
        const EvalKey &ek)
    {
        const int n = (int)G.size();
        std::vector<GGP> vals(n);
        for (int i = 0; i < n; ++i) vals[i] = {G[i], GP[i]};
        const std::vector<GGP> pref = bk_r3_prefix(vals, ek);

        std::vector<TLWE<lvl2param>> Cin(n + 1);
        Cin[0] = carry_one ? CT_ONE : CT_ZERO;
        for (int i = 0; i < n; ++i) {
            Cin[i + 1] = carry_one ? pref[i].GP : pref[i].G;
        }
        return Cin;
    }
} // namespace

namespace TFHEpp
{
    // ============================================================
    
    // ============================================================
    HomAddResult HomAdd_unsigned(const std::vector<TLWE<lvl2param>> &A_in,
                                 const std::vector<TLWE<lvl2param>> &B_in,
                                 int W, const EvalKey &ek, bool compute_carry)
    {
        if ((int)A_in.size() != W || (int)B_in.size() != W)
            throw std::runtime_error("HomAdd_unsigned: input size must equal W");

        const int Wpad = ((W + 2) / 3) * 3;
        const int pad = Wpad - W;

        std::vector<TLWE<lvl2param>> A(Wpad, CT_ZERO), B(Wpad, CT_ZERO);
        for (int i = 0; i < W; ++i) { A[pad + i] = A_in[i]; B[pad + i] = B_in[i]; }

        std::vector<TLWE<lvl2param>> SUMpad(Wpad);
        TLWE<lvl2param> carry = CT_ZERO;

        for (int blk = 0; blk * 3 < Wpad; ++blk) {
            hom_trace_next_layer("add_serial_block");
            int hi = Wpad - 1 - blk * 3;
            std::vector<TLWE<lvl2param>> in7 = {A[hi-2], A[hi-1], A[hi], B[hi-2], B[hi-1], B[hi], carry}, out4;
            HomFA3_unsigned_PBS(out4, in7, ek);
            for (int k = 0; k < 3; k++) SUMpad[hi - k] = out4[k];
            carry = out4[3];
        }

        HomAddResult res;
        res.SUM.resize(W);
        for (int i = 0; i < W; i++) res.SUM[i] = SUMpad[pad + i];
        if (compute_carry) {
            
            res.ovf_ct = (pad > 0) ? SUMpad[pad - 1] : carry;
        }

        return res;
    }

    // ============================================================
    
    // ============================================================
    HomAddResult HomSub_unsigned(const std::vector<TLWE<lvl2param>> &A_in,
                                 const std::vector<TLWE<lvl2param>> &B_in,
                                 int W, const EvalKey &ek, bool compute_borrow)
    {

        if ((int)A_in.size() != W || (int)B_in.size() != W)
            throw std::runtime_error("HomSub_unsigned: input size must equal W");
        const int Wpad = ((W + 2) / 3) * 3;
        const int pad = Wpad - W;

        std::vector<TLWE<lvl2param>> A(Wpad, CT_ZERO), B_not(Wpad, CT_ONE);
        for (int i = 0; i < W; ++i) {
            A[pad + i] = A_in[i];
            HomNOT<lvl2param>(B_not[pad + i], B_in[i]);
        }

        std::vector<TLWE<lvl2param>> SUMpad(Wpad);
        TLWE<lvl2param> carry = CT_ONE; // A - B = A + NOT(B) + 1

        for (int blk = 0; blk * 3 < Wpad; ++blk) {
            hom_trace_next_layer("sub_serial_block");
            int hi = Wpad - 1 - blk * 3;
            std::vector<TLWE<lvl2param>> in7 = {A[hi-2], A[hi-1], A[hi], B_not[hi-2], B_not[hi-1], B_not[hi], carry}, out4;
            HomFA3_unsigned_PBS(out4, in7, ek);
            for (int k = 0; k < 3; k++) SUMpad[hi - k] = out4[k];
            carry = out4[3];
        }

        HomAddResult res;
        res.SUM.resize(W);
        for (int i = 0; i < W; i++) res.SUM[i] = SUMpad[pad + i];
        if (compute_borrow) {
            
            
            if (pad > 0) {
                res.ovf_ct = SUMpad[pad - 1];
            } else {
                TLWE<lvl2param> borrow;
                HomNOT<lvl2param>(borrow, carry);
                res.ovf_ct = borrow;
            }
        }

        return res;
    }

    // ============================================================
    
    // ============================================================
    HomAddResult HomAdd_unsigned_wide(const std::vector<TLWE<lvl2param>> &A_in,
                                     const std::vector<TLWE<lvl2param>> &B_in,
                                     int W, const EvalKey &ek, bool compute_carry,
                                     int threads_hint)
    {

        if ((int)A_in.size() != W || (int)B_in.size() != W)
            throw std::runtime_error("HomAdd_unsigned_pal: input size must equal W");
        const int Wpad = ((W + 2) / 3) * 3;
        const int pad = Wpad - W;
        std::vector<TLWE<lvl2param>> A(Wpad, CT_ZERO), B(Wpad, CT_ZERO);
        for (int i = 0; i < W; ++i) { A[pad + i] = A_in[i]; B[pad + i] = B_in[i]; }

        const int nBlocks = Wpad / 3;
        const int useT = std::max(1, std::min(nBlocks, pick_thread_count_hint(threads_hint)));
        std::vector<TLWE<lvl2param>> Pblk(nBlocks), Gblk(nBlocks), GPblk(nBlocks);

        
        {
            hom_trace_next_layer("add_pg");
            std::vector<TLWE<lvl2param>> linears(nBlocks);
#pragma omp parallel for schedule(static) if (useT > 1) num_threads(useT)
            for (int blk = 0; blk < nBlocks; ++blk) {
                const int lo = Wpad - 3 * (blk + 1);
                std::vector<TLWE<lvl2param>> in6 = {A[lo], A[lo+1], A[lo+2], B[lo], B[lo+1], B[lo+2]};
                Linear<lvl2param>(linears[blk], in6, pg_weight(), 15);
            }
            std::vector<const std::vector<std::vector<int>>*> luts(nBlocks, &pg_lut());
            std::vector<std::vector<TLWE<lvl2param>>> outs;
            PBS_Multi_batch(outs, linears, luts, ek, useT > 1, useT);
            for (int blk = 0; blk < nBlocks; ++blk) {
                Pblk[blk] = outs[blk][0]; Gblk[blk] = outs[blk][1]; GPblk[blk] = outs[blk][2];
            }
        }

        std::vector<TLWE<lvl2param>> Cin =
            carry_from_block_pg_prefix(Gblk, GPblk, false, ek);
        std::vector<TLWE<lvl2param>> SUMpad(Wpad);

        // HomTrace::begin_layer("add_HomFA3_unsigned_PBS");
        {
            hom_trace_next_layer("add_sum");
            std::vector<TLWE<lvl2param>> linears(nBlocks);
#pragma omp parallel for schedule(static) if (useT > 1) num_threads(useT)
            for (int blk = 0; blk < nBlocks; ++blk) {
                const int lo = Wpad - 3 * (blk + 1);
                std::vector<TLWE<lvl2param>> in7 = {A[lo], A[lo+1], A[lo+2], B[lo], B[lo+1], B[lo+2], Cin[blk]};
                Linear<lvl2param>(linears[blk], in7, fa3_weight(), 16);
            }
            std::vector<const std::vector<std::vector<int>>*> luts(nBlocks, &fa3_lut());
            std::vector<std::vector<TLWE<lvl2param>>> outs;
            PBS_Multi_batch(outs, linears, luts, ek, useT > 1, useT);
            for (int blk = 0; blk < nBlocks; ++blk) {
                const int lo = Wpad - 3 * (blk + 1);
                for (int k = 0; k < 3; ++k) SUMpad[lo+2-k] = outs[blk][k];
            }
        }

        HomAddResult res;
        res.SUM.resize(W);
        for (int i = 0; i < W; ++i) res.SUM[i] = SUMpad[pad + i];
        if (compute_carry) res.ovf_ct = (pad > 0) ? SUMpad[pad - 1] : Cin[nBlocks];

        return res;
    }

    // ============================================================
    
    // ============================================================
    HomAddResult HomSub_unsigned_wide(const std::vector<TLWE<lvl2param>> &A_in,
                                     const std::vector<TLWE<lvl2param>> &B_in,
                                     int W, const EvalKey &ek, bool compute_borrow,
                                     int threads_hint)
    {

        if ((int)A_in.size() != W || (int)B_in.size() != W)
            throw std::runtime_error("HomSub_unsigned_pal: input size must equal W");
        const int Wpad = ((W + 2) / 3) * 3;
        const int pad = Wpad - W;
        std::vector<TLWE<lvl2param>> A(Wpad, CT_ZERO), Bneg(Wpad, CT_ONE);
        for (int i = 0; i < W; ++i) {
            A[pad + i] = A_in[i];
            HomNOT<lvl2param>(Bneg[pad + i], B_in[i]);
        }

        const int nBlocks = Wpad / 3;
        const int useT = std::max(1, std::min(nBlocks, pick_thread_count_hint(threads_hint)));
        std::vector<TLWE<lvl2param>> Pblk(nBlocks), Gblk(nBlocks), GPblk(nBlocks);

        
        {
            hom_trace_next_layer("sub_pg");
            std::vector<TLWE<lvl2param>> linears(nBlocks);
#pragma omp parallel for schedule(static) if (useT > 1) num_threads(useT)
            for (int blk = 0; blk < nBlocks; ++blk) {
                const int lo = Wpad - 3 * (blk + 1);
                std::vector<TLWE<lvl2param>> in6 = {A[lo], A[lo+1], A[lo+2], Bneg[lo], Bneg[lo+1], Bneg[lo+2]};
                Linear<lvl2param>(linears[blk], in6, pg_weight(), 15);
            }
            std::vector<const std::vector<std::vector<int>>*> luts(nBlocks, &pg_lut());
            std::vector<std::vector<TLWE<lvl2param>>> outs;
            PBS_Multi_batch(outs, linears, luts, ek, useT > 1, useT);
            for (int blk = 0; blk < nBlocks; ++blk) {
                Pblk[blk] = outs[blk][0]; Gblk[blk] = outs[blk][1]; GPblk[blk] = outs[blk][2];
            }
        }

        std::vector<TLWE<lvl2param>> Cin =
            carry_from_block_pg_prefix(Gblk, GPblk, true, ek);
        std::vector<TLWE<lvl2param>> SUMpad(Wpad);

        // HomTrace::begin_layer("sub_HomFA3_unsigned_PBS");
        {
            hom_trace_next_layer("sub_sum");
            std::vector<TLWE<lvl2param>> linears(nBlocks);
#pragma omp parallel for schedule(static) if (useT > 1) num_threads(useT)
            for (int blk = 0; blk < nBlocks; ++blk) {
                const int lo = Wpad - 3 * (blk + 1);
                std::vector<TLWE<lvl2param>> in7 = {A[lo], A[lo+1], A[lo+2], Bneg[lo], Bneg[lo+1], Bneg[lo+2], Cin[blk]};
                Linear<lvl2param>(linears[blk], in7, fa3_weight(), 16);
            }
            std::vector<const std::vector<std::vector<int>>*> luts(nBlocks, &fa3_lut());
            std::vector<std::vector<TLWE<lvl2param>>> outs;
            PBS_Multi_batch(outs, linears, luts, ek, useT > 1, useT);
            for (int blk = 0; blk < nBlocks; ++blk) {
                const int lo = Wpad - 3 * (blk + 1);
                for (int k = 0; k < 3; ++k) SUMpad[lo+2-k] = outs[blk][k];
            }
        }

        HomAddResult res;
        res.SUM.resize(W);
        for (int i = 0; i < W; ++i) res.SUM[i] = SUMpad[pad + i];
        if (compute_borrow) {
            
            if (pad > 0) {
                res.ovf_ct = SUMpad[pad - 1];
            } else {
                TLWE<lvl2param> borrow;
                HomNOT<lvl2param>(borrow, Cin[nBlocks]);
                res.ovf_ct = borrow;
            }
        }

        return res;
    }

    
    
    HomAddResult HomAdd_unsigned_pal(const std::vector<TLWE<lvl2param>> &A_in,
                                     const std::vector<TLWE<lvl2param>> &B_in,
                                     int W, const EvalKey &ek, bool compute_carry,
                                     int threads_hint)
    {
        if (wide_parallel_enabled(threads_hint))
            return HomAdd_unsigned_wide(A_in, B_in, W, ek, compute_carry, threads_hint);
        return HomAdd_unsigned(A_in, B_in, W, ek, compute_carry);
    }

    HomAddResult HomSub_unsigned_pal(const std::vector<TLWE<lvl2param>> &A_in,
                                     const std::vector<TLWE<lvl2param>> &B_in,
                                     int W, const EvalKey &ek, bool compute_borrow,
                                     int threads_hint)
    {
        if (wide_parallel_enabled(threads_hint))
            return HomSub_unsigned_wide(A_in, B_in, W, ek, compute_borrow, threads_hint);
        return HomSub_unsigned(A_in, B_in, W, ek, compute_borrow);
    }
} // namespace TFHEpp
