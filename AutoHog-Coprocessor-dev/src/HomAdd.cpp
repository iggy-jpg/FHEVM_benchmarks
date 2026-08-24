#include "HomOP.h"
#include "HomUtil.h"
#include "mem_probe.hpp" 

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{

    using namespace TFHEpp;

    static void Overflow_PBS(TLWE<lvl2param> &ovf_out,
                             const TLWE<lvl2param> &a_m,
                             const TLWE<lvl2param> &b_m,
                             const TLWE<lvl2param> &y_m,
                             const EvalKey &ek)
    {
        static const std::vector<int> weight = {1, 2, 4};
        static const std::vector<std::vector<int>> lut = {{0, 0, 0, 1, 1, 0, 0, 0}};
        std::vector<TLWE<lvl2param>> in3 = {a_m, b_m, y_m};
        TLWE<lvl2param> linear;

        
        // print_mem("[Overflow_PBS] before Linear");
        Linear<lvl2param>(linear, in3, weight, 8);
        // print_mem("[Overflow_PBS] before PBS_Single");
        PBS_Single(ovf_out, linear, ek, lut);
        // print_mem("[Overflow_PBS] after PBS_Single");
    }

    static void OverflowSub_PBS(
        TLWE<lvl2param> &ovf_out,
        const TLWE<lvl2param> &a_m,
        const TLWE<lvl2param> &b_m,
        const TLWE<lvl2param> &y_m,
        const EvalKey &ek)
    {
        static const std::vector<int> weight = {1, 2, 4};
        static const std::vector<std::vector<int>> lut = {{0, 1, 0, 0, 0, 0, 1, 0}};
        std::vector<TLWE<lvl2param>> in3 = {a_m, b_m, y_m};

        TLWE<lvl2param> linear;
        // print_mem("[OverflowSub_PBS] before Linear");
        Linear<lvl2param>(linear, in3, weight, 8);
        // print_mem("[OverflowSub_PBS] before PBS_Single");
        PBS_Single(ovf_out, linear, ek, lut);
        // print_mem("[OverflowSub_PBS] after PBS_Single");
    }

    static void HomFA3_PBS(std::vector<TLWE<lvl2param>> &cipher_out,
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
        // print_mem("[HomFA3_PBS] before Linear");
        Linear<lvl2param>(linear, cipher_in, weight, 16);
        // print_mem("[HomFA3_PBS] before PBS_Multi");
        PBS_Multi(cipher_out, linear, ek, lut);
        // print_mem("[HomFA3_PBS] after PBS_Multi");
    }

    static void pgblock_PBS(std::vector<TLWE<lvl2param>> &outs,
                            const std::vector<TLWE<lvl2param>> &cipher_in,
                            const EvalKey &ek)
    {
        if (cipher_in.size() != 6)
            throw std::runtime_error("pgblock_PBS: cipher_in.size() must be 6 (A3+B3).");

        outs.resize(3); 

        static const std::vector<int> weight = {4, 2, 1, 4, 2, 1};
        static const std::vector<std::vector<int>> lut = {
            {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1},
            {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}};

        TLWE<lvl2param> linear;
        // print_mem("[pgblock_PBS] before Linear");
        Linear<lvl2param>(linear, cipher_in, weight, 15);
        // print_mem("[pgblock_PBS] before PBS_Multi");
        PBS_Multi(outs, linear, ek, lut);
        // print_mem("[pgblock_PBS] after PBS_Multi");
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

    
    static std::vector<GGP> bk_r3_prefix(const std::vector<GGP> &vals, const EvalKey &ek)
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

        const std::vector<GGP> pf = bk_r3_prefix(groups, ek);

        std::vector<GGP> out(n);
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
        for (int i = 0; i < n; ++i)
        {
            Cin[i + 1] = carry_one ? pref[i].GP : pref[i].G;
        }
        return Cin;
    }
} // namespace

namespace TFHEpp
{

    HomAddResult HomAdd(const std::vector<TLWE<lvl2param>> &A_in,
                        const std::vector<TLWE<lvl2param>> &B_in,
                        int W,
                        const EvalKey &ek,
                        bool compute_overflow)
    {
        // [MOD]
        // MemScope ms("[HomAdd serial]");
        // print_mem("[HomAdd] enter");

        if ((int)A_in.size() != W || (int)B_in.size() != W)
            throw std::runtime_error("HomAdd: input size must equal W");

        HomAddResult res;

        const int Wpad = ((W + 2) / 3) * 3;
        const int pad = Wpad - W;

        std::vector<TLWE<lvl2param>> A(Wpad), B(Wpad);
        for (int i = 0; i < pad; ++i)
            A[i] = CT_ZERO, B[i] = CT_ZERO;
        for (int i = 0; i < W; ++i)
            A[pad + i] = A_in[i], B[pad + i] = B_in[i];

        std::vector<TLWE<lvl2param>> SUMpad(Wpad);
        TLWE<lvl2param> carry = CT_ZERO;
        TLWE<lvl2param> ovf_tmp;

        // print_mem("[HomAdd] before FA3 loop");
        for (int blk = 0; blk * 3 < Wpad; ++blk)
        {
            int hi = Wpad - 1 - blk * 3;
            int lo = hi - 2;

            std::vector<TLWE<lvl2param>> in7(7);
            for (int i = 0; i < 3; i++)
                in7[i] = A[lo + i];
            for (int i = 0; i < 3; i++)
                in7[3 + i] = B[lo + i];
            in7[6] = carry;

            std::vector<TLWE<lvl2param>> out4;
            HomFA3_PBS(out4, in7, ek);
            for (int k = 0; k < 3; k++)
                SUMpad[hi - k] = out4[k];
            carry = out4[3];
        }
        // print_mem("[HomAdd] after FA3 loop");

        if (compute_overflow)
        {
            const int msb_idx = pad;
            // print_mem("[HomAdd] before Overflow_PBS");
            Overflow_PBS(ovf_tmp, A[msb_idx], B[msb_idx], SUMpad[msb_idx], ek);
            // print_mem("[HomAdd] after Overflow_PBS");
            res.ovf_ct = ovf_tmp;
        }

        res.SUM.resize(W);
        for (int i = 0; i < W; i++)
            res.SUM[i] = SUMpad[Wpad - W + i];

        // print_mem("[HomAdd] exit");
        return res;
    }

    HomAddResult HomAdd_wide(const std::vector<TLWE<lvl2param>> &A_in,
                            const std::vector<TLWE<lvl2param>> &B_in,
                            int W,
                            const EvalKey &ek,
                            bool compute_overflow,
                            int threads_hint)
    {
        // [MOD]
        // MemScope ms("[HomAdd_par]");
        // print_mem("[HomAdd_par] enter");

        if ((int)A_in.size() != W || (int)B_in.size() != W)
            throw std::runtime_error("HomAdd_pal: input size must equal W");

        HomAddResult res;

        const int Wpad = ((W + 2) / 3) * 3;
        const int pad = Wpad - W;

        std::vector<TLWE<lvl2param>> A(Wpad), B(Wpad);
        for (int i = 0; i < pad; ++i)
            A[i] = CT_ZERO, B[i] = CT_ZERO;
        for (int i = 0; i < W; ++i)
            A[pad + i] = A_in[i], B[pad + i] = B_in[i];

        const int nBlocks = Wpad / 3;

        std::vector<TLWE<lvl2param>> Pblk(nBlocks), Gblk(nBlocks), GPblk(nBlocks);
        {
            // print_mem("[HomAdd_par] PG blocks: before workers");
#pragma omp parallel for if (nBlocks > 1)
            for (int blk = 0; blk < nBlocks; ++blk)
            {
                const int lo = Wpad - 3 * (blk + 1);
                const int hi = lo + 2;

                std::vector<TLWE<lvl2param>> in6(6), pg(3);
                in6[0] = A[lo];
                in6[1] = A[lo + 1];
                in6[2] = A[hi];
                in6[3] = B[lo];
                in6[4] = B[lo + 1];
                in6[5] = B[hi];

                pgblock_PBS(pg, in6, ek);
                Pblk[blk] = pg[0];
                Gblk[blk] = pg[1];
                GPblk[blk] = pg[2];
            }
            // print_mem("[HomAdd_par] PG blocks: after workers");
        }

        // print_mem("[HomAdd_par] before carry_prefix");
        std::vector<TLWE<lvl2param>> Cin =
            carry_from_block_pg_prefix(Gblk, GPblk, false, ek);
        // print_mem("[HomAdd_par] after carry_prefix");

        std::vector<TLWE<lvl2param>> SUMpad(Wpad);
        {
            // print_mem("[HomAdd_par] SUM combine: before workers");
#pragma omp parallel for if (nBlocks > 1)
            for (int blk = 0; blk < nBlocks; ++blk)
            {
                const int lo = Wpad - 3 * (blk + 1);
                const int hi = lo + 2;

                std::vector<TLWE<lvl2param>> in7(7), out4;
                in7[0] = A[lo];
                in7[1] = A[lo + 1];
                in7[2] = A[hi];
                in7[3] = B[lo];
                in7[4] = B[lo + 1];
                in7[5] = B[hi];
                in7[6] = Cin[blk];

                HomFA3_PBS(out4, in7, ek);
                SUMpad[hi] = out4[0];
                SUMpad[hi - 1] = out4[1];
                SUMpad[hi - 2] = out4[2];
            }
            // print_mem("[HomAdd_par] SUM combine: after workers");
        }

        if (compute_overflow)
        {
            const int msb_idx = pad;
            TLWE<lvl2param> ovf;
            // print_mem("[HomAdd_par] before Overflow_PBS");
            Overflow_PBS(ovf, A[msb_idx], B[msb_idx], SUMpad[msb_idx], ek);
            // print_mem("[HomAdd_par] after Overflow_PBS");
            res.ovf_ct = ovf;
        }

        res.SUM.resize(W);
        for (int i = 0; i < W; ++i)
            res.SUM[i] = SUMpad[Wpad - W + i];

        // print_mem("[HomAdd_par] exit");
        return res;
    }

    HomAddResult HomSub(
        const std::vector<TLWE<lvl2param>> &A_in,
        const std::vector<TLWE<lvl2param>> &B_in,
        int W,
        const EvalKey &ek,
        bool compute_overflow)
    {
        // [MOD]
        // MemScope ms("[HomSub serial]");
        // print_mem("[HomSub] enter");

        if ((int)A_in.size() != W || (int)B_in.size() != W)
            throw std::runtime_error("HomSub: A_in/B_in size must equal W");

        HomAddResult res;

        const int Wpad = ((W + 2) / 3) * 3;
        const int pad = Wpad - W;

        const auto &b_m_ct = B_in[0];
        std::vector<TLWE<lvl2param>> A(Wpad), B(Wpad);
        for (int i = 0; i < pad; ++i)
        {
            A[i] = CT_ZERO;
            B[i] = CT_ZERO;
        }
        for (int i = 0; i < W; ++i)
        {
            A[pad + i] = A_in[i];
            B[pad + i] = B_in[i];
        }
        for (size_t i = 0; i < B.size(); ++i)
            HomNOT<lvl2param>(B[i], B[i]);

        std::vector<TLWE<lvl2param>> SUMpad(Wpad);
        TLWE<lvl2param> carry = CT_ONE;
        TLWE<lvl2param> ovf_tmp;

        // print_mem("[HomSub] before FA3 loop");
        for (int blk = 0; blk * 3 < Wpad; ++blk)
        {
            const int hi = Wpad - 1 - blk * 3;
            const int lo = hi - 2;

            std::vector<TLWE<lvl2param>> in7(7);
            for (int i = 0; i < 3; ++i)
                in7[i] = A[lo + i];
            for (int i = 0; i < 3; ++i)
                in7[3 + i] = B[lo + i];
            in7[6] = carry;

            std::vector<TLWE<lvl2param>> out4;
            HomFA3_PBS(out4, in7, ek);

            for (int k = 0; k < 3; ++k)
                SUMpad[hi - k] = out4[k];
            carry = out4[3];
        }
        // print_mem("[HomSub] after FA3 loop");

        if (compute_overflow)
        {
            const int msb_idx = pad;
            const auto &a_m_ct = A[msb_idx];
            const auto &y_m_ct = SUMpad[msb_idx];
            // print_mem("[HomSub] before OverflowSub_PBS");
            OverflowSub_PBS(ovf_tmp, a_m_ct, b_m_ct, y_m_ct, ek);
            // print_mem("[HomSub] after OverflowSub_PBS");
            res.ovf_ct = ovf_tmp;
        }

        res.SUM.resize(W);
        for (int i = 0; i < W; ++i)
            res.SUM[i] = SUMpad[Wpad - W + i];
        // print_mem("[HomSub] exit");
        return res;
    }

    HomAddResult HomSub_wide(const std::vector<TLWE<lvl2param>> &A_in,
                            const std::vector<TLWE<lvl2param>> &B_in,
                            int W,
                            const EvalKey &ek,
                            bool compute_overflow,
                            int threads_hint)
    {
        // [MOD]
        // MemScope ms("[HomSub_par]");
        // print_mem("[HomSub_par] enter");

        if ((int)A_in.size() != W || (int)B_in.size() != W)
            throw std::runtime_error("HomSub_pal: input size must equal W");

        HomAddResult res;

        const int Wpad = ((W + 2) / 3) * 3;
        const int pad = Wpad - W;

        std::vector<TLWE<lvl2param>> A(Wpad), B(Wpad), Bneg(Wpad);
        for (int i = 0; i < pad; ++i)
        {
            A[i] = CT_ZERO;
            B[i] = CT_ZERO;
        }
        for (int i = 0; i < W; ++i)
        {
            A[pad + i] = A_in[i];
            B[pad + i] = B_in[i];
        }
        for (int i = 0; i < Wpad; ++i)
        {
            HomNOT<lvl2param>(Bneg[i], B[i]);
        }

        const int nBlocks = Wpad / 3;

        std::vector<TLWE<lvl2param>> Pblk(nBlocks), Gblk(nBlocks), GPblk(nBlocks);
        {
            // print_mem("[HomSub_par] PG blocks: before workers");
#pragma omp parallel for if (nBlocks > 1)
            for (int blk = 0; blk < nBlocks; ++blk)
            {
                const int lo = Wpad - 3 * (blk + 1);
                const int hi = lo + 2;

                std::vector<TLWE<lvl2param>> in6(6), pg(3);
                in6[0] = A[lo];
                in6[1] = A[lo + 1];
                in6[2] = A[hi];
                in6[3] = Bneg[lo];
                in6[4] = Bneg[lo + 1];
                in6[5] = Bneg[hi];
                pgblock_PBS(pg, in6, ek);
                Pblk[blk] = pg[0];
                Gblk[blk] = pg[1];
                GPblk[blk] = pg[2];
            }
            // print_mem("[HomSub_par] PG blocks: after workers");
        }

        // print_mem("[HomSub_par] before carry_prefix");
        std::vector<TLWE<lvl2param>> Cin =
            carry_from_block_pg_prefix(Gblk, GPblk, true, ek);
        // print_mem("[HomSub_par] after carry_prefix");

        std::vector<TLWE<lvl2param>> SUMpad(Wpad);
        {
            // print_mem("[HomSub_par] SUM combine: before workers");
#pragma omp parallel for if (nBlocks > 1)
            for (int blk = 0; blk < nBlocks; ++blk)
            {
                const int lo = Wpad - 3 * (blk + 1);
                const int hi = lo + 2;
                std::vector<TLWE<lvl2param>> in7(7), out4;
                in7[0] = A[lo];
                in7[1] = A[lo + 1];
                in7[2] = A[hi];
                in7[3] = Bneg[lo];
                in7[4] = Bneg[lo + 1];
                in7[5] = Bneg[hi];
                in7[6] = Cin[blk];
                HomFA3_PBS(out4, in7, ek);
                SUMpad[hi] = out4[0];
                SUMpad[hi - 1] = out4[1];
                SUMpad[hi - 2] = out4[2];
            }
            // print_mem("[HomSub_par] SUM combine: after workers");
        }

        if (compute_overflow)
        {
            const int msb_idx = pad;
            TLWE<lvl2param> ovf;
            // print_mem("[HomSub_par] before OverflowSub_PBS");
            OverflowSub_PBS(ovf, A[msb_idx], B_in[0], SUMpad[msb_idx], ek);
            // print_mem("[HomSub_par] after OverflowSub_PBS");
            res.ovf_ct = ovf;
        }

        res.SUM.resize(W);
        for (int i = 0; i < W; ++i)
            res.SUM[i] = SUMpad[Wpad - W + i];

        // print_mem("[HomSub_par] exit");
        return res;
    }

    HomAddResult HomAdd_pal(const std::vector<TLWE<lvl2param>> &A_in,
                            const std::vector<TLWE<lvl2param>> &B_in,
                            int W, const EvalKey &ek, bool compute_overflow,
                            int threads_hint)
    {
        if (wide_parallel_enabled(threads_hint))
            return HomAdd_wide(A_in, B_in, W, ek, compute_overflow, threads_hint);
        return HomAdd(A_in, B_in, W, ek, compute_overflow);
    }

    HomAddResult HomSub_pal(const std::vector<TLWE<lvl2param>> &A_in,
                            const std::vector<TLWE<lvl2param>> &B_in,
                            int W, const EvalKey &ek, bool compute_overflow,
                            int threads_hint)
    {
        if (wide_parallel_enabled(threads_hint))
            return HomSub_wide(A_in, B_in, W, ek, compute_overflow, threads_hint);
        return HomSub(A_in, B_in, W, ek, compute_overflow);
    }
} // namespace TFHEpp
