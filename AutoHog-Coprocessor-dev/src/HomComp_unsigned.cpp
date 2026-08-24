#include "HomOP.h"
#include "HomUtil.h"
#include <omp.h>

namespace
{
    using namespace TFHEpp;

    // ==========================================
    
    // ==========================================

    
    static const std::vector<int> &gt2_weight()
    { static const std::vector<int> w = {1, 2, 4, 8}; return w; }
    static const std::vector<std::vector<int>> &gt2_lut()
    {
        static const std::vector<std::vector<int>> lut = [] {
            std::vector<int> G(16), GP(16);
            for (int i = 0; i < 16; ++i) {
                int ah = (i >> 0) & 1, al = (i >> 1) & 1, bh = (i >> 2) & 1, bl = (i >> 3) & 1;
                int a = (ah << 1) | al, b = (bh << 1) | bl;
                G[i] = (a > b); GP[i] = (a >= b);
            }
            return std::vector<std::vector<int>>{G, GP};
        }();
        return lut;
    }
    static const std::vector<int> &gt3_weight()
    { static const std::vector<int> w = {4, 2, 1, 4, 2, 1}; return w; }
    static const std::vector<std::vector<int>> &gt3_lut()
    {
        static const std::vector<std::vector<int>> lut = {
            {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1},   // G  = [v>=8]
            {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}};  // GP = [v>=7]
        return lut;
    }
    static const std::vector<int> &gt1_weight()
    { static const std::vector<int> w = {1, 2}; return w; }
    static const std::vector<std::vector<int>> &gt1_lut()
    {   // idx = a + 2b：G = [a>b]；GP = [a>=b]
        static const std::vector<std::vector<int>> lut = {{0, 1, 0, 0}, {1, 1, 0, 1}};
        return lut;
    }
    static const std::vector<int> &cmpc2_weight()
    { static const std::vector<int> w = {2, 2, 1, 1}; return w; }
    static const std::vector<std::vector<int>> &cmpc2_lut()
    {
        static const std::vector<std::vector<int>> lut = {
            {0, 0, 0, 0, 1, 1, 1},
            {0, 0, 0, 1, 1, 1, 1}};
        return lut;
    }
    static const std::vector<int> &cmpc3_weight()
    { static const std::vector<int> w = {4, 4, 2, 2, 1, 1}; return w; }
    static const std::vector<std::vector<int>> &cmpc3_lut()
    {
        static const std::vector<std::vector<int>> lut = {
            {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1},
            {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}};
        return lut;
    }

    
    static void GT2_unsigned_PBS(std::vector<TLWE<lvl2param>> &out,
                                 const std::vector<TLWE<lvl2param>> &bits,
                                 const EvalKey &ek)
    {
        out.resize(2);
        TLWE<lvl2param> linear;
        Linear<lvl2param>(linear, bits, gt2_weight(), 16);
        PBS_Multi(out, linear, ek, gt2_lut());
        for (auto& c : out) sanitize_tlwe(c);
    }

    
    
    static void GT3_unsigned_PBS(std::vector<TLWE<lvl2param>> &out,
                                 const std::vector<TLWE<lvl2param>> &a3_msb_first,
                                 const std::vector<TLWE<lvl2param>> &b3_msb_first,
                                 const EvalKey &ek)
    {
        out.resize(2);
        std::vector<TLWE<lvl2param>> in(6);
        in[0] = a3_msb_first[0];
        in[1] = a3_msb_first[1];
        in[2] = a3_msb_first[2];
        HomNOT<lvl2param>(in[3], b3_msb_first[0]);
        HomNOT<lvl2param>(in[4], b3_msb_first[1]);
        HomNOT<lvl2param>(in[5], b3_msb_first[2]);
        TLWE<lvl2param> linear;
        Linear<lvl2param>(linear, in, gt3_weight(), 15);
        PBS_Multi(out, linear, ek, gt3_lut());
        for (auto& c : out) sanitize_tlwe(c);
    }

    
    static void GT1_unsigned_PBS(std::vector<TLWE<lvl2param>> &out,
                                 const TLWE<lvl2param> &a,
                                 const TLWE<lvl2param> &b,
                                 const EvalKey &ek)
    {
        out.resize(2);
        TLWE<lvl2param> linear;
        std::vector<TLWE<lvl2param>> in = {a, b};
        Linear<lvl2param>(linear, in, gt1_weight(), 4);
        PBS_Multi(out, linear, ek, gt1_lut());
        for (auto& c : out) sanitize_tlwe(c);
    }

    
    // C2: S = 2·s_hi + s_lo，G' = [S>=4]，GP' = [S>=3]（t7, Σw²=10）
    static void CmpC2_PBS(std::vector<TLWE<lvl2param>> &outs,
                          const TLWE<lvl2param> &Ghi, const TLWE<lvl2param> &GPhi,
                          const TLWE<lvl2param> &Glo, const TLWE<lvl2param> &GPlo,
                          const EvalKey &ek)
    {
        outs.resize(2);
        std::vector<TLWE<lvl2param>> in = {Ghi, GPhi, Glo, GPlo};
        TLWE<lvl2param> lin;
        Linear<lvl2param>(lin, in, cmpc2_weight(), 7);
        PBS_Multi(outs, lin, ek, cmpc2_lut());
        for (auto& c : outs) sanitize_tlwe(c);
    }

    // C3: S = 4·s_hi + 2·s_mid + s_lo，G' = [S>=8]，GP' = [S>=7]（t15, Σw²=42）
    static void CmpC3_PBS(std::vector<TLWE<lvl2param>> &outs,
                          const TLWE<lvl2param> &Ghi, const TLWE<lvl2param> &GPhi,
                          const TLWE<lvl2param> &Gmi, const TLWE<lvl2param> &GPmi,
                          const TLWE<lvl2param> &Glo, const TLWE<lvl2param> &GPlo,
                          const EvalKey &ek)
    {
        outs.resize(2);
        std::vector<TLWE<lvl2param>> in = {Ghi, GPhi, Gmi, GPmi, Glo, GPlo};
        TLWE<lvl2param> lin;
        Linear<lvl2param>(lin, in, cmpc3_weight(), 15);
        PBS_Multi(outs, lin, ek, cmpc3_lut());
        for (auto& c : outs) sanitize_tlwe(c);
    }

    
    static void GGP_reduce(
        std::vector<std::array<TLWE<lvl2param>, 2>> &nodes,
        TLWE<lvl2param> &G_final,
        const EvalKey &ek,
        int threads_hint = -1)
    {
        while (nodes.size() > 1) {
            const int n = (int)nodes.size();
            std::vector<std::array<int, 3>> plan;
            int i = 0;
            while (i + 3 <= n) { plan.push_back({i, i + 1, i + 2}); i += 3; }
            if (n - i == 2) plan.push_back({i, i + 1, -1});
            else if (n - i == 1) plan.push_back({i, -1, -1});
            std::vector<std::array<TLWE<lvl2param>, 2>> next(plan.size());
            const int useT = std::max(1, std::min<int>((int)plan.size(),
                                                       pick_thread_count_hint(threads_hint)));
            
            {
                const int np = (int)plan.size();
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
                        Linear<lvl2param>(lin, in, cmpc3_weight(), 15);
                        luts.push_back(&cmpc3_lut());
                    } else {
                        std::vector<TLWE<lvl2param>> in = {
                            nodes[d[0]][0], nodes[d[0]][1],
                            nodes[d[1]][0], nodes[d[1]][1]};
                        Linear<lvl2param>(lin, in, cmpc2_weight(), 7);
                        luts.push_back(&cmpc2_lut());
                    }
                    linears.push_back(lin);
                    slot.push_back(g);
                }
                if (!linears.empty()) {
                    std::vector<std::vector<TLWE<lvl2param>>> outs;
                    PBS_Multi_batch(outs, linears, luts, ek, useT > 1, useT);
                    for (size_t t = 0; t < slot.size(); ++t) {
                        sanitize_tlwe(outs[t][0]); sanitize_tlwe(outs[t][1]);
                        next[slot[t]] = {outs[t][0], outs[t][1]};
                    }
                }
            }
            nodes.swap(next);
        }
        G_final = nodes[0][0];
    }

    // ==========================================
    
    // ==========================================

    
    
    
    static void HomFA3_carry_PBS(
        TLWE<lvl2param> &carry_out,
        const std::vector<TLWE<lvl2param>> &cipher_in,
        const EvalKey &ek)
    {
        
        static const std::vector<int> weight = {4, 2, 1, 4, 2, 1, 1};
        
        static const std::vector<std::vector<int>> lut = {{0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}};
        TLWE<lvl2param> linear;
        Linear<lvl2param>(linear, cipher_in, weight, 16);
        PBS_Single(carry_out, linear, ek, lut);
        sanitize_tlwe(carry_out);
    }

} // namespace

namespace TFHEpp
{
    // ==========================================
    
    // ==========================================
    
    
    
    
    TLWE<lvl2param> HomGT_unsigned(const std::vector<TLWE<lvl2param>> &A, 
                                   const std::vector<TLWE<lvl2param>> &B, 
                                   int W, const EvalKey &ek) 
    {
        if ((int)A.size() != W || (int)B.size() != W)
            throw std::runtime_error("HomGT_unsigned: size mismatch");

        
        const int Wpad = ((W + 2) / 3) * 3;
        const int pad = Wpad - W;

        
        std::vector<TLWE<lvl2param>> B_in(Wpad, CT_ZERO);
        std::vector<TLWE<lvl2param>> A_inv(Wpad, CT_ONE); 

        for(int i=0; i<W; ++i) {
            B_in[pad + i] = B[i];
            HomNOT<lvl2param>(A_inv[pad + i], A[i]);
        }

        TLWE<lvl2param> carry = CT_ONE; 

        
        for(int blk=0; blk < Wpad/3; ++blk) {
            int hi = Wpad - 1 - blk*3;
            
            
            std::vector<TLWE<lvl2param>> in7(7);
            in7[0] = B_in[hi-2];   // Weight 4
            in7[1] = B_in[hi-1];   // Weight 2
            in7[2] = B_in[hi];     // Weight 1
            in7[3] = A_inv[hi-2];  // Weight 4
            in7[4] = A_inv[hi-1];  // Weight 2
            in7[5] = A_inv[hi];    // Weight 1
            in7[6] = carry;        // Weight 1

            TLWE<lvl2param> cnext;
            HomFA3_carry_PBS(cnext, in7, ek);
            carry = cnext;
        }

        
        TLWE<lvl2param> res;
        HomNOT<lvl2param>(res, carry);
        return res;
    }

    // ==========================================
    
    // ==========================================
    TLWE<lvl2param> HomGT_unsigned_wide(const std::vector<TLWE<lvl2param>> &A_in,
                                       const std::vector<TLWE<lvl2param>> &B_in,
                                       int W, const EvalKey &ek, int threads_hint)
    {
        if ((int)A_in.size() != W || (int)B_in.size() != W)
            throw std::runtime_error("HomGT_u_pal: size mismatch");

        
        int nblk = (W + 2) / 3;
        std::vector<std::array<TLWE<lvl2param>, 2>> nodes(nblk);

        const int useT = std::max(1, std::min(nblk, pick_thread_count_hint(threads_hint)));
        
        {
            std::vector<TLWE<lvl2param>> linears(nblk);
            std::vector<const std::vector<std::vector<int>>*> luts(nblk);
#pragma omp parallel for schedule(static) num_threads(useT)
            for (int i = 0; i < nblk; ++i) {
                const int base = 3 * i;
                const int nb = std::min(3, W - base);
                if (nb == 3) {
                    std::vector<TLWE<lvl2param>> in(6);
                    in[0] = A_in[base]; in[1] = A_in[base + 1]; in[2] = A_in[base + 2];
                    HomNOT<lvl2param>(in[3], B_in[base]);
                    HomNOT<lvl2param>(in[4], B_in[base + 1]);
                    HomNOT<lvl2param>(in[5], B_in[base + 2]);
                    Linear<lvl2param>(linears[i], in, gt3_weight(), 15);
                    luts[i] = &gt3_lut();
                } else if (nb == 2) {
                    std::vector<TLWE<lvl2param>> in4 = {A_in[base], A_in[base + 1],
                                                        B_in[base], B_in[base + 1]};
                    Linear<lvl2param>(linears[i], in4, gt2_weight(), 16);
                    luts[i] = &gt2_lut();
                } else {
                    std::vector<TLWE<lvl2param>> in2 = {A_in[base], B_in[base]};
                    Linear<lvl2param>(linears[i], in2, gt1_weight(), 4);
                    luts[i] = &gt1_lut();
                }
            }
            std::vector<std::vector<TLWE<lvl2param>>> outs;
            PBS_Multi_batch(outs, linears, luts, ek, useT > 1, useT);
            for (int i = 0; i < nblk; ++i) {
                sanitize_tlwe(outs[i][0]); sanitize_tlwe(outs[i][1]);
                nodes[i] = {outs[i][0], outs[i][1]};  // (G, GP)
            }
        }

        TLWE<lvl2param> G_final;
        GGP_reduce(nodes, G_final, ek, threads_hint);
        return G_final;
    }

    // ==========================================
    
    // ==========================================

    TLWE<lvl2param> HomLT_unsigned_pal(const std::vector<TLWE<lvl2param>> &A, const std::vector<TLWE<lvl2param>> &B, int W, const EvalKey &ek, int threads_hint) {
        return HomGT_unsigned_pal(B, A, W, ek, threads_hint);
    }

    TLWE<lvl2param> HomLT_unsigned(const std::vector<TLWE<lvl2param>> &A, const std::vector<TLWE<lvl2param>> &B, int W, const EvalKey &ek) {
        return HomGT_unsigned(B, A, W, ek);
    }

    TLWE<lvl2param> HomLE_unsigned_pal(const std::vector<TLWE<lvl2param>> &A, const std::vector<TLWE<lvl2param>> &B, int W, const EvalKey &ek, int threads_hint) {
        TLWE<lvl2param> gt = HomGT_unsigned_pal(A, B, W, ek, threads_hint);
        TLWE<lvl2param> res;
        HomNOT<lvl2param>(res, gt);
        return res;
    }

    TLWE<lvl2param> HomLE_unsigned(const std::vector<TLWE<lvl2param>> &A, const std::vector<TLWE<lvl2param>> &B, int W, const EvalKey &ek) {
        TLWE<lvl2param> gt = HomGT_unsigned(A, B, W, ek);
        TLWE<lvl2param> res;
        HomNOT<lvl2param>(res, gt);
        return res;
    }

    TLWE<lvl2param> HomGE_unsigned_pal(const std::vector<TLWE<lvl2param>> &A, const std::vector<TLWE<lvl2param>> &B, int W, const EvalKey &ek, int threads_hint) {
        TLWE<lvl2param> lt = HomLT_unsigned_pal(A, B, W, ek, threads_hint);
        TLWE<lvl2param> res;
        HomNOT<lvl2param>(res, lt);
        return res;
    }

    TLWE<lvl2param> HomGE_unsigned(const std::vector<TLWE<lvl2param>> &A, const std::vector<TLWE<lvl2param>> &B, int W, const EvalKey &ek) {
        TLWE<lvl2param> lt = HomLT_unsigned(A, B, W, ek);
        TLWE<lvl2param> res;
        HomNOT<lvl2param>(res, lt);
        return res;
    }

    
    TLWE<lvl2param> HomGT_unsigned_pal(const std::vector<TLWE<lvl2param>> &A_in,
                                       const std::vector<TLWE<lvl2param>> &B_in,
                                       int W, const EvalKey &ek, int threads_hint)
    {
        
        if (wide_parallel_enabled(threads_hint, 2))
            return HomGT_unsigned_wide(A_in, B_in, W, ek, threads_hint);
        return HomGT_unsigned(A_in, B_in, W, ek);
    }
}

