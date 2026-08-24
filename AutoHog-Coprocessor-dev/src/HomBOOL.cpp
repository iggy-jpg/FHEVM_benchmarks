#include "HomOP.h"
#include "HomUtil.h"
#include <omp.h> 

namespace
{

    using namespace TFHEpp;

    static void AND2_PBS(
        std::vector<TLWE<lvl2param>> &out,
        const std::vector<TLWE<lvl2param>> &bits,
        const EvalKey &ek)
    {
        static const std::vector<int> weight = {1, 2, 4, 8};
        static const std::vector<std::vector<int>>
            lut = {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1},
                   {0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1}};
        TLWE<lvl2param> linear;
        Linear<lvl2param>(linear, bits, weight, 16);
        PBS_Multi(out, linear, ek, lut);
    }

    static void OR2_PBS(
        std::vector<TLWE<lvl2param>> &out,
        const std::vector<TLWE<lvl2param>> &bits,
        const EvalKey &ek)
    {
        static const std::vector<int> weight = {1, 2, 4, 8};
        static const std::vector<std::vector<int>>
            lut = {{0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                   {0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1}};
        TLWE<lvl2param> linear;
        Linear<lvl2param>(linear, bits, weight, 16);
        PBS_Multi(out, linear, ek, lut);
    }

    static void XOR2_PBS(
        std::vector<TLWE<lvl2param>> &out,
        const std::vector<TLWE<lvl2param>> &bits,
        const EvalKey &ek)
    {
        static const std::vector<int> weight = {1, 2, 4, 8};
        static const std::vector<std::vector<int>>
            lut = {{0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0},
                   {0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0}};
        TLWE<lvl2param> linear;
        Linear<lvl2param>(linear, bits, weight, 16);
        PBS_Multi(out, linear, ek, lut);
    }

    static const std::vector<std::vector<int>> &select_lut()
    {
        static const std::vector<std::vector<int>> lut = {{0, 1, 0, 1, 0, 0, 1, 1}};
        return lut;
    }

    static void SELECT_PBS(
    std::vector<TLWE<lvl2param>>& out,
    const std::vector<TLWE<lvl2param>>& bits,
    const EvalKey& ek)
    {
    static const std::vector<int> weight = {1, 2, 4};

    if (bits.size() != 3) throw std::runtime_error("SELECT_PBS: bits.size()!=3");
    if (out.size() != 1) out.resize(1);
    TLWE<lvl2param> linear;
    Linear<lvl2param>(linear, bits, weight, 8); 
    PBS_Multi(out, linear, ek, select_lut());
    }








} // namespace

namespace TFHEpp
{

    std::vector<TLWE<lvl2param>> HomNOT(
        const std::vector<TLWE<lvl2param>> &A_in,
        int W,
        const EvalKey &)
    {
        if (W < 0 || static_cast<size_t>(W) > A_in.size())
            throw std::runtime_error("HomNOT: invalid W");

        std::vector<TLWE<lvl2param>> out(W);
        for (int i = 0; i < W; i++)
        {
            HomNOT<lvl2param>(out[i], A_in[i]);
        }
        return out;
    }

    std::vector<TLWE<lvl2param>> HomAND(
        const std::vector<TLWE<lvl2param>> &A_in,
        const std::vector<TLWE<lvl2param>> &B_in,
        int W,
        const EvalKey &ek,
        int threads_hint)
    {

        if ((int)A_in.size() != W || (int)B_in.size() != W)
            throw std::runtime_error("HomAND: input size must equal W");
        const int Wpad = ((W + 1) / 2) * 2;
        const int pad = Wpad - W;
        const int nblk = Wpad / 2;

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

        std::vector<TLWE<lvl2param>> out_bits(Wpad);

        const int useT = std::max(1, std::min(nblk, pick_thread_count_hint(threads_hint)));


#pragma omp parallel num_threads(useT)
        {
            std::vector<TLWE<lvl2param>> in4(4), out2(2);

#pragma omp for schedule(static)
            for (int blk = 0; blk < nblk; ++blk)
            {
                const int hi = Wpad - 1 - blk * 2;
                const int lo = hi - 1;
                in4[0] = A[hi];
                in4[1] = A[lo];
                in4[2] = B[hi];
                in4[3] = B[lo];
                AND2_PBS(out2, in4, ek);
                out_bits[hi] = out2[1];
                out_bits[lo] = out2[0];
            }
        }

        std::vector<TLWE<lvl2param>> trimmed(W);
        for (int i = 0; i < W; ++i)
            trimmed[i] = out_bits[pad + i];
        return trimmed;
    }

    std::vector<TLWE<lvl2param>> HomOR(
        const std::vector<TLWE<lvl2param>> &A_in,
        const std::vector<TLWE<lvl2param>> &B_in,
        int W,
        const EvalKey &ek,
        int threads_hint)
    {

        if ((int)A_in.size() != W || (int)B_in.size() != W)
            throw std::runtime_error("HomOR: input size must equal W");
        const int Wpad = ((W + 1) / 2) * 2;
        const int pad = Wpad - W;
        const int nblk = Wpad / 2;

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

        std::vector<TLWE<lvl2param>> out_bits(Wpad);

        const int useT = std::max(1, std::min(nblk, pick_thread_count_hint(threads_hint)));

#pragma omp parallel num_threads(useT)
        {
            std::vector<TLWE<lvl2param>> in4(4), out2(2);

#pragma omp for schedule(static)
            for (int blk = 0; blk < nblk; ++blk)
            {
                const int hi = Wpad - 1 - blk * 2;
                const int lo = hi - 1;
                in4[0] = A[hi];
                in4[1] = A[lo];
                in4[2] = B[hi];
                in4[3] = B[lo];
                OR2_PBS(out2, in4, ek);
                out_bits[hi] = out2[1];
                out_bits[lo] = out2[0];
            }
        }

        std::vector<TLWE<lvl2param>> trimmed(W);
        for (int i = 0; i < W; ++i)
            trimmed[i] = out_bits[pad + i];
        return trimmed;
    }

    std::vector<TLWE<lvl2param>> HomXOR(
        const std::vector<TLWE<lvl2param>> &A_in,
        const std::vector<TLWE<lvl2param>> &B_in,
        int W,
        const EvalKey &ek,
        int threads_hint)
    {

        if ((int)A_in.size() != W || (int)B_in.size() != W)
            throw std::runtime_error("HomXOR: input size must equal W");
        const int Wpad = ((W + 1) / 2) * 2;
        const int pad = Wpad - W;
        const int nblk = Wpad / 2;

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

        std::vector<TLWE<lvl2param>> out_bits(Wpad);

        const int useT = std::max(1, std::min(nblk, pick_thread_count_hint(threads_hint)));

#pragma omp parallel num_threads(useT)
        {
            std::vector<TLWE<lvl2param>> in4(4), out2(2);

#pragma omp for schedule(static)
            for (int blk = 0; blk < nblk; ++blk)
            {
                const int hi = Wpad - 1 - blk * 2;
                const int lo = hi - 1;
                in4[0] = A[hi];
                in4[1] = A[lo];
                in4[2] = B[hi];
                in4[3] = B[lo];
                XOR2_PBS(out2, in4, ek);
                out_bits[hi] = out2[1];
                out_bits[lo] = out2[0];
            }
        }

        std::vector<TLWE<lvl2param>> trimmed(W);
        for (int i = 0; i < W; ++i)
            trimmed[i] = out_bits[pad + i];
        return trimmed;
    }

    std::vector<TLWE<lvl2param>> HomSelect(
        const std::vector<TLWE<lvl2param>> &A_in,
        const std::vector<TLWE<lvl2param>> &B_in,
        const TLWE<lvl2param> &sel,
        int W,
        const EvalKey &ek,
        int threads_hint)
    {

        if ((int)A_in.size() != W || (int)B_in.size() != W)
            throw std::runtime_error("HomSelect: input size must equal W");
        hom_trace_next_layer("select");
        
        const int useT = std::max(1, std::min(W, pick_thread_count_hint(threads_hint)));
        static const std::vector<int> weight = {1, 2, 4};
        std::vector<TLWE<lvl2param>> linears(W);
#pragma omp parallel for schedule(static) num_threads(useT)
        for (int i = 0; i < W; ++i)
        {
            std::vector<TLWE<lvl2param>> in3 = {A_in[i], B_in[i], sel};
            Linear<lvl2param>(linears[i], in3, weight, 8);
        }
        std::vector<const std::vector<std::vector<int>>*> luts(W, &select_lut());
        std::vector<std::vector<TLWE<lvl2param>>> outs;
        PBS_Multi_batch(outs, linears, luts, ek, true, useT);

        std::vector<TLWE<lvl2param>> out_bits(W);
        for (int i = 0; i < W; ++i) out_bits[i] = outs[i][0];
        return out_bits;
    }

} // namespace TFHEpp
