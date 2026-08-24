#include "HomOP.h"
#include "HomUtil.h"







namespace
{
    using namespace TFHEpp;

    static void flip_msb(std::vector<TLWE<lvl2param>> &v)
    {
        if (!v.empty()) HomNOT<lvl2param>(v[0], v[0]);
    }
} // anonymous namespace

namespace TFHEpp
{

    TLWE<lvl2param> HomGT(
        const std::vector<TLWE<lvl2param>> &A_in,
        const std::vector<TLWE<lvl2param>> &B_in,
        int W,
        const EvalKey &ek)
    {
        if ((int)A_in.size() != W || (int)B_in.size() != W)
            throw std::runtime_error("HomGT: A_in/B_in size must equal W");
        std::vector<TLWE<lvl2param>> A = A_in, B = B_in;
        flip_msb(A);
        flip_msb(B);
        return HomGT_unsigned(A, B, W, ek);
    }

    TLWE<lvl2param> HomLE(
        const std::vector<TLWE<lvl2param>> &A_in,
        const std::vector<TLWE<lvl2param>> &B_in,
        int W,
        const EvalKey &ek)
    {
        TLWE<lvl2param> gt = HomGT(A_in, B_in, W, ek);
        TLWE<lvl2param> res;
        HomNOT<lvl2param>(res, gt);
        return res;
    }

    TLWE<lvl2param> HomGT_pal(
        const std::vector<TLWE<lvl2param>> &A_in,
        const std::vector<TLWE<lvl2param>> &B_in,
        int W,
        const EvalKey &ek,
        int threads_hint)
    {
        if ((int)A_in.size() != W || (int)B_in.size() != W)
            throw std::runtime_error("HomGT_pal: A_in/B_in size must equal W");
        std::vector<TLWE<lvl2param>> A = A_in, B = B_in;
        flip_msb(A);
        flip_msb(B);
        return HomGT_unsigned_pal(A, B, W, ek, threads_hint);
    }

    TLWE<lvl2param> HomLE_pal(
        const std::vector<TLWE<lvl2param>> &A_in,
        const std::vector<TLWE<lvl2param>> &B_in,
        int W,
        const EvalKey &ek,
        int threads_hint)
    {
        TLWE<lvl2param> gt = HomGT_pal(A_in, B_in, W, ek, threads_hint);
        TLWE<lvl2param> res;
        HomNOT<lvl2param>(res, gt);
        return res;
    }

    TLWE<lvl2param> HomLT(const std::vector<TLWE<lvl2param>> &A,
                          const std::vector<TLWE<lvl2param>> &B,
                          int W, const EvalKey &ek)
    {
        return HomGT(B, A, W, ek);
    }

    TLWE<lvl2param> HomGE(const std::vector<TLWE<lvl2param>> &A,
                          const std::vector<TLWE<lvl2param>> &B,
                          int W, const EvalKey &ek)
    {
        return HomLE(B, A, W, ek);
    }

    TLWE<lvl2param> HomLT_pal(const std::vector<TLWE<lvl2param>> &A,
                              const std::vector<TLWE<lvl2param>> &B,
                              int W, const EvalKey &ek, int threads_hint)
    {
        return HomGT_pal(B, A, W, ek, threads_hint);
    }

    TLWE<lvl2param> HomGE_pal(const std::vector<TLWE<lvl2param>> &A,
                              const std::vector<TLWE<lvl2param>> &B,
                              int W, const EvalKey &ek, int threads_hint)
    {
        return HomLE_pal(B, A, W, ek, threads_hint);
    }
} // namespace TFHEpp
