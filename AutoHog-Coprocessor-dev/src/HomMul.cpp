#include "HomOP.h"
#include "HomUtil.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using namespace TFHEpp;

static int bit_at_msb_vector(int width, int lsb_index)
{
    return width - 1 - lsb_index;
}

static const std::vector<int> &mul2x2_weight()
{
    static const std::vector<int> w = {1, 2, 4, 8};
    return w;
}

static const std::vector<std::vector<int>> &mul2x2_lut();

static void Mul2x2_PBS(
    std::vector<TLWE<lvl2param>> &out,
    const std::vector<TLWE<lvl2param>> &bits,
    const EvalKey &ek)
{
    if (bits.size() != 4) throw std::runtime_error("Mul2x2_PBS: bits.size() must be 4");
    TLWE<lvl2param> linear;
    Linear<lvl2param>(linear, bits, mul2x2_weight(), 16);
    PBS_Multi(out, linear, ek, mul2x2_lut());
}

static const std::vector<std::vector<int>> &mul2x2_lut()
{
    static const std::vector<std::vector<int>> lut = [] {
        std::vector<int> p0(16), p1(16), p2(16), p3(16);
        for (int idx = 0; idx < 16; ++idx) {
            const int a0 = (idx >> 0) & 1;
            const int a1 = (idx >> 1) & 1;
            const int b0 = (idx >> 2) & 1;
            const int b1 = (idx >> 3) & 1;
            const int product = ((a1 << 1) | a0) * ((b1 << 1) | b0);
            p0[idx] = (product >> 0) & 1;
            p1[idx] = (product >> 1) & 1;
            p2[idx] = (product >> 2) & 1;
            p3[idx] = (product >> 3) & 1;
        }
        return std::vector<std::vector<int>>{p0, p1, p2, p3};
    }();
    return lut;
}




static const std::vector<std::vector<int>> &popcount_lut(int t_size, int nbits)
{
    static const auto make = [](int ts, int nb) {
        std::vector<std::vector<int>> lut(nb, std::vector<int>(ts));
        for (int j = 0; j < nb; ++j)
            for (int idx = 0; idx < ts; ++idx) lut[j][idx] = (idx >> j) & 1;
        return lut;
    };
    static const std::vector<std::vector<int>> lut4 = make(4, 2);
    static const std::vector<std::vector<int>> lut8 = make(8, 3);
    static const std::vector<std::vector<int>> lut16 = make(16, 4);
    if (t_size == 4) return lut4;
    if (t_size == 8) return lut8;
    return lut16;
}


static void PopcountG_PBS(
    std::vector<TLWE<lvl2param>> &out,
    const std::vector<TLWE<lvl2param>> &bits,
    const EvalKey &ek)
{
    const int g = static_cast<int>(bits.size());
    if (g < 2 || g > 15) throw std::runtime_error("PopcountG_PBS: need 2..15 bits");

    const int t_size = (g >= 8) ? 16 : (g >= 4 ? 8 : 4);
    int nbits = 1;
    while ((1 << nbits) <= g) ++nbits;

    const std::vector<int> weight(g, 1);
    TLWE<lvl2param> linear;
    Linear<lvl2param>(linear, bits, weight, t_size);
    PBS_Multi(out, linear, ek, popcount_lut(t_size, nbits));
    out.resize(nbits);
}

static int max_column_height(const std::vector<std::vector<TLWE<lvl2param>>> &columns)
{
    size_t h = 0;
    for (const auto &c : columns) h = std::max(h, c.size());
    return static_cast<int>(h);
}

    
    
    static std::vector<TLWE<lvl2param>> compress_columns_add(
        std::vector<std::vector<TLWE<lvl2param>>> &columns,
        int outW,
        const EvalKey &ek,
        bool parallel,
        int threads_hint)
    {
    
    
    
    int total_bits = 0;
    for (const auto &col : columns) total_bits += (int)col.size();
    const int guard_limit = total_bits + outW + 64;
    for (int guard = 0;; ++guard) {
        if (guard > guard_limit)
            throw std::runtime_error("mul_engine: compression did not converge");

        struct Group {
            int col, start, size;
        };
        std::vector<Group> groups;

        
        
        
        
        
        
        
        const bool flush_tail = max_column_height(columns) == 3;
        for (int c = 0; c < outW; ++c) {
            int idx = 0;
            int remaining = static_cast<int>(columns[c].size());
            while (remaining > 2 || (flush_tail && remaining == 2)) {
                const int g = std::min(15, remaining);
                groups.push_back({c, idx, g});
                idx += g;
                remaining -= g;
            }
        }
        if (groups.empty()) break;

        std::vector<std::vector<TLWE<lvl2param>>> group_out(groups.size());
        hom_trace_next_layer("mul_column_compress_" + std::to_string(guard));
        const int gT = parallel
            ? std::max(1, std::min(static_cast<int>(groups.size()),
                                   pick_thread_count_hint(threads_hint)))
            : 1;
        {
            std::vector<TLWE<lvl2param>> linears(groups.size());
            std::vector<const std::vector<std::vector<int>>*> luts(groups.size());
#pragma omp parallel for schedule(static) if (parallel) num_threads(gT)
            for (int gi = 0; gi < static_cast<int>(groups.size()); ++gi) {
                const auto &g = groups[gi];
                std::vector<TLWE<lvl2param>> in(
                    columns[g.col].begin() + g.start,
                    columns[g.col].begin() + g.start + g.size);
                const int t_size = (g.size >= 8) ? 16 : (g.size >= 4 ? 8 : 4);
                int nbits = 1;
                while ((1 << nbits) <= g.size) ++nbits;
                const std::vector<int> weight(g.size, 1);
                Linear<lvl2param>(linears[gi], in, weight, t_size);
                luts[gi] = &popcount_lut(t_size, nbits);
            }
            PBS_Multi_batch(group_out, linears, luts, ek, parallel, gT);
        }

        std::vector<std::vector<TLWE<lvl2param>>> next(outW);
        std::vector<int> consumed(outW, 0);
        for (const auto &g : groups) consumed[g.col] += g.size;
        for (int c = 0; c < outW; ++c)
            for (size_t k = consumed[c]; k < columns[c].size(); ++k)
                next[c].push_back(columns[c][k]);
        for (size_t gi = 0; gi < groups.size(); ++gi)
            for (size_t j = 0; j < group_out[gi].size(); ++j) {
                const int c = groups[gi].col + static_cast<int>(j);
                if (c < outW) next[c].push_back(group_out[gi][j]);
            }
        columns.swap(next);
    }

    bool any_second_row = false;
    for (int c = 0; c < outW; ++c)
        if (columns[c].size() >= 2) any_second_row = true;

    std::vector<TLWE<lvl2param>> row0(outW, CT_ZERO);
    std::vector<TLWE<lvl2param>> row1(outW, CT_ZERO);
    for (int c = 0; c < outW; ++c) {
        if (!columns[c].empty()) row0[bit_at_msb_vector(outW, c)] = columns[c][0];
        if (columns[c].size() >= 2) row1[bit_at_msb_vector(outW, c)] = columns[c][1];
    }
    if (!any_second_row) return row0;

    if (parallel) {
        
        
        
        if (pick_thread_count_hint(threads_hint) >= 4) {
            return HomAdd_unsigned_wide(row0, row1, outW, ek, false, threads_hint).SUM;
        }
        return HomAdd_unsigned(row0, row1, outW, ek, false).SUM;
    }
    return HomAdd_unsigned(row0, row1, outW, ek, false).SUM;
}





static std::vector<TLWE<lvl2param>> mul_engine_rect(
    const std::vector<TLWE<lvl2param>> &A_in,
    int Wa,
    const std::vector<TLWE<lvl2param>> &B_in,
    int Wb,
    int lo_col,
    int hi_col,
    const EvalKey &ek,
    bool parallel,
    int threads_hint)
{
    const int Ha = (Wa + 1) / 2;
    const int Hb = (Wb + 1) / 2;
    const int outW = hi_col - lo_col;

    struct TileTask {
        int I, J, shift;
    };
    std::vector<TileTask> tiles;
    for (int I = 0; I < Ha; ++I)
        for (int J = 0; J < Hb; ++J) {
            const int shift = 2 * (I + J);
            if (shift < hi_col && shift + 3 >= lo_col) tiles.push_back({I, J, shift});
        }

    const int ntiles = static_cast<int>(tiles.size());
    const int useT = parallel
        ? std::max(1, std::min(std::max(ntiles, 1), pick_thread_count_hint(threads_hint)))
        : 1;

    auto bit_or_zero = [](const std::vector<TLWE<lvl2param>> &V, int width, int lsb) {
        return (lsb < width) ? V[width - 1 - lsb] : CT_ZERO;
    };

    
    
    std::vector<std::vector<TLWE<lvl2param>>> tile_out(tiles.size());
    {
        hom_trace_next_layer("mul_partial_products");
        std::vector<TLWE<lvl2param>> linears(tiles.size());
#pragma omp parallel for schedule(static) if (parallel) num_threads(useT)
        for (int t = 0; t < static_cast<int>(tiles.size()); ++t) {
            const auto &tk = tiles[t];
            std::vector<TLWE<lvl2param>> in = {
                bit_or_zero(A_in, Wa, 2 * tk.I), bit_or_zero(A_in, Wa, 2 * tk.I + 1),
                bit_or_zero(B_in, Wb, 2 * tk.J), bit_or_zero(B_in, Wb, 2 * tk.J + 1)};
            Linear<lvl2param>(linears[t], in, mul2x2_weight(), 16);
        }
        std::vector<const std::vector<std::vector<int>>*> luts(tiles.size(), &mul2x2_lut());
        PBS_Multi_batch(tile_out, linears, luts, ek, parallel, useT);
    }

    
    std::vector<std::vector<TLWE<lvl2param>>> columns(outW);
    for (size_t t = 0; t < tiles.size(); ++t)
        for (int k = 0; k < 4; ++k) {
            const int c = tiles[t].shift + k - lo_col;
            if (c >= 0 && c < outW) columns[c].push_back(tile_out[t][k]);
        }
    tile_out.clear();

    return compress_columns_add(columns, outW, ek, parallel, threads_hint);
}

} // namespace

namespace TFHEpp {

std::vector<TLWE<lvl2param>> HomMulSchoolbook2x2(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<TLWE<lvl2param>> &B_in,
    int W,
    const EvalKey &ek)
{
    if (W <= 0) throw std::runtime_error("HomMul: W must be positive");
    if (static_cast<int>(A_in.size()) != W || static_cast<int>(B_in.size()) != W) {
        throw std::runtime_error("HomMul: input size must equal W");
    }

    const int outW = W;
    std::vector<TLWE<lvl2param>> acc(outW, CT_ZERO);

    for (int ai = 0; ai < W; ai += 2) {
        const TLWE<lvl2param> a0 = A_in[bit_at_msb_vector(W, ai)];
        const TLWE<lvl2param> a1 = (ai + 1 < W) ? A_in[bit_at_msb_vector(W, ai + 1)] : CT_ZERO;

        for (int bi = 0; bi < W; bi += 2) {
            const int shift = ai + bi;
            if (shift >= outW) continue;

            const TLWE<lvl2param> b0 = B_in[bit_at_msb_vector(W, bi)];
            const TLWE<lvl2param> b1 = (bi + 1 < W) ? B_in[bit_at_msb_vector(W, bi + 1)] : CT_ZERO;

            std::vector<TLWE<lvl2param>> partial4;
            Mul2x2_PBS(partial4, {a0, a1, b0, b1}, ek);

            std::vector<TLWE<lvl2param>> shifted(outW, CT_ZERO);
            for (int k = 0; k < 4; ++k) {
                const int out_lsb = shift + k;
                if (out_lsb < outW) {
                    shifted[bit_at_msb_vector(outW, out_lsb)] = partial4[k];
                }
            }

            acc = HomAdd_unsigned(acc, shifted, outW, ek, false).SUM;
        }
    }

    return acc;
}

std::vector<TLWE<lvl2param>> HomMulCompressorLow(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<TLWE<lvl2param>> &B_in,
    int W,
    const EvalKey &ek)
{
    if (W <= 0) throw std::runtime_error("HomMulCompressorLow: W must be positive");
    if (static_cast<int>(A_in.size()) != W || static_cast<int>(B_in.size()) != W) {
        throw std::runtime_error("HomMulCompressorLow: input size must equal W");
    }
    return mul_engine_rect(A_in, W, B_in, W, 0, W, ek, false, 1);
}

std::vector<TLWE<lvl2param>> HomMulCompressorLow_pal(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<TLWE<lvl2param>> &B_in,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    if (W <= 0) throw std::runtime_error("HomMulCompressorLow_pal: W must be positive");
    if (static_cast<int>(A_in.size()) != W || static_cast<int>(B_in.size()) != W) {
        throw std::runtime_error("HomMulCompressorLow_pal: input size must equal W");
    }
    return mul_engine_rect(A_in, W, B_in, W, 0, W, ek, true, threads_hint);
}

std::vector<TLWE<lvl2param>> HomMul_pal(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<TLWE<lvl2param>> &B_in,
    int W,
    const EvalKey &ek,
    int threads_hint)
{
    if (W <= 0) throw std::runtime_error("HomMul_pal: W must be positive");
    if (static_cast<int>(A_in.size()) != W || static_cast<int>(B_in.size()) != W) {
        throw std::runtime_error("HomMul_pal: input size must equal W");
    }
    auto result = mul_engine_rect(A_in, W, B_in, W, 0, W * 2, ek, true, threads_hint);
    if (hom_trace_schedule_only()) {
        for (const auto &signal : result)
            hom_trace_mark_output_level(hom_trace_decode_level(signal[0]));
    } else {
        hom_trace_mark_outputs();
    }
    return result;
}

std::vector<TLWE<lvl2param>> HomMul(
    const std::vector<TLWE<lvl2param>> &A_in,
    const std::vector<TLWE<lvl2param>> &B_in,
    int W,
    const EvalKey &ek)
{
    if (W <= 0) throw std::runtime_error("HomMul: W must be positive");
    if (static_cast<int>(A_in.size()) != W || static_cast<int>(B_in.size()) != W) {
        throw std::runtime_error("HomMul: input size must equal W");
    }
    return mul_engine_rect(A_in, W, B_in, W, 0, W * 2, ek, false, 1);
}

std::vector<TLWE<lvl2param>> HomMulRectWindow(
    const std::vector<TLWE<lvl2param>> &A_in,
    int Wa,
    const std::vector<TLWE<lvl2param>> &B_in,
    int Wb,
    int lo_col,
    int hi_col,
    const EvalKey &ek,
    int threads_hint)
{
    if (Wa <= 0 || Wb <= 0) throw std::runtime_error("HomMulRectWindow: widths must be positive");
    if (static_cast<int>(A_in.size()) != Wa || static_cast<int>(B_in.size()) != Wb) {
        throw std::runtime_error("HomMulRectWindow: input sizes must match widths");
    }
    if (lo_col < 0 || hi_col <= lo_col || hi_col > Wa + Wb) {
        throw std::runtime_error("HomMulRectWindow: invalid column window");
    }
    return mul_engine_rect(A_in, Wa, B_in, Wb, lo_col, hi_col, ek, true, threads_hint);
}

std::vector<TLWE<lvl2param>> HomColumnsCompressAdd(
    std::vector<std::vector<TLWE<lvl2param>>> &columns,
    int outW,
    const EvalKey &ek,
    bool parallel,
    int threads_hint)
{
    if (outW <= 0 || (int)columns.size() != outW) {
        throw std::runtime_error("HomColumnsCompressAdd: columns.size() must equal outW");
    }
    return compress_columns_add(columns, outW, ek, parallel, threads_hint);
}

} // namespace TFHEpp
