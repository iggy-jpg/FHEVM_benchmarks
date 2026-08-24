#include "HomOP.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace TFHEpp {

std::vector<TLWE<lvl2param>> HomCast(
    const std::vector<TLWE<lvl2param>> &A_in,
    int fromW,
    int toW,
    bool signed_extend)
{
    if (fromW <= 0) throw std::runtime_error("HomCast: fromW must be positive");
    if (toW <= 0) throw std::runtime_error("HomCast: toW must be positive");
    if (static_cast<int>(A_in.size()) != fromW) {
        throw std::runtime_error("HomCast: input size must equal fromW");
    }

    std::vector<TLWE<lvl2param>> out(toW);
    const int copied = std::min(fromW, toW);
    const int src_offset = fromW - copied;
    const int dst_offset = toW - copied;

    for (int i = 0; i < dst_offset; ++i) {
        out[i] = signed_extend ? A_in[0] : CT_ZERO;
    }
    for (int i = 0; i < copied; ++i) {
        out[dst_offset + i] = A_in[src_offset + i];
    }

    return out;
}

} // namespace TFHEpp
