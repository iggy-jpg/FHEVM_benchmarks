#include "tpch_boolean_ops.h"

#include "adder_ks_gpu.h"

namespace tpch16 {

Ciphertext EvalAdd(
    const CryptoContext& cc,
    const Ciphertext& lhs,
    const Ciphertext& rhs,
    double* elapsedMs) {
    constexpr int directionSign = -1;
    auto masks = ks256gpu::BuildAndLoadMasks(cc, directionSign);
    ks256gpu::TimingReport timing;
    auto result =
        ks256gpu::EvalAddKS(cc, lhs, rhs, directionSign, masks, &timing);
    if (elapsedMs) {
        *elapsedMs = timing.totalMs;
    }
    return result;
}

}  // namespace tpch16

