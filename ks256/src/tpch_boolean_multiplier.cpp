#include "tpch_boolean_ops.h"

#include "multiplier_ks_gpu.h"

namespace tpch16 {

Ciphertext EvalMultiply(
    const CryptoContext& cc,
    const Ciphertext& lhs,
    const Ciphertext& rhs,
    double* elapsedMs) {
    constexpr int directionSign = -1;
    auto masks =
        ks256gpu::BuildAndLoadMultiplierMasks(cc, directionSign);
    ks256gpu::MultiplierTimingReport timing;
    auto result = ks256gpu::EvalMulDaddaKS(
        cc, lhs, rhs, directionSign, masks, &timing);
    if (elapsedMs) {
        *elapsedMs = timing.totalMs;
    }
    return result;
}

}  // namespace tpch16

