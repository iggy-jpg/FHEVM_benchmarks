#include "tpch_boolean_ops.h"

#include "greaterthan_ks_gpu.h"

namespace tpch16 {

Ciphertext EvalGreaterThan(
    const CryptoContext& cc,
    const Ciphertext& lhs,
    const Ciphertext& rhs,
    double* elapsedMs) {
    auto masks = ks256gpu::BuildAndLoadMasks(cc);
    ks256gpu::TimingReport timing;
    auto result =
        ks256gpu::EvalGreaterThanPOnly(cc, lhs, rhs, masks, &timing);
    if (elapsedMs) {
        *elapsedMs = timing.totalMs;
    }
    return result;
}

}  // namespace tpch16

