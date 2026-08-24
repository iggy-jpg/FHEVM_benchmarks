#include "tpch_boolean_ops.h"

#include "and_gpu.h"

namespace tpch16 {

Ciphertext EvalAnd(
    const CryptoContext& cc,
    const Ciphertext& lhs,
    const Ciphertext& rhs,
    double* elapsedMs) {
    ks256gpu::TimingReport timing;
    auto result = ks256gpu::EvalAnd(cc, lhs, rhs, &timing);
    if (elapsedMs) {
        *elapsedMs = timing.totalMs;
    }
    return result;
}

}  // namespace tpch16

