#pragma once

#include <openfhe.h>

namespace tpch16 {

using Ciphertext = lbcrypto::Ciphertext<lbcrypto::DCRTPoly>;
using CryptoContext = lbcrypto::CryptoContext<lbcrypto::DCRTPoly>;

constexpr std::uint32_t kWordBits = 16;
constexpr std::uint32_t kGuardBits = 16;
constexpr std::uint32_t kStride = 32;
constexpr std::uint32_t kRecordsPerBatch = 2048;
constexpr std::uint32_t kSlots = 65536;

Ciphertext EvalGreaterThan(
    const CryptoContext& cc,
    const Ciphertext& lhs,
    const Ciphertext& rhs,
    double* elapsedMs = nullptr);

Ciphertext EvalAnd(
    const CryptoContext& cc,
    const Ciphertext& lhs,
    const Ciphertext& rhs,
    double* elapsedMs = nullptr);

Ciphertext EvalMultiply(
    const CryptoContext& cc,
    const Ciphertext& lhs,
    const Ciphertext& rhs,
    double* elapsedMs = nullptr);

Ciphertext EvalAdd(
    const CryptoContext& cc,
    const Ciphertext& lhs,
    const Ciphertext& rhs,
    double* elapsedMs = nullptr);

}  // namespace tpch16

