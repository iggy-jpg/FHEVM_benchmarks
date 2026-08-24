#pragma once

#include <openfhe.h>

#include <chrono>
#include <cstdint>
#include <stdexcept>

namespace ks256gpu {

using lbcrypto::Ciphertext;
using lbcrypto::CryptoContext;
using lbcrypto::DCRTPoly;

constexpr uint32_t kWordBits   = 16;
constexpr uint32_t kGuardBits  = 16;
constexpr uint32_t kStride     = kWordBits + kGuardBits;  // 32 slots/word
constexpr uint32_t kNumWords   = 2048;
constexpr uint32_t kTotalSlots = kNumWords * kStride;     // 65536 slots

static_assert(
    kTotalSlots <= 65536,
    "This test is configured for ring dimension 131072, so the slot count must not exceed 65536");
static_assert(
    (kTotalSlots & (kTotalSlots - 1U)) == 0U,
    "The CKKS slot count must be a power of two");

struct TimingReport {
    double productMs = 0.0;
    double totalMs   = 0.0;
};

namespace detail {

using Clock = std::chrono::steady_clock;

inline double Milliseconds(const Clock::time_point& begin,
                           const Clock::time_point& end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

inline void Synchronize(const CryptoContext<DCRTPoly>& cc) {
}

template <class Fn>
inline auto TimedGpuCall(const CryptoContext<DCRTPoly>& cc,
                         double* elapsedMs,
                         Fn&& fn) {
    Synchronize(cc);
    const auto begin = Clock::now();
    auto result = fn();
    Synchronize(cc);
    *elapsedMs = Milliseconds(begin, Clock::now());
    return result;
}

}  // namespace detail

/**
 * Evaluates 2,048 independent 16-bit Boolean AND operations.
 *
 * Input/output layout for every packed word:
 *   base+0     = Boolean value (0 or 1),
 *   base+1...511 = zero/inactive slots.
 *
 * This layout is directly compatible with the comparison operators, whose
 * Boolean result is stored in bit 0 of every 512-slot word.
 *
 * For Boolean inputs x,y in {0,1}:
 *
 *   x AND y = xy.
 *
 * The circuit therefore uses exactly one ciphertext-ciphertext
 * multiplication and has ciphertext multiplication depth 1. No rotations or
 * plaintext masks are required. The caller performs exactly one bootstrap
 * after the complete AND circuit.
 */
inline Ciphertext<DCRTPoly> EvalAnd(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& cA,
    const Ciphertext<DCRTPoly>& cB,
    TimingReport* timing = nullptr) {

    if (!cc || !cA || !cB) {
        throw std::invalid_argument(
            "EvalAnd received a null context or ciphertext");
    }

    TimingReport local;
    TimingReport& report = timing ? *timing : local;
    report = TimingReport{};

    detail::Synchronize(cc);
    const auto totalBegin = detail::Clock::now();

    // AND = A * B. This is the only ciphertext-ciphertext multiplication.
    auto result = detail::TimedGpuCall(
        cc,
        &report.productMs,
        [&]() { return cc->EvalMult(cA, cB); });

    report.totalMs =
        detail::Milliseconds(totalBegin, detail::Clock::now());
    return result;
}

}  // namespace ks256gpu




