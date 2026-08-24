#pragma once

#include <openfhe.h>

#include <chrono>
#include <cstdint>
#include <stdexcept>

namespace ks256gpu {

using lbcrypto::Ciphertext;
using lbcrypto::CryptoContext;
using lbcrypto::DCRTPoly;

constexpr uint32_t kWordBits   = 256;
constexpr uint32_t kGuardBits  = 256;
constexpr uint32_t kStride     = kWordBits + kGuardBits;  // 512 slots/word
constexpr uint32_t kNumWords   = 128;
constexpr uint32_t kTotalSlots = kNumWords * kStride;     // 65536 slots

static_assert(
    kTotalSlots <= 65536,
    "This test is configured for ring dimension 131072, so the slot count must not exceed 65536");
static_assert(
    (kTotalSlots & (kTotalSlots - 1U)) == 0U,
    "The CKKS slot count must be a power of two");

struct TimingReport {
    double inputAddMs       = 0.0;
    double productMs        = 0.0;
    double doubleProductMs  = 0.0;
    double subtractMs       = 0.0;
    double totalMs          = 0.0;
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
 * Evaluates 128 independent bitwise XOR operations on packed uint256 values.
 *
 * Layout for both ciphertexts:
 *   [256 Boolean data slots][256 zero guard slots] x 128.
 *
 * Bit order inside every word:
 *   base+0   = bit 0  (LSB),
 *   base+255 = bit 255 (MSB).
 *
 * For Boolean inputs x,y in {0,1}, XOR is evaluated as
 *
 *   x XOR y = x + y - 2xy.
 *
 * The circuit uses exactly one ciphertext-ciphertext multiplication and has
 * ciphertext multiplication depth 1. The product xy is doubled with a
 * ciphertext addition, so no extra multiplication level is consumed.
 *
 * Because all guard slots in both inputs are zero, the formula also leaves
 * every guard slot at zero. No rotations or plaintext masks are required.
 * The caller performs exactly one bootstrap after the complete XOR circuit.
 */
inline Ciphertext<DCRTPoly> EvalXor(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& cA,
    const Ciphertext<DCRTPoly>& cB,
    TimingReport* timing = nullptr) {

    if (!cc || !cA || !cB) {
        throw std::invalid_argument(
            "EvalXor received a null context or ciphertext");
    }

    TimingReport local;
    TimingReport& report = timing ? *timing : local;
    report = TimingReport{};

    detail::Synchronize(cc);
    const auto totalBegin = detail::Clock::now();

    // sum = A + B.
    auto sum = detail::TimedGpuCall(
        cc,
        &report.inputAddMs,
        [&]() { return cc->EvalAdd(cA, cB); });

    // product = A * B. This is the only ciphertext-ciphertext multiplication.
    auto product = detail::TimedGpuCall(
        cc,
        &report.productMs,
        [&]() { return cc->EvalMult(cA, cB); });

    // twiceProduct = 2AB, implemented as AB + AB without increasing depth.
    auto twiceProduct = detail::TimedGpuCall(
        cc,
        &report.doubleProductMs,
        [&]() { return cc->EvalAdd(product, product); });

    // XOR = A + B - 2AB.
    auto result = detail::TimedGpuCall(
        cc,
        &report.subtractMs,
        [&]() { return cc->EvalSub(sum, twiceProduct); });

    report.totalMs =
        detail::Milliseconds(totalBegin, detail::Clock::now());
    return result;
}

}  // namespace ks256gpu




