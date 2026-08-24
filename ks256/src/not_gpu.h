#pragma once

#include <openfhe.h>

#include <chrono>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace ks256gpu {

using lbcrypto::Ciphertext;
using lbcrypto::CryptoContext;
using lbcrypto::DCRTPoly;
using lbcrypto::Plaintext;

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
    double subtractMs = 0.0;
    double totalMs    = 0.0;
};

struct Masks {
    Plaintext resultMask;  // 1 only at bit 0 of every packed word
};

namespace detail {

using Clock = std::chrono::steady_clock;

inline double Milliseconds(const Clock::time_point& begin,
                           const Clock::time_point& end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

inline std::vector<double> MakeResultMaskValues() {
    std::vector<double> mask(kTotalSlots, 0.0);
    for (uint32_t word = 0; word < kNumWords; ++word) {
        mask[word * kStride] = 1.0;
    }
    return mask;
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

/** Creates and explicitly loads the bit-0 Boolean result mask. */
inline Masks BuildAndLoadMasks(const CryptoContext<DCRTPoly>& cc) {
    if (!cc) {
        throw std::invalid_argument(
            "BuildAndLoadMasks received a null context");
    }

    Masks masks;
    masks.resultMask = cc->MakeCKKSPackedPlaintext(
        detail::MakeResultMaskValues(), 1, 0, nullptr, kTotalSlots);
    masks.resultMask->SetLength(kTotalSlots);
detail::Synchronize(cc);
    return masks;
}

/**
 * Evaluates 128 independent Boolean NOT operations.
 *
 * Input/output layout for every packed word:
 *   base+0       = Boolean value (0 or 1),
 *   base+1...511 = zero/inactive slots.
 *
 * For a Boolean input x in {0,1}:
 *
 *   NOT x = 1 - x.
 *
 * The plaintext result mask supplies 1 only in bit 0 of every packed word, so
 * all inactive slots remain zero. The circuit uses no ciphertext-ciphertext
 * multiplications, has multiplication depth 0, and requires no rotations. The
 * caller performs exactly one bootstrap after the complete NOT circuit.
 */
inline Ciphertext<DCRTPoly> EvalNot(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& cA,
    Masks& masks,
    TimingReport* timing = nullptr) {

    if (!cc || !cA) {
        throw std::invalid_argument(
            "EvalNot received a null context or ciphertext");
    }
    if (!masks.resultMask) {
        throw std::invalid_argument(
            "EvalNot received an invalid result mask");
    }

    TimingReport local;
    TimingReport& report = timing ? *timing : local;
    report = TimingReport{};

    detail::Synchronize(cc);
    const auto totalBegin = detail::Clock::now();

    auto result = detail::TimedGpuCall(
        cc,
        &report.subtractMs,
        [&]() { return cc->EvalSub(masks.resultMask, cA); });

    report.totalMs =
        detail::Milliseconds(totalBegin, detail::Clock::now());
    return result;
}

}  // namespace ks256gpu




