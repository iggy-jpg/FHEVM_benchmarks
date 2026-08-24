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
    Plaintext activeMask;  // 1 in all 256 active bits, 0 in all guard slots
};

namespace detail {

using Clock = std::chrono::steady_clock;

inline double Milliseconds(const Clock::time_point& begin,
                           const Clock::time_point& end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

inline std::vector<double> MakeActiveMaskValues() {
    std::vector<double> mask(kTotalSlots, 0.0);
    for (uint32_t word = 0; word < kNumWords; ++word) {
        const uint32_t base = word * kStride;
        for (uint32_t bit = 0; bit < kWordBits; ++bit) {
            mask[base + bit] = 1.0;
        }
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

/** Creates and explicitly loads the 256-active-bit plaintext mask. */
inline Masks BuildAndLoadBitwiseNotMasks(
    const CryptoContext<DCRTPoly>& cc) {

    if (!cc) {
        throw std::invalid_argument(
            "BuildAndLoadBitwiseNotMasks received a null context");
    }

    Masks masks;
    masks.activeMask = cc->MakeCKKSPackedPlaintext(
        detail::MakeActiveMaskValues(), 1, 0, nullptr, kTotalSlots);
    masks.activeMask->SetLength(kTotalSlots);
detail::Synchronize(cc);
    return masks;
}

/**
 * Evaluates 128 independent uint256 bitwise NOT operations.
 *
 * Layout:
 *   [256 Boolean data slots][256 zero guard slots] x 128.
 *
 * For every active Boolean bit x:
 *
 *   NOT x = 1 - x.
 *
 * The active mask supplies 1 in each of the 256 data slots and 0 in every
 * guard slot, so the guard region remains zero. The circuit uses no
 * ciphertext-ciphertext multiplication, has multiplication depth 0, and
 * requires no rotations. The caller performs exactly one bootstrap after the
 * complete circuit.
 */
inline Ciphertext<DCRTPoly> EvalBitwiseNot(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& cA,
    Masks& masks,
    TimingReport* timing = nullptr) {

    if (!cc || !cA) {
        throw std::invalid_argument(
            "EvalBitwiseNot received a null context or ciphertext");
    }
    if (!masks.activeMask) {
        throw std::invalid_argument(
            "EvalBitwiseNot received an invalid active mask");
    }

    TimingReport local;
    TimingReport& report = timing ? *timing : local;
    report = TimingReport{};

    detail::Synchronize(cc);
    const auto totalBegin = detail::Clock::now();

    auto result = detail::TimedGpuCall(
        cc,
        &report.subtractMs,
        [&]() { return cc->EvalSub(masks.activeMask, cA); });

    report.totalMs =
        detail::Milliseconds(totalBegin, detail::Clock::now());
    return result;
}

}  // namespace ks256gpu




