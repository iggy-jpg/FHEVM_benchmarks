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
    double rotateMs = 0.0;
    double maskMs   = 0.0;
    double zeroMs   = 0.0;
    double totalMs  = 0.0;
};

struct Masks {
    // 1 on the 256 Boolean data slots and 0 on all guard slots.
    Plaintext activeMask;

    // For a public shift s<256:
    //   1 on result bits [0,255-s],
    //   0 on result bits [256-s,255],
    //   0 on all guard slots.
    Plaintext shiftMask;
};

namespace detail {

using Clock = std::chrono::steady_clock;

inline double Milliseconds(const Clock::time_point& begin,
                           const Clock::time_point& end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

inline std::vector<double> MakeActiveMaskValues() {
    std::vector<double> values(kTotalSlots, 0.0);

    for (uint32_t word = 0; word < kNumWords; ++word) {
        const uint32_t base = word * kStride;

        for (uint32_t bit = 0; bit < kWordBits; ++bit) {
            values[base + bit] = 1.0;
        }
    }

    return values;
}

inline std::vector<double> MakeShiftMaskValues(uint32_t shift) {
    std::vector<double> values(kTotalSlots, 0.0);

    if (shift >= kWordBits) {
        return values;
    }

    for (uint32_t word = 0; word < kNumWords; ++word) {
        const uint32_t base = word * kStride;

        for (uint32_t bit = 0; bit < kWordBits - shift; ++bit) {
            values[base + bit] = 1.0;
        }
    }

    return values;
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
 * Builds the active-slot and public-right-shift masks and explicitly loads
 * them to the selected GPUs. Call this only after cc->LoadContext(...).
 */
inline Masks BuildAndLoadMasks(const CryptoContext<DCRTPoly>& cc,
                               uint32_t shift) {
    if (!cc) {
        throw std::invalid_argument(
            "BuildAndLoadMasks received a null context");
    }

    Masks masks;

    masks.activeMask = cc->MakeCKKSPackedPlaintext(
        detail::MakeActiveMaskValues(),
        1,
        0,
        nullptr,
        kTotalSlots);

    masks.shiftMask = cc->MakeCKKSPackedPlaintext(
        detail::MakeShiftMaskValues(shift),
        1,
        0,
        nullptr,
        kTotalSlots);

    masks.activeMask->SetLength(kTotalSlots);
    masks.shiftMask->SetLength(kTotalSlots);
detail::Synchronize(cc);

    return masks;
}

/**
 * Evaluates 128 independent Solidity-compatible uint256 logical right shifts
 * with one public shift amount.
 *
 * Layout:
 *   [256 Boolean data slots][256 zero guard slots] x 128.
 *
 * Bit order inside every word:
 *   base+0   = bit 0  (LSB),
 *   base+255 = bit 255 (MSB).
 *
 * Solidity uint256 logical-right-shift semantics:
 *
 *   result = input >> shift,  if shift < 256;
 *   result = 0,                           if shift >= 256.
 *
 * Optimized execution:
 *
 *   shift == 0:
 *     return the input ciphertext directly; no rotation and no mask multiply.
 *
 *   1 <= shift < 256:
 *     one CKKS rotation by +shift, followed by one plaintext mask multiply.
 *     The mask clears the cyclic wraparound, the high result bits, and
 *     all guard slots.
 *
 *   shift >= 256:
 *     create an encrypted zero as input-input.
 *
 * No ciphertext-ciphertext multiplication is used. For shifts 1...255,
 * one ciphertext-plaintext mask multiplication clears cyclic wraparound.
 * No bootstrap is performed or required by this operator.
 */
inline Ciphertext<DCRTPoly> EvalShiftRightPublic(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& input,
    uint32_t shift,
    Masks& masks,
    TimingReport* timing = nullptr) {

    if (!cc || !input) {
        throw std::invalid_argument(
            "EvalShiftRightPublic received a null context or ciphertext");
    }

    if (!masks.activeMask || !masks.shiftMask) {
        throw std::invalid_argument(
            "EvalShiftRightPublic received invalid masks");
    }

    TimingReport local;
    TimingReport& report = timing ? *timing : local;
    report = TimingReport{};

    detail::Synchronize(cc);
    const auto totalBegin = detail::Clock::now();

    if (shift >= kWordBits) {
        auto zero = detail::TimedGpuCall(
            cc,
            &report.zeroMs,
            [&]() {
                return cc->EvalSub(input, input);
            });

        report.totalMs =
            detail::Milliseconds(totalBegin, detail::Clock::now());

        return zero;
    }

    if (shift == 0U) {
        // The input layout contract already guarantees zero guard slots.
        report.totalMs =
            detail::Milliseconds(totalBegin, detail::Clock::now());

        return input;
    }

    auto rotated = detail::TimedGpuCall(
        cc,
        &report.rotateMs,
        [&]() {
            return cc->EvalRotate(
                input,
                static_cast<int32_t>(shift));
        });

    auto result = detail::TimedGpuCall(
        cc,
        &report.maskMs,
        [&]() {
            return cc->EvalMult(
                rotated,
                masks.shiftMask);
        });

    report.totalMs =
        detail::Milliseconds(totalBegin, detail::Clock::now());

    return result;
}

}  // namespace ks256gpu




