#pragma once

#include <openfhe.h>

#include <array>
#include <chrono>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace ks256gpu {

using lbcrypto::Ciphertext;
using lbcrypto::CryptoContext;
using lbcrypto::DCRTPoly;
using lbcrypto::Plaintext;

constexpr uint32_t kWordBits          = 256;
constexpr uint32_t kGuardBits         = 256;
constexpr uint32_t kStride            = kWordBits + kGuardBits;  // 512 slots/word
constexpr uint32_t kNumWords          = 128;
constexpr uint32_t kTotalSlots        = kNumWords * kStride;     // 65536 slots
constexpr uint32_t kShiftControlBits  = 8;                       // low 8 bits only
constexpr uint32_t kBarrelLayers      = 8;                       // shifts 1,2,...,128
constexpr uint32_t kPatternLayers     = 5;                       // repeat 8-bit pattern to 256 bits
constexpr uint32_t kBroadcastLayers   = 3;                       // fill each 8-slot block

static_assert(
    kTotalSlots <= 65536,
    "This test is configured for ring dimension 131072, so the slot count must not exceed 65536");
static_assert(
    (kTotalSlots & (kTotalSlots - 1U)) == 0U,
    "The CKKS slot count must be a power of two");
static_assert(
    (1U << kShiftControlBits) == kWordBits,
    "The encrypted barrel shifter expects a power-of-two word width");

struct PatternLayerTiming {
    double rotateMs = 0.0;
    double addMs    = 0.0;
    double wallMs   = 0.0;
};

struct BarrelLayerTiming {
    double controlLaneMaskMs = 0.0;
    double controlAnchorMs   = 0.0;
    std::array<double, kBroadcastLayers> controlBroadcastRotateMs{};
    std::array<double, kBroadcastLayers> controlBroadcastAddMs{};
    double controlBuildWallMs = 0.0;

    double dataRotateMs       = 0.0;
    double shiftMaskMs        = 0.0;
    double candidateDiffMs    = 0.0;
    double selectMultMs       = 0.0;
    double selectAddMs        = 0.0;
    double wallMs             = 0.0;
};

struct TimingReport {
    double lowByteMaskMs = 0.0;
    std::array<PatternLayerTiming, kPatternLayers> patternLayers{};
    std::array<BarrelLayerTiming, kBarrelLayers> barrelLayers{};
    double totalMs = 0.0;
};

struct Masks {
    // 1 on shift-amount bits 0...7 of every packed word.
    Plaintext lowByteMask;

    // For stage k, 1 at active positions whose bit index is congruent to k
    // modulo 8. These masks extract s_k from the repeated low-byte pattern.
    std::array<Plaintext, kShiftControlBits> controlLaneMasks{};

    // For stage k with distance d=2^k:
    //   1 on active result bits [0,255-d],
    //   0 on active result bits [256-d,255],
    //   0 on all guard slots.
    std::array<Plaintext, kBarrelLayers> shiftMasks{};
};

namespace detail {

using Clock = std::chrono::steady_clock;

inline double Milliseconds(const Clock::time_point& begin,
                           const Clock::time_point& end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

inline std::vector<double> MakeLowByteMaskValues() {
    std::vector<double> values(kTotalSlots, 0.0);

    for (uint32_t word = 0; word < kNumWords; ++word) {
        const uint32_t base = word * kStride;

        for (uint32_t bit = 0; bit < kShiftControlBits; ++bit) {
            values[base + bit] = 1.0;
        }
    }

    return values;
}

inline std::vector<double> MakeControlLaneMaskValues(uint32_t lane) {
    if (lane >= kShiftControlBits) {
        throw std::invalid_argument(
            "MakeControlLaneMaskValues received an invalid lane");
    }

    std::vector<double> values(kTotalSlots, 0.0);

    for (uint32_t word = 0; word < kNumWords; ++word) {
        const uint32_t base = word * kStride;

        for (uint32_t bit = lane; bit < kWordBits; bit += kShiftControlBits) {
            values[base + bit] = 1.0;
        }
    }

    return values;
}

inline std::vector<double> MakeShiftMaskValues(uint32_t distance) {
    if (distance == 0U || distance >= kWordBits) {
        throw std::invalid_argument(
            "MakeShiftMaskValues received an invalid distance");
    }

    std::vector<double> values(kTotalSlots, 0.0);

    for (uint32_t word = 0; word < kNumWords; ++word) {
        const uint32_t base = word * kStride;

        for (uint32_t bit = 0; bit < kWordBits - distance; ++bit) {
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
 * Builds and explicitly loads all plaintext masks used by the encrypted
 * uint256 barrel shifter. Call this only after cc->LoadContext(...).
 */
inline Masks BuildAndLoadMasks(const CryptoContext<DCRTPoly>& cc) {
    if (!cc) {
        throw std::invalid_argument(
            "BuildAndLoadMasks received a null context");
    }

    Masks masks;

    masks.lowByteMask = cc->MakeCKKSPackedPlaintext(
        detail::MakeLowByteMaskValues(),
        1,
        0,
        nullptr,
        kTotalSlots);

    masks.lowByteMask->SetLength(kTotalSlots);
for (uint32_t lane = 0; lane < kShiftControlBits; ++lane) {
        masks.controlLaneMasks[lane] = cc->MakeCKKSPackedPlaintext(
            detail::MakeControlLaneMaskValues(lane),
            1,
            0,
            nullptr,
            kTotalSlots);

        masks.controlLaneMasks[lane]->SetLength(kTotalSlots);
}

    for (uint32_t layer = 0; layer < kBarrelLayers; ++layer) {
        const uint32_t distance = 1U << layer;

        masks.shiftMasks[layer] = cc->MakeCKKSPackedPlaintext(
            detail::MakeShiftMaskValues(distance),
            1,
            0,
            nullptr,
            kTotalSlots);

        masks.shiftMasks[layer]->SetLength(kTotalSlots);
}

    detail::Synchronize(cc);
    return masks;
}

/**
 * Evaluates 128 independent uint256 logical right shifts with encrypted
 * uint256 shift amounts.
 *
 * Layout for both input ciphertexts:
 *   [256 Boolean data slots][256 zero guard slots] x 128.
 *
 * Bit order inside every word:
 *   base+0   = bit 0  (LSB),
 *   base+255 = bit 255 (MSB).
 *
 * Only the low eight bits of each encrypted shift amount are used:
 *
 *   effectiveShift = encryptedShift mod 256.
 *
 * Thus the operator implements the same modulo-word-width semantics as a
 * power-of-two barrel shifter. Shift-amount bits 8...255 remain encrypted but
 * do not participate in the circuit.
 *
 * Control preparation:
 *   1. Keep shift bits 0...7.
 *   2. Repeat that 8-bit pattern across the 256 active slots using five
 *      rotate-add doubling layers (distances 8,16,32,64,128).
 *   3. At barrel layer k, extract lane k and broadcast s_k within every
 *      8-slot block using rotations 1,2,4.
 *
 * Barrel layer k, distance d=2^k:
 *
 *   candidate = current >> d
 *   current   = current + s_k * (candidate - current)
 *
 * Each layer therefore uses exactly one ciphertext-ciphertext multiplication.
 * The eight serial layers have ciphertext multiplication depth 8. Plaintext
 * masks clear cyclic wraparound and all guard slots after every candidate
 * rotation. The caller performs one bootstrap after the complete circuit.
 */
inline Ciphertext<DCRTPoly> EvalShiftRightEncrypted(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& input,
    const Ciphertext<DCRTPoly>& encryptedShift,
    Masks& masks,
    TimingReport* timing = nullptr) {

    if (!cc || !input || !encryptedShift) {
        throw std::invalid_argument(
            "EvalShiftRightEncrypted received a null context or ciphertext");
    }

    if (!masks.lowByteMask) {
        throw std::invalid_argument(
            "EvalShiftRightEncrypted received invalid masks");
    }

    for (uint32_t layer = 0; layer < kBarrelLayers; ++layer) {
        if (!masks.controlLaneMasks[layer] || !masks.shiftMasks[layer]) {
            throw std::invalid_argument(
                "EvalShiftRightEncrypted received incomplete masks");
        }
    }

    TimingReport local;
    TimingReport& report = timing ? *timing : local;
    report = TimingReport{};

    detail::Synchronize(cc);
    const auto totalBegin = detail::Clock::now();

    // Keep only s_0...s_7 in every packed shift word.
    auto lowBytePattern = detail::TimedGpuCall(
        cc,
        &report.lowByteMaskMs,
        [&]() {
            return cc->EvalMult(
                encryptedShift,
                masks.lowByteMask);
        });

    // Repeat [s0,...,s7] 32 times across each 256-slot active region.
    for (uint32_t layer = 0; layer < kPatternLayers; ++layer) {
        auto& item = report.patternLayers[layer];

        detail::Synchronize(cc);
        const auto layerBegin = detail::Clock::now();

        const int32_t distance =
            -static_cast<int32_t>(8U << layer);

        const auto oldPattern = lowBytePattern;

        auto rotated = detail::TimedGpuCall(
            cc,
            &item.rotateMs,
            [&]() {
                return cc->EvalRotate(
                    oldPattern,
                    distance);
            });

        lowBytePattern = detail::TimedGpuCall(
            cc,
            &item.addMs,
            [&]() {
                return cc->EvalAdd(
                    oldPattern,
                    rotated);
            });

        item.wallMs =
            detail::Milliseconds(layerBegin, detail::Clock::now());
    }

    auto current = input;

    // Eight conditional shift stages controlled by s0...s7.
    for (uint32_t layer = 0; layer < kBarrelLayers; ++layer) {
        auto& item = report.barrelLayers[layer];

        detail::Synchronize(cc);
        const auto layerBegin = detail::Clock::now();
        const auto controlBegin = layerBegin;

        // Keep s_layer at positions whose bit index is layer modulo 8.
        auto control = detail::TimedGpuCall(
            cc,
            &item.controlLaneMaskMs,
            [&]() {
                return cc->EvalMult(
                    lowBytePattern,
                    masks.controlLaneMasks[layer]);
            });

        // Move each selected control to the first slot of its 8-slot block.
        if (layer != 0U) {
            control = detail::TimedGpuCall(
                cc,
                &item.controlAnchorMs,
                [&]() {
                    return cc->EvalRotate(
                        control,
                        static_cast<int32_t>(layer));
                });
        }

        // Broadcast the block anchor to all eight slots in each block.
        for (uint32_t broadcastLayer = 0;
             broadcastLayer < kBroadcastLayers;
             ++broadcastLayer) {
            const int32_t distance =
                -static_cast<int32_t>(1U << broadcastLayer);

            const auto oldControl = control;

            auto rotated = detail::TimedGpuCall(
                cc,
                &item.controlBroadcastRotateMs[broadcastLayer],
                [&]() {
                    return cc->EvalRotate(
                        oldControl,
                        distance);
                });

            control = detail::TimedGpuCall(
                cc,
                &item.controlBroadcastAddMs[broadcastLayer],
                [&]() {
                    return cc->EvalAdd(
                        oldControl,
                        rotated);
                });
        }

        item.controlBuildWallMs =
            detail::Milliseconds(controlBegin, detail::Clock::now());

        const uint32_t distance = 1U << layer;
        const auto oldCurrent = current;

        auto rotatedData = detail::TimedGpuCall(
            cc,
            &item.dataRotateMs,
            [&]() {
                return cc->EvalRotate(
                    oldCurrent,
                    static_cast<int32_t>(distance));
            });

        auto candidate = detail::TimedGpuCall(
            cc,
            &item.shiftMaskMs,
            [&]() {
                return cc->EvalMult(
                    rotatedData,
                    masks.shiftMasks[layer]);
            });

        auto difference = detail::TimedGpuCall(
            cc,
            &item.candidateDiffMs,
            [&]() {
                return cc->EvalSub(
                    candidate,
                    oldCurrent);
            });

        auto selectedDifference = detail::TimedGpuCall(
            cc,
            &item.selectMultMs,
            [&]() {
                return cc->EvalMult(
                    control,
                    difference);
            });

        current = detail::TimedGpuCall(
            cc,
            &item.selectAddMs,
            [&]() {
                return cc->EvalAdd(
                    oldCurrent,
                    selectedDifference);
            });

        item.wallMs =
            detail::Milliseconds(layerBegin, detail::Clock::now());
    }

    report.totalMs =
        detail::Milliseconds(totalBegin, detail::Clock::now());

    return current;
}

}  // namespace ks256gpu




