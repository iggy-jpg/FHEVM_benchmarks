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

constexpr uint32_t kWordBits   = 16;
constexpr uint32_t kGuardBits  = 16;
constexpr uint32_t kStride     = kWordBits + kGuardBits;  // 32 slots/word
constexpr uint32_t kNumWords   = 2048;
constexpr uint32_t kTotalSlots = kNumWords * kStride;     // 65536 slots
constexpr uint32_t kNumLayers  = 4;                       // log2(16)

static_assert(kTotalSlots <= 65536,
              "This test is configured for a ring dimension of 131072, so the slot count must not exceed 65536");
static_assert((kTotalSlots & (kTotalSlots - 1)) == 0,
              "The CKKS slot count must be a power of two");

struct LayerTiming {
    double rotateGMs    = 0.0;
    double rotatePMs    = 0.0;
    double termMultMs   = 0.0;
    double prefixMultMs = 0.0;
    double updateGMs    = 0.0;
    double wallMs       = 0.0;
};

struct TimingReport {
    double preprocessMs  = 0.0;
    std::array<LayerTiming, kNumLayers> layers{};
    double carryRotateMs = 0.0;
    double carryMaskMs   = 0.0;
    double groupMaskMs   = 0.0;
    double reconstructMs = 0.0;
    double totalMs       = 0.0;
};

struct Masks {
    Plaintext activeMask;  // 1 on the 16 data slots, 0 on guard slots
    Plaintext carryMask;   // additionally clears the carry-in of bit 0
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

inline std::vector<double> MakeCarryMaskValues(int directionSign) {
    if (directionSign != -1 && directionSign != 1) {
        throw std::invalid_argument("directionSign must be -1 or +1");
    }

    std::vector<double> mask(kTotalSlots, 0.0);
    for (uint32_t word = 0; word < kNumWords; ++word) {
        const uint32_t base = word * kStride;

        if (directionSign == -1) {
            // Layout used by the test program:
            // base+0 is bit 0 (LSB), so C[base+0] must be zero.
            for (uint32_t bit = 1; bit < kWordBits; ++bit) {
                mask[base + bit] = 1.0;
            }
        }
        else {
            // Reverse layout: base+255 is bit 0.
            for (uint32_t bit = 0; bit + 1 < kWordBits; ++bit) {
                mask[base + bit] = 1.0;
            }
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

/**
 * Creates and explicitly loads the two masks on the selected CPU.
 * Call this only after cc->LoadContext(...).
 */
inline Masks BuildAndLoadMasks(const CryptoContext<DCRTPoly>& cc,
                               int directionSign) {
    if (!cc) {
        throw std::invalid_argument("BuildAndLoadMasks received a null context");
    }

    Masks masks;
    masks.activeMask = cc->MakeCKKSPackedPlaintext(
        detail::MakeActiveMaskValues(), 1, 0, nullptr, kTotalSlots);
    masks.carryMask = cc->MakeCKKSPackedPlaintext(
        detail::MakeCarryMaskValues(directionSign), 1, 0, nullptr, kTotalSlots);

    masks.activeMask->SetLength(kTotalSlots);
    masks.carryMask->SetLength(kTotalSlots);
detail::Synchronize(cc);
    return masks;
}

/**
 * Evaluates 2,048 independent 16-bit additions in one CKKS ciphertext.
 *
 * Default layout for directionSign = -1:
 *   [16 Boolean data slots][16 zero guard slots] x 2,048.
 *
 * After the prefix network, G[i] = c_{i+1}.  Rotating G by one slot gives
 * C[i] = c_i.  The final sum uses the full-adder identity
 *
 *   s_i = a_i + b_i + c_i - 2c_{i+1},
 *
 * so the final ciphertext-ciphertext multiplication required by
 * S = P XOR C is removed.
 *
 * Timing is synchronized at every reported CPU operation.  This produces a
 * detailed profile, but the synchronization overhead means that this function
 * should not be used as the final high-throughput benchmark without profiling.
 */
inline Ciphertext<DCRTPoly> EvalAddKS(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& cA,
    const Ciphertext<DCRTPoly>& cB,
    int directionSign,
    Masks& masks,
    TimingReport* timing = nullptr) {

    if (!cc || !cA || !cB) {
        throw std::invalid_argument(
            "EvalAddKS received a null context or ciphertext");
    }
    if (!masks.activeMask || !masks.carryMask) {
        throw std::invalid_argument("EvalAddKS received invalid masks");
    }
    if (directionSign != -1 && directionSign != 1) {
        throw std::invalid_argument("directionSign must be -1 or +1");
    }

    TimingReport local;
    TimingReport& report = timing ? *timing : local;
    report = TimingReport{};

    detail::Synchronize(cc);
    const auto totalBegin = detail::Clock::now();

    // 1. Bit generate and propagate:
    //    g = A AND B = AB
    //    p = A XOR B = A + B - 2AB
    detail::Synchronize(cc);
    auto begin = detail::Clock::now();
    auto g = cc->EvalMult(cA, cB);
    auto p = cc->EvalSub(cc->EvalAdd(cA, cB), cc->EvalAdd(g, g));
    detail::Synchronize(cc);
    report.preprocessMs = detail::Milliseconds(begin, detail::Clock::now());

    auto G  = g;
    auto Pk = p;

    // 2. Four Kogge-Stone stages: distances 1,2,4,8.
    for (uint32_t layer = 0; layer < kNumLayers; ++layer) {
        detail::Synchronize(cc);
        const auto layerBegin = detail::Clock::now();

        const int32_t step = static_cast<int32_t>(1U << layer);
        const int32_t rotation = directionSign * step;

        const auto Gold = G;
        const auto Pold = Pk;

        auto Grot = detail::TimedGpuCall(
            cc, &report.layers[layer].rotateGMs,
            [&]() { return cc->EvalRotate(Gold, rotation); });

        auto Prot = detail::TimedGpuCall(
            cc, &report.layers[layer].rotatePMs,
            [&]() { return cc->EvalRotate(Pold, rotation); });

        auto term = detail::TimedGpuCall(
            cc, &report.layers[layer].termMultMs,
            [&]() { return cc->EvalMult(Grot, Pold); });

        if (layer + 1 < kNumLayers) {
            Pk = detail::TimedGpuCall(
                cc, &report.layers[layer].prefixMultMs,
                [&]() { return cc->EvalMult(Prot, Pold); });
        }

        G = detail::TimedGpuCall(
            cc, &report.layers[layer].updateGMs,
            [&]() { return cc->EvalAdd(Gold, term); });

        detail::Synchronize(cc);
        report.layers[layer].wallMs =
            detail::Milliseconds(layerBegin, detail::Clock::now());
    }

    // 3. C[i] = G[i-1], then clear word boundaries and guard slots.
    auto C = detail::TimedGpuCall(
        cc, &report.carryRotateMs,
        [&]() { return cc->EvalRotate(G, directionSign); });

    C = detail::TimedGpuCall(
        cc, &report.carryMaskMs,
        [&]() { return cc->EvalMult(C, masks.carryMask); });

    // Keeping G zero in all guard slots prevents the reconstruction term -2G
    // from placing -2 in guard slots before bootstrapping.
    G = detail::TimedGpuCall(
        cc, &report.groupMaskMs,
        [&]() { return cc->EvalMult(G, masks.activeMask); });

    // 4. Exact full-adder reconstruction with additions/subtraction only.
    detail::Synchronize(cc);
    begin = detail::Clock::now();
    auto result = cc->EvalAdd(cA, cB);
    result = cc->EvalAdd(result, C);
    result = cc->EvalSub(result, cc->EvalAdd(G, G));
    detail::Synchronize(cc);
    report.reconstructMs = detail::Milliseconds(begin, detail::Clock::now());

    report.totalMs = detail::Milliseconds(totalBegin, detail::Clock::now());
    return result;
}

}  // namespace ks256gpu




