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

// The 40 GB A100 configuration that has already been verified to fit:
// 64 words x (256 active bits + 256 guard slots) = 32768 CKKS slots.
constexpr uint32_t kWordBits   = 256;
constexpr uint32_t kGuardBits  = 256;
constexpr uint32_t kStride     = kWordBits + kGuardBits;  // 512 slots/word
constexpr uint32_t kNumWords   = 128;
constexpr uint32_t kTotalSlots = kNumWords * kStride;     // 32768 slots
constexpr uint32_t kNumLayers  = 8;                       // log2(256)

static_assert(kTotalSlots <= 65536,
              "For ring dimension 131072, CKKS slot count cannot exceed 65536");
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

struct SubTimingReport {
    double preprocessMs   = 0.0;
    std::array<LayerTiming, kNumLayers> layers{};
    double borrowRotateMs = 0.0;
    double borrowMaskMs   = 0.0;
    double groupMaskMs    = 0.0;
    double reconstructMs  = 0.0;
    double totalMs        = 0.0;
};

struct Masks {
    Plaintext activeMask;  // 1 on the 256 data slots, 0 on guard slots
    Plaintext carryMask;   // 0 on bit 0 and guards; 1 on bits 1..255
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
            // Default layout: base+0 is bit 0 (LSB).
            // The carry/borrow into bit 0 must be zero.
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
 * Creates and explicitly loads the two plaintext masks on the selected CPU.
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
 * Evaluates kNumWords independent 256-bit additions in one CKKS ciphertext.
 *
 * After the prefix network, G[i] = c_{i+1}. Rotating G by one slot gives
 * C[i] = c_i. The final sum uses
 *
 *   s_i = a_i + b_i + c_i - 2c_{i+1}.
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

    // g = A AND B = AB
    // p = A XOR B = A + B - 2AB
    detail::Synchronize(cc);
    auto begin = detail::Clock::now();
    auto g = cc->EvalMult(cA, cB);
    auto p = cc->EvalSub(cc->EvalAdd(cA, cB), cc->EvalAdd(g, g));
    detail::Synchronize(cc);
    report.preprocessMs = detail::Milliseconds(begin, detail::Clock::now());

    auto G  = g;
    auto Pk = p;

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

    auto C = detail::TimedGpuCall(
        cc, &report.carryRotateMs,
        [&]() { return cc->EvalRotate(G, directionSign); });

    C = detail::TimedGpuCall(
        cc, &report.carryMaskMs,
        [&]() { return cc->EvalMult(C, masks.carryMask); });

    G = detail::TimedGpuCall(
        cc, &report.groupMaskMs,
        [&]() { return cc->EvalMult(G, masks.activeMask); });

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

/**
 * Evaluates kNumWords independent unsigned 256-bit subtractions using a
 * borrow Kogge-Stone prefix network.
 *
 * Per bit:
 *   borrow-generate  g_i = (~a_i) & b_i = b_i - a_i b_i
 *   borrow-propagate p_i = ~(a_i XOR b_i)
 *                         = 1 - a_i - b_i + 2a_i b_i
 *
 * Since the guard slots must have p=0, activeMask is used as the vector of
 * ones over active bits:
 *   p = activeMask - (A XOR B).
 *
 * The borrow recurrence is
 *   r_{i+1} = g_i OR (p_i AND r_i),
 * and the same eight-stage Kogge-Stone prefix operator is used.
 *
 * After the prefix network, G[i] = r_{i+1}. Rotating G by one slot gives
 * R[i] = r_i. The difference is reconstructed without another ciphertext
 * multiplication:
 *
 *   d_i = a_i - b_i - r_i + 2r_{i+1}.
 *
 * Therefore the core ciphertext multiplication depth is the same as the
 * optimized adder: one preprocessing multiplication plus eight prefix levels.
 *
 * If finalBorrowPrefix is non-null, it receives the active-slot-masked G.
 * The bit-255 slot of each word is the unsigned underflow flag.
 */
inline Ciphertext<DCRTPoly> EvalSubKS(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& cA,
    const Ciphertext<DCRTPoly>& cB,
    int directionSign,
    Masks& masks,
    SubTimingReport* timing = nullptr,
    Ciphertext<DCRTPoly>* finalBorrowPrefix = nullptr) {

    if (!cc || !cA || !cB) {
        throw std::invalid_argument(
            "EvalSubKS received a null context or ciphertext");
    }
    if (!masks.activeMask || !masks.carryMask) {
        throw std::invalid_argument("EvalSubKS received invalid masks");
    }
    if (directionSign != -1 && directionSign != 1) {
        throw std::invalid_argument("directionSign must be -1 or +1");
    }

    SubTimingReport local;
    SubTimingReport& report = timing ? *timing : local;
    report = SubTimingReport{};

    detail::Synchronize(cc);
    const auto totalBegin = detail::Clock::now();

    // 1. Borrow generate and propagate.
    //    ab  = A*B
    //    g   = B - AB                 = (~A)&B
    //    xor = A + B - 2AB
    //    p   = activeMask - xor       = ~(A XOR B) on active slots
    detail::Synchronize(cc);
    auto begin = detail::Clock::now();
    auto ab = cc->EvalMult(cA, cB);
    auto xorAB = cc->EvalSub(
        cc->EvalAdd(cA, cB), cc->EvalAdd(ab, ab));
    auto g = cc->EvalSub(cB, ab);
    auto p = cc->EvalSub(masks.activeMask, xorAB);
    detail::Synchronize(cc);
    report.preprocessMs = detail::Milliseconds(begin, detail::Clock::now());

    auto G  = g;
    auto Pk = p;

    // 2. Eight borrow-prefix stages: distances 1,2,...,128.
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

        // G and P*Grot are mutually exclusive Boolean states, so arithmetic
        // addition implements OR without an extra multiplication.
        G = detail::TimedGpuCall(
            cc, &report.layers[layer].updateGMs,
            [&]() { return cc->EvalAdd(Gold, term); });

        detail::Synchronize(cc);
        report.layers[layer].wallMs =
            detail::Milliseconds(layerBegin, detail::Clock::now());
    }

    // 3. R[i] = G[i-1]. Clear the borrow-in of each word's LSB and guards.
    auto R = detail::TimedGpuCall(
        cc, &report.borrowRotateMs,
        [&]() { return cc->EvalRotate(G, directionSign); });

    R = detail::TimedGpuCall(
        cc, &report.borrowMaskMs,
        [&]() { return cc->EvalMult(R, masks.carryMask); });

    // Keep only active G slots. G[base+255] is the final borrow/underflow flag.
    G = detail::TimedGpuCall(
        cc, &report.groupMaskMs,
        [&]() { return cc->EvalMult(G, masks.activeMask); });

    if (finalBorrowPrefix != nullptr) {
        *finalBorrowPrefix = G;
    }

    // 4. Exact full-subtractor reconstruction using additions/subtractions.
    //    D = A - B - R + 2G
    detail::Synchronize(cc);
    begin = detail::Clock::now();
    auto result = cc->EvalSub(cA, cB);
    result = cc->EvalSub(result, R);
    result = cc->EvalAdd(result, cc->EvalAdd(G, G));
    detail::Synchronize(cc);
    report.reconstructMs = detail::Milliseconds(begin, detail::Clock::now());

    report.totalMs = detail::Milliseconds(totalBegin, detail::Clock::now());
    return result;
}

}  // namespace ks256gpu




