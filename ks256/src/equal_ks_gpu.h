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

constexpr uint32_t kWordBits   = 256;
constexpr uint32_t kGuardBits  = 256;
constexpr uint32_t kStride     = kWordBits + kGuardBits;  // 512 slots/word
constexpr uint32_t kNumWords   = 128;
constexpr uint32_t kTotalSlots = kNumWords * kStride;     // 65536 slots
constexpr uint32_t kNumLayers  = 8;                       // log2(256)

static_assert(
    kTotalSlots <= 65536,
    "This test is configured for ring dimension 131072, so the slot count must not exceed 65536");
static_assert(
    (kTotalSlots & (kTotalSlots - 1)) == 0,
    "The CKKS slot count must be a power of two");

// The public timing layout is intentionally retained so the existing CPP file
// remains source-compatible.  rotatePMs now covers the G and P rotations in
// one Kogge-Stone stage, while prefixMultMs covers the G/P stage update.
struct LayerTiming {
    double rotatePMs    = 0.0;
    double prefixMultMs = 0.0;
    double wallMs       = 0.0;
};

struct TimingReport {
    double preprocessMs   = 0.0;
    std::array<LayerTiming, kNumLayers> layers{};
    double resultRotateMs = 0.0;
    double resultMaskMs   = 0.0;
    double totalMs        = 0.0;
};

struct Masks {
    Plaintext activeMask;  // 1 on 256 data slots, 0 on guard slots
    Plaintext resultMask;  // 1 only at bit 0 of every packed word
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

/**
 * Creates and explicitly loads the equality masks on the selected CPU.
 * Call this only after cc->LoadContext(...).
 */
inline Masks BuildAndLoadMasks(const CryptoContext<DCRTPoly>& cc) {
    if (!cc) {
        throw std::invalid_argument(
            "BuildAndLoadMasks received a null context");
    }

    Masks masks;
    masks.activeMask = cc->MakeCKKSPackedPlaintext(
        detail::MakeActiveMaskValues(), 1, 0, nullptr, kTotalSlots);
    masks.resultMask = cc->MakeCKKSPackedPlaintext(
        detail::MakeResultMaskValues(), 1, 0, nullptr, kTotalSlots);

    masks.activeMask->SetLength(kTotalSlots);
    masks.resultMask->SetLength(kTotalSlots);
detail::Synchronize(cc);
    return masks;
}

/**
 * Evaluates 128 independent 256-bit equality operations.
 *
 * Layout:
 *   [256 Boolean data slots][256 zero guard slots] x 128.
 *
 * Instead of directly multiplying 256 XNOR values, this circuit propagates a
 * single "a difference exists" state through an eight-stage Kogge-Stone
 * network.
 *
 * Per bit:
 *
 *   g_i = XOR(a_i,b_i)
 *       = a_i + b_i - 2a_i b_i,
 *
 *   p_i = XNOR(a_i,b_i)
 *       = 1 - a_i - b_i + 2a_i b_i.
 *
 * For a higher bit/interval H and a lower interval L:
 *
 *   G(H:L) = G(H) OR (P(H) AND G(L)),
 *   P(H:L) = P(H) AND P(L).
 *
 * G(H) and P(H)G(L) are mutually exclusive Boolean states, so arithmetic
 * addition implements OR without an additional ciphertext multiplication:
 *
 *   G(H:L) = G(H) + P(H)G(L).
 *
 * After eight stages, bit 255 of each word is NE(A,B).  It is rotated to bit
 * 0, masked, and negated:
 *
 *   EQ(A,B) = 1 - NE(A,B).
 *
 * The caller performs exactly one bootstrap after the complete circuit.
 */
inline Ciphertext<DCRTPoly> EvalEqualKS(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& cA,
    const Ciphertext<DCRTPoly>& cB,
    Masks& masks,
    TimingReport* timing = nullptr) {

    if (!cc || !cA || !cB) {
        throw std::invalid_argument(
            "EvalEqualKS received a null context or ciphertext");
    }
    if (!masks.activeMask || !masks.resultMask) {
        throw std::invalid_argument(
            "EvalEqualKS received invalid masks");
    }

    TimingReport local;
    TimingReport& report = timing ? *timing : local;
    report = TimingReport{};

    detail::Synchronize(cc);
    const auto totalBegin = detail::Clock::now();

    // 1. Build the local difference-generate state G=XOR(A,B) and the
    // equality-propagate state P=XNOR(A,B).  AB is computed only once.
    detail::Synchronize(cc);
    auto begin = detail::Clock::now();

    auto ab = cc->EvalMult(cA, cB);
    auto xorAB = cc->EvalSub(
        cc->EvalAdd(cA, cB),
        cc->EvalAdd(ab, ab));

    auto G = xorAB;
    auto Pk = cc->EvalSub(masks.activeMask, xorAB);

    detail::Synchronize(cc);
    report.preprocessMs =
        detail::Milliseconds(begin, detail::Clock::now());

    // 2. Eight LSB-to-MSB Kogge-Stone stages.  With the existing rotation
    // keys, Rotate(...,-step) supplies the lower interval at each destination:
    //
    //   Gnew = Gold + Pold*Grot
    //   Pnew = Pold*Prot
    //
    // P is not required after the final stage.
    for (uint32_t layer = 0; layer < kNumLayers; ++layer) {
        detail::Synchronize(cc);
        const auto layerBegin = detail::Clock::now();

        const int32_t step = static_cast<int32_t>(1U << layer);
        const int32_t rotation = -step;

        const auto Gold = G;
        const auto Pold = Pk;

        detail::Synchronize(cc);
        begin = detail::Clock::now();

        auto Grot = cc->EvalRotate(Gold, rotation);
        auto Prot = cc->EvalRotate(Pold, rotation);

        detail::Synchronize(cc);
        report.layers[layer].rotatePMs =
            detail::Milliseconds(begin, detail::Clock::now());

        detail::Synchronize(cc);
        begin = detail::Clock::now();

        auto lowerDifference = cc->EvalMult(Pold, Grot);

        if (layer + 1 < kNumLayers) {
            Pk = cc->EvalMult(Pold, Prot);
        }

        G = cc->EvalAdd(Gold, lowerDifference);

        detail::Synchronize(cc);
        report.layers[layer].prefixMultMs =
            detail::Milliseconds(begin, detail::Clock::now());

        report.layers[layer].wallMs =
            detail::Milliseconds(layerBegin, detail::Clock::now());
    }

    // 3. G[base+255] is NE(A,B).  Rotate it to bit 0.
    auto result = detail::TimedGpuCall(
        cc,
        &report.resultRotateMs,
        [&]() {
            return cc->EvalRotate(
                G, static_cast<int32_t>(kWordBits - 1U));
        });

    // 4. Keep one result slot per word and compute EQ=1-NE.  Both operations
    // are grouped in resultMaskMs to preserve the original public timing API.
    result = detail::TimedGpuCall(
        cc,
        &report.resultMaskMs,
        [&]() {
            auto notEqual = cc->EvalMult(result, masks.resultMask);
            return cc->EvalSub(masks.resultMask, notEqual);
        });

    report.totalMs =
        detail::Milliseconds(totalBegin, detail::Clock::now());
    return result;
}

}  // namespace ks256gpu




