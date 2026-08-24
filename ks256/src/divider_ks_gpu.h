#pragma once

#include <openfhe.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ks256gpu {

using lbcrypto::Ciphertext;
using lbcrypto::CryptoContext;
using lbcrypto::DCRTPoly;
using lbcrypto::Plaintext;

// -----------------------------------------------------------------------------
// Packed Boolean layout
// -----------------------------------------------------------------------------
// Each uint256 division owns 1024 CKKS slots:
//   [0,255]   : Boolean uint256 input/output bits, LSB first
//   [256,511] : wide remainder headroom (exact radix-16 recurrence)
//   [512,1023]: guard slots
//
// The ciphertext remains full-packed at 65536 slots for N=131072.
constexpr uint32_t kDividerWordBits = 256;
constexpr uint32_t kDividerWideBits = 512;
constexpr uint32_t kDividerGuardBits = 512;
constexpr uint32_t kDividerStride = 1024;
constexpr uint32_t kDividerPackedWordCapacity = 64;
constexpr uint32_t kDividerNumWords = 64;
constexpr uint32_t kDividerTotalSlots =
    kDividerPackedWordCapacity * kDividerStride;

constexpr uint32_t kDividerKs256Layers = 8;
constexpr uint32_t kDividerKs512Layers = 9;
constexpr uint32_t kRadixBits = 4;
constexpr uint32_t kRadix = 1U << kRadixBits;
constexpr uint32_t kRadixRounds = kDividerWordBits / kRadixBits;
constexpr uint32_t kNonZeroCandidates = kRadix - 1U;

// Each radix-16 trial uses one KS subtraction to obtain both the trial
// difference and its borrow bit. The borrow selects either the trial
// difference or the old state through one mux.
//
// Fixed optimized path:
//   * G/P: serial Bootstrap + Boolean projection.
//   * borrow: Boolean projection only.
//   * state: Bootstrap + Boolean projection after every trial.
//   * quotient: no per-round Bootstrap; projected qBits are written once.
constexpr uint32_t kKsIntermediateBootstrapLayer = 3U;
constexpr uint32_t kCandidateKsDepth = 12;
constexpr uint32_t kCandidateMuxDepth = 1;
constexpr uint32_t kExactRadix16RoundDepth =
    kCandidateKsDepth + kCandidateMuxDepth;

static_assert(kDividerTotalSlots == 65536);
static_assert(kRadixRounds == 64);
static_assert(kDividerWideBits + kDividerGuardBits == kDividerStride);

enum class DividerDiagnosticLayout : uint8_t {
    BitZero,
    Active256,
    Active512,
};

using DividerDiagnosticSink = std::function<void(
    const std::string&,
    const Ciphertext<DCRTPoly>&,
    DividerDiagnosticLayout)>;

struct DividerOptions {
    uint32_t bootstrapInputMaxLevel = 39;
    bool conciseProgress = true;
    bool debugProgress = false;
    uint32_t debugMaxRounds = kRadixRounds;

    // Diagnostic mode executes one complete radix-16 round (8D,4D,2D,D),
    // checks every suspicious ciphertext, then returns before round two.
    bool diagnosticStopAfterFirstRound = false;
    DividerDiagnosticSink diagnosticSink{};
};

struct DividerTimingReport {
    uint64_t radixRounds = 0;
    uint64_t candidateSubtractions = 0;
    uint64_t ksAdditions = 0;
    uint64_t ksSubtractions = 0;
    uint64_t candidateMuxGates = 0;
    uint64_t ordinaryBootstraps = 0;
    uint64_t booleanProjections = 0;
    uint64_t roundStateBootstraps = 0;
    uint64_t multiplePreparationBootstraps = 0;
    uint64_t finalOutputBootstraps = 0;
    double totalMs = 0.0;
};

struct DividerResult {
    Ciphertext<DCRTPoly> quotient;
    Ciphertext<DCRTPoly> remainder;
};

struct DividerMasks {
    Plaintext bitZeroMask;
    Plaintext active256Mask;
    Plaintext active512Mask;
};

namespace divider_detail {

using Clock = std::chrono::steady_clock;

inline double Milliseconds(const Clock::time_point& begin,
                           const Clock::time_point& end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

inline void Synchronize(const CryptoContext<DCRTPoly>& cc) {
}

inline DividerDiagnosticLayout DiagnosticLayoutForWidth(uint32_t width) {
    if (width == kDividerWordBits) {
        return DividerDiagnosticLayout::Active256;
    }
    if (width == kDividerWideBits) {
        return DividerDiagnosticLayout::Active512;
    }
    throw std::invalid_argument("diagnostic width must be 256 or 512");
}

inline void EmitDiagnostic(
    const CryptoContext<DCRTPoly>& cc,
    const DividerOptions& options,
    const std::string& name,
    const Ciphertext<DCRTPoly>& value,
    DividerDiagnosticLayout layout) {

    if (!options.diagnosticSink || !value) {
        return;
    }

    Synchronize(cc);
    options.diagnosticSink(name, value, layout);
}

inline std::vector<double> MakeBitZeroMaskValues() {
    std::vector<double> values(kDividerTotalSlots, 0.0);
    for (uint32_t word = 0; word < kDividerNumWords; ++word) {
        values[word * kDividerStride] = 1.0;
    }
    return values;
}

inline std::vector<double> MakeActiveMaskValues(uint32_t width) {
    if (width != kDividerWordBits && width != kDividerWideBits) {
        throw std::invalid_argument("active width must be 256 or 512");
    }
    std::vector<double> values(kDividerTotalSlots, 0.0);
    for (uint32_t word = 0; word < kDividerNumWords; ++word) {
        const uint32_t base = word * kDividerStride;
        for (uint32_t bit = 0; bit < width; ++bit) {
            values[base + bit] = 1.0;
        }
    }
    return values;
}

inline Plaintext& ActiveMaskForWidth(uint32_t width, DividerMasks& masks) {
    if (width == kDividerWordBits) {
        return masks.active256Mask;
    }
    if (width == kDividerWideBits) {
        return masks.active512Mask;
    }
    throw std::invalid_argument("active width must be 256 or 512");
}

inline Ciphertext<DCRTPoly> RotateByPowers(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& input,
    int32_t shift) {

    if (!input || shift == 0) {
        return input;
    }

    Ciphertext<DCRTPoly> result = input;
    std::vector<Ciphertext<DCRTPoly>> keepAlive;
    keepAlive.reserve(10);

    uint32_t amount = static_cast<uint32_t>(shift < 0 ? -shift : shift);
    int32_t step = 1;
    while (amount != 0U) {
        if ((amount & 1U) != 0U) {
            auto next = cc->EvalRotate(result, shift < 0 ? -step : step);
            keepAlive.push_back(result);
            result = std::move(next);
        }
        amount >>= 1U;
        step <<= 1;
    }
    Synchronize(cc);
    return result;
}

inline Ciphertext<DCRTPoly> ShiftLeft(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& value,
    uint32_t amount) {
    return RotateByPowers(cc, value, -static_cast<int32_t>(amount));
}

inline Ciphertext<DCRTPoly> ShiftRight(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& value,
    uint32_t amount) {
    return RotateByPowers(cc, value, static_cast<int32_t>(amount));
}

inline Ciphertext<DCRTPoly> MaskWidth(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& value,
    uint32_t width,
    DividerMasks& masks) {
    auto masked = cc->EvalMult(value, ActiveMaskForWidth(width, masks));
    Synchronize(cc);
    return masked;
}

inline Ciphertext<DCRTPoly> ExtractBitToZero(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& value,
    uint32_t bit,
    DividerMasks& masks) {
    auto rotated = ShiftRight(cc, value, bit);
    auto extracted = cc->EvalMult(rotated, masks.bitZeroMask);
    Synchronize(cc);
    return extracted;
}

// Starting with one encrypted bit at slot zero of each 1024-slot word, fill
// slots [0,width-1] by an addition-only doubling network. No final plaintext
// mask is required: width is at most 512, so the expansion never reaches the
// neighbouring packed word.
inline Ciphertext<DCRTPoly> BroadcastBitZeroRaw(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& bitAtZero,
    uint32_t width) {

    const uint32_t layers = width == kDividerWordBits
                                ? kDividerKs256Layers
                                : (width == kDividerWideBits
                                       ? kDividerKs512Layers
                                       : 0U);
    if (layers == 0U) {
        throw std::invalid_argument("broadcast width must be 256 or 512");
    }

    auto result = bitAtZero;
    std::vector<Ciphertext<DCRTPoly>> keepAlive;
    keepAlive.reserve(2U * layers);
    for (uint32_t layer = 0; layer < layers; ++layer) {
        const int32_t step = -static_cast<int32_t>(1U << layer);
        auto rotated = cc->EvalRotate(result, step);
        auto next = cc->EvalAdd(result, rotated);
        keepAlive.push_back(result);
        keepAlive.push_back(std::move(rotated));
        result = std::move(next);
    }
    Synchronize(cc);
    return result;
}

inline Ciphertext<DCRTPoly> ClearWordBitZero(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& value,
    DividerMasks& masks) {
    auto bitZero = cc->EvalMult(value, masks.bitZeroMask);
    auto cleared = cc->EvalSub(value, bitZero);
    Synchronize(cc);
    return cleared;
}

inline Ciphertext<DCRTPoly> BootstrapRaw(
    const CryptoContext<DCRTPoly>& cc,
    Ciphertext<DCRTPoly> value,
    const DividerOptions& options,
    DividerTimingReport* timing) {

    if (!value) {
        throw std::invalid_argument("null ciphertext passed to bootstrap");
    }
    const uint32_t level = value->GetLevel();
    if (level > options.bootstrapInputMaxLevel) {
        throw std::runtime_error(
            "ordinary CKKS bootstrap requested too late: level=" +
            std::to_string(level) + ", latest=" +
            std::to_string(options.bootstrapInputMaxLevel));
    }
    
    
        return value;
    }
    auto refreshed = cc->EvalBootstrap(value, 2, 8);
    Synchronize(cc);  // keep the source alive until bootstrap completes
    if (timing) {
        ++timing->ordinaryBootstraps;
    }
    return refreshed;
}

inline Ciphertext<DCRTPoly> BootstrapResetWidth(
    const CryptoContext<DCRTPoly>& cc,
    Ciphertext<DCRTPoly> value,
    uint32_t width,
    DividerMasks& masks,
    const DividerOptions& options,
    DividerTimingReport* timing) {

    auto refreshed = BootstrapRaw(
        cc, std::move(value), options, timing);

    // Reusable packed state must have clean guard slots before rotations.
    // No cubic projection is applied: it adds ciphertext multiplications and
    // can increase the decryption-noise estimate.
    return MaskWidth(cc, refreshed, width, masks);
}

inline Ciphertext<DCRTPoly> BootstrapResetBitZero(
    const CryptoContext<DCRTPoly>& cc,
    Ciphertext<DCRTPoly> value,
    DividerMasks& masks,
    const DividerOptions& options,
    DividerTimingReport* timing) {

    auto refreshed = BootstrapRaw(
        cc, std::move(value), options, timing);
    auto bitZero = cc->EvalMult(refreshed, masks.bitZeroMask);
    Synchronize(cc);
    return bitZero;
}

inline Ciphertext<DCRTPoly> ProjectBooleanRaw(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& value,
    const DividerOptions& options,
    DividerTimingReport* timing,
    const std::string& label,
    DividerDiagnosticLayout layout) {

    if (!value) {
        throw std::invalid_argument(
            "null ciphertext passed to Boolean projection");
    }

    // Smoothstep Boolean projection:
    //     f(x) = 3x^2 - 2x^3
    // It fixes 0 and 1 and has zero derivative at both endpoints.
    auto x2 = cc->EvalMult(value, value);
    auto x3 = cc->EvalMult(x2, value);
    auto twoX2 = cc->EvalAdd(x2, x2);
    auto threeX2 = cc->EvalAdd(twoX2, x2);
    auto twoX3 = cc->EvalAdd(x3, x3);
    auto projected = cc->EvalSub(threeX2, twoX3);
    Synchronize(cc);

    if (timing) {
        ++timing->booleanProjections;
    }

    EmitDiagnostic(
        cc, options, label + ".x2", x2, layout);
    EmitDiagnostic(
        cc, options, label + ".x3", x3, layout);
    EmitDiagnostic(
        cc, options, label + ".output", projected, layout);

    return projected;
}

inline Ciphertext<DCRTPoly> ProjectBitZeroOnly(
    const CryptoContext<DCRTPoly>& cc,
    Ciphertext<DCRTPoly> value,
    DividerMasks& masks,
    const DividerOptions& options,
    DividerTimingReport* timing,
    const std::string& label) {

    if (!value) {
        throw std::invalid_argument(
            "null ciphertext passed to bit-only Boolean projection");
    }

    EmitDiagnostic(
        cc,
        options,
        label + ".input",
        value,
        DividerDiagnosticLayout::BitZero);

    auto projected = ProjectBooleanRaw(
        cc,
        value,
        options,
        timing,
        label + ".projection",
        DividerDiagnosticLayout::BitZero);

    auto bitZero = cc->EvalMult(
        projected,
        masks.bitZeroMask);
    Synchronize(cc);

    EmitDiagnostic(
        cc,
        options,
        label + ".masked",
        bitZero,
        DividerDiagnosticLayout::BitZero);

    return bitZero;
}

inline Ciphertext<DCRTPoly> BootstrapProjectVectorAlways(
    const CryptoContext<DCRTPoly>& cc,
    Ciphertext<DCRTPoly> value,
    uint32_t width,
    DividerMasks& masks,
    const DividerOptions& options,
    DividerTimingReport* timing,
    const std::string& label) {

    if (!value) {
        throw std::invalid_argument(
            "null ciphertext passed to projected vector bootstrap");
    }

    const uint32_t level = value->GetLevel();
    if (level > options.bootstrapInputMaxLevel) {
        throw std::runtime_error(
            "projected vector bootstrap requested too late: level=" +
            std::to_string(level) + ", latest=" +
            std::to_string(options.bootstrapInputMaxLevel));
    }

    const auto layout = DiagnosticLayoutForWidth(width);
    EmitDiagnostic(cc, options, label + ".input", value, layout);

    // Force a real single-iteration bootstrap.
    auto refreshed = cc->EvalBootstrap(value, 1, 0);
    Synchronize(cc);

    if (timing) {
        ++timing->ordinaryBootstraps;
    }

    EmitDiagnostic(
        cc, options, label + ".bootstrap_raw", refreshed, layout);

    auto projected = ProjectBooleanRaw(
        cc,
        refreshed,
        options,
        timing,
        label + ".projection",
        layout);

    auto masked = MaskWidth(
        cc,
        projected,
        width,
        masks);

    EmitDiagnostic(
        cc, options, label + ".masked", masked, layout);

    return masked;
}

struct KsResult {
    Ciphertext<DCRTPoly> value;
    Ciphertext<DCRTPoly> carryOrBorrowAtZero;
};

inline KsResult EvalAddKSWidth(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& a,
    const Ciphertext<DCRTPoly>& b,
    uint32_t width,
    DividerMasks& masks,
    DividerTimingReport* timing) {

    Plaintext* active = nullptr;
    uint32_t layers = 0;
    if (width == kDividerWordBits) {
        active = &masks.active256Mask;
        layers = kDividerKs256Layers;
    }
    else if (width == kDividerWideBits) {
        active = &masks.active512Mask;
        layers = kDividerKs512Layers;
    }
    else {
        throw std::invalid_argument("KS width must be 256 or 512");
    }

    auto ab = cc->EvalMult(a, b);
    auto aPlusB = cc->EvalAdd(a, b);
    auto twiceAB = cc->EvalAdd(ab, ab);
    auto P = cc->EvalSub(aPlusB, twiceAB);
    auto G = ab;

    std::vector<Ciphertext<DCRTPoly>> keepAlive;
    keepAlive.reserve(layers * 5U + 8U);
    for (uint32_t layer = 0; layer < layers; ++layer) {
        const int32_t rotation = -static_cast<int32_t>(1U << layer);
        const auto Gold = G;
        const auto Pold = P;
        auto Grot = cc->EvalRotate(Gold, rotation);
        auto Prot = cc->EvalRotate(Pold, rotation);
        auto term = cc->EvalMult(Grot, Pold);
        auto nextG = cc->EvalAdd(Gold, term);
        Ciphertext<DCRTPoly> nextP;
        if (layer + 1U < layers) {
            nextP = cc->EvalMult(Prot, Pold);
        }
        keepAlive.push_back(Gold);
        keepAlive.push_back(Pold);
        keepAlive.push_back(std::move(Grot));
        keepAlive.push_back(std::move(Prot));
        keepAlive.push_back(std::move(term));
        G = std::move(nextG);
        if (layer + 1U < layers) {
            P = std::move(nextP);
        }
    }

    auto carryInRaw = cc->EvalRotate(G, -1);
    auto carryIn = ClearWordBitZero(cc, carryInRaw, masks);
    auto maskedG = cc->EvalMult(G, *active);
    auto sum0 = cc->EvalAdd(a, b);
    auto sum1 = cc->EvalAdd(sum0, carryIn);
    auto twiceG = cc->EvalAdd(maskedG, maskedG);
    auto sum = cc->EvalSub(sum1, twiceG);
    auto carryOut = ExtractBitToZero(cc, maskedG, width - 1U, masks);
    Synchronize(cc);

    if (timing) {
        ++timing->ksAdditions;
    }
    return {std::move(sum), std::move(carryOut)};
}

inline KsResult EvalSubKSWidth(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& a,
    const Ciphertext<DCRTPoly>& b,
    uint32_t width,
    DividerMasks& masks,
    const DividerOptions& options,
    DividerTimingReport* timing,
    const std::string& diagnosticPrefix) {

    Plaintext* active = nullptr;
    uint32_t layers = 0;
    if (width == kDividerWordBits) {
        active = &masks.active256Mask;
        layers = kDividerKs256Layers;
    }
    else if (width == kDividerWideBits) {
        active = &masks.active512Mask;
        layers = kDividerKs512Layers;
    }
    else {
        throw std::invalid_argument("KS width must be 256 or 512");
    }

    const auto layout = DiagnosticLayoutForWidth(width);
    EmitDiagnostic(cc, options, diagnosticPrefix + ".input.a", a, layout);
    EmitDiagnostic(cc, options, diagnosticPrefix + ".input.b", b, layout);

    auto ab = cc->EvalMult(a, b);
    auto aPlusB = cc->EvalAdd(a, b);
    auto twiceAB = cc->EvalAdd(ab, ab);
    auto xorAB = cc->EvalSub(aPlusB, twiceAB);
    auto G = cc->EvalSub(b, ab);
    auto P = cc->EvalSub(*active, xorAB);
    Synchronize(cc);

    EmitDiagnostic(cc, options, diagnosticPrefix + ".init.ab", ab, layout);
    EmitDiagnostic(cc, options, diagnosticPrefix + ".init.xorAB", xorAB, layout);
    EmitDiagnostic(cc, options, diagnosticPrefix + ".init.G", G, layout);
    EmitDiagnostic(cc, options, diagnosticPrefix + ".init.P", P, layout);

    std::vector<Ciphertext<DCRTPoly>> keepAlive;
    keepAlive.reserve(layers * 5U + 8U);

    for (uint32_t layer = 0; layer < layers; ++layer) {
        const int32_t rotation = -static_cast<int32_t>(1U << layer);
        const auto Gold = G;
        const auto Pold = P;

        auto Grot = cc->EvalRotate(Gold, rotation);
        auto Prot = cc->EvalRotate(Pold, rotation);
        auto term = cc->EvalMult(Grot, Pold);
        auto nextG = cc->EvalAdd(Gold, term);

        Ciphertext<DCRTPoly> nextP;
        if (layer + 1U < layers) {
            nextP = cc->EvalMult(Prot, Pold);
        }

        keepAlive.push_back(Gold);
        keepAlive.push_back(Pold);
        keepAlive.push_back(std::move(Grot));
        keepAlive.push_back(std::move(Prot));
        keepAlive.push_back(std::move(term));

        G = std::move(nextG);
        if (layer + 1U < layers) {
            P = std::move(nextP);
        }

        if (width == kDividerWideBits &&
            layer == kKsIntermediateBootstrapLayer) {
            Synchronize(cc);

            EmitDiagnostic(
                cc, options, diagnosticPrefix + ".mid.pre.G", G, layout);
            EmitDiagnostic(
                cc, options, diagnosticPrefix + ".mid.pre.P", P, layout);

            G = BootstrapProjectVectorAlways(
                cc,
                std::move(G),
                width,
                masks,
                options,
                timing,
                diagnosticPrefix + ".mid.G");

            P = BootstrapProjectVectorAlways(
                cc,
                std::move(P),
                width,
                masks,
                options,
                timing,
                diagnosticPrefix + ".mid.P");

            keepAlive.clear();
        }
    }

    EmitDiagnostic(cc, options, diagnosticPrefix + ".final.G", G, layout);
    EmitDiagnostic(cc, options, diagnosticPrefix + ".final.P", P, layout);

    auto borrowInRaw = cc->EvalRotate(G, -1);
    auto borrowIn = ClearWordBitZero(cc, borrowInRaw, masks);
    auto maskedG = cc->EvalMult(G, *active);
    auto difference0 = cc->EvalSub(a, b);
    auto difference1 = cc->EvalSub(difference0, borrowIn);
    auto twiceG = cc->EvalAdd(maskedG, maskedG);
    auto difference = cc->EvalAdd(difference1, twiceG);
    auto borrowOut = ExtractBitToZero(cc, maskedG, width - 1U, masks);
    Synchronize(cc);

    EmitDiagnostic(
        cc, options, diagnosticPrefix + ".reconstruct.borrowInRaw",
        borrowInRaw, layout);
    EmitDiagnostic(
        cc, options, diagnosticPrefix + ".reconstruct.borrowIn",
        borrowIn, layout);
    EmitDiagnostic(
        cc, options, diagnosticPrefix + ".reconstruct.maskedG",
        maskedG, layout);
    EmitDiagnostic(
        cc, options, diagnosticPrefix + ".reconstruct.difference0",
        difference0, layout);
    EmitDiagnostic(
        cc, options, diagnosticPrefix + ".reconstruct.difference1",
        difference1, layout);
    EmitDiagnostic(
        cc, options, diagnosticPrefix + ".reconstruct.twiceG",
        twiceG, layout);
    EmitDiagnostic(
        cc, options, diagnosticPrefix + ".output.difference",
        difference, layout);
    EmitDiagnostic(
        cc, options, diagnosticPrefix + ".output.borrow",
        borrowOut, DividerDiagnosticLayout::BitZero);

    if (timing) {
        ++timing->ksSubtractions;
    }

    return {std::move(difference), std::move(borrowOut)};
}

inline Ciphertext<DCRTPoly> EvalSuffixOr256(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& input) {

    auto suffix = input;
    std::vector<Ciphertext<DCRTPoly>> keepAlive;
    keepAlive.reserve(kDividerKs256Layers * 4U);
    for (uint32_t layer = 0; layer < kDividerKs256Layers; ++layer) {
        const int32_t step = static_cast<int32_t>(1U << layer);
        auto old = suffix;
        auto shifted = cc->EvalRotate(old, step);
        auto both = cc->EvalMult(old, shifted);
        auto sum = cc->EvalAdd(old, shifted);
        auto next = cc->EvalSub(sum, both);
        keepAlive.push_back(std::move(old));
        keepAlive.push_back(std::move(shifted));
        keepAlive.push_back(std::move(both));
        keepAlive.push_back(std::move(sum));
        suffix = std::move(next);
    }
    Synchronize(cc);
    return suffix;
}

struct SafeDivisor {
    Ciphertext<DCRTPoly> valueWide;
    Ciphertext<DCRTPoly> divisorWasZeroAtZero;
};

inline SafeDivisor BuildSafeDivisor(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& divisor,
    DividerMasks& masks,
    const DividerOptions& options,
    DividerTimingReport* timing) {

    auto suffix = EvalSuffixOr256(cc, divisor);
    auto anyBit = cc->EvalMult(suffix, masks.bitZeroMask);
    auto oneAtZero = cc->EvalAdd(cc->EvalSub(divisor, divisor),
                                 masks.bitZeroMask);
    auto zeroBit = cc->EvalSub(oneAtZero, anyBit);
    Synchronize(cc);

    // The zero flag is reused throughout the division, so reset its modulus
    // chain once.  Ordinary bootstrap precision is sufficient because the
    // input remains on the correct side of the Boolean threshold.
    zeroBit = BootstrapResetBitZero(
        cc, std::move(zeroBit), masks, options, timing);

    // D=0 becomes the 257-bit value 2^256.
    // For every 256-bit dividend A, this naturally gives q=0 and r=A.
    auto zeroAtBit256 = ShiftLeft(cc, zeroBit, kDividerWordBits);
    auto safeWide = cc->EvalAdd(divisor, zeroAtBit256);
    Synchronize(cc);
    return {std::move(safeWide), std::move(zeroBit)};
}

// Build only the binary-weighted divisor multiples required by the
// sequential radix-16 digit extractor: D, 2D, 4D and 8D.
inline std::vector<Ciphertext<DCRTPoly>> BuildRadix16Multiples(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& safeDivisorWide,
    DividerMasks& masks,
    const DividerOptions& options,
    DividerTimingReport* timing) {

    // The sequential radix-16 digit extractor only needs the binary-weighted
    // divisor multiples D, 2D, 4D and 8D.
    std::vector<Ciphertext<DCRTPoly>> m(kRadix);
    m[0] = cc->EvalSub(safeDivisorWide, safeDivisorWide);

    m[1] = BootstrapResetWidth(
        cc, safeDivisorWide, kDividerWideBits, masks, options, timing);
    if (timing) {
        ++timing->multiplePreparationBootstraps;
    }

    m[2] = ShiftLeft(cc, m[1], 1U);
    m[4] = ShiftLeft(cc, m[1], 2U);
    m[8] = ShiftLeft(cc, m[1], 3U);
    Synchronize(cc);

    return m;
}

inline Ciphertext<DCRTPoly> BuildInputNibbleAtLowBits(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& dividend,
    uint32_t lowSourceBit,
    DividerMasks& masks) {

    auto nibble = cc->EvalSub(dividend, dividend);
    for (uint32_t bit = 0; bit < kRadixBits; ++bit) {
        auto atZero = ExtractBitToZero(
            cc, dividend, lowSourceBit + bit, masks);
        auto placed = bit == 0U ? atZero : ShiftLeft(cc, atZero, bit);
        auto next = cc->EvalAdd(nibble, placed);
        Synchronize(cc);
        nibble = std::move(next);
    }
    return nibble;
}

struct RadixRoundResult {
    Ciphertext<DCRTPoly> remainder;
    std::array<Ciphertext<DCRTPoly>, kRadixBits> quotientDigitBitsAtZero;
};

inline Ciphertext<DCRTPoly> MuxWideByBitAtZero(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& low,
    const Ciphertext<DCRTPoly>& high,
    const Ciphertext<DCRTPoly>& selectorAtZero,
    const DividerOptions& options,
    DividerTimingReport* timing,
    const std::string& diagnosticPrefix) {

    EmitDiagnostic(
        cc, options, diagnosticPrefix + ".low",
        low, DividerDiagnosticLayout::Active512);
    EmitDiagnostic(
        cc, options, diagnosticPrefix + ".high",
        high, DividerDiagnosticLayout::Active512);
    EmitDiagnostic(
        cc, options, diagnosticPrefix + ".selector",
        selectorAtZero, DividerDiagnosticLayout::BitZero);

    auto selectorWide = BroadcastBitZeroRaw(
        cc, selectorAtZero, kDividerWideBits);

    EmitDiagnostic(
        cc, options, diagnosticPrefix + ".selectorWide",
        selectorWide, DividerDiagnosticLayout::Active512);

    auto delta = cc->EvalSub(high, low);
    auto selectedDelta = cc->EvalMult(delta, selectorWide);
    auto result = cc->EvalAdd(low, selectedDelta);
    Synchronize(cc);

    EmitDiagnostic(
        cc, options, diagnosticPrefix + ".delta",
        delta, DividerDiagnosticLayout::Active512);
    EmitDiagnostic(
        cc, options, diagnosticPrefix + ".selectedDelta",
        selectedDelta, DividerDiagnosticLayout::Active512);
    EmitDiagnostic(
        cc, options, diagnosticPrefix + ".output",
        result, DividerDiagnosticLayout::Active512);

    if (timing) {
        ++timing->candidateMuxGates;
    }

    return result;
}

// Exact radix-16 round:
//   1. Form T = 16R + the next dividend nibble.
//   2. Test 8D, 4D, 2D and D from high bit to low bit.
//   3. One KS subtraction yields both the candidate difference and borrow.
//   4. A projected borrow drives a mux between the candidate and old state.
//   5. The selected state is bootstrapped and Boolean-projected before reuse.
inline RadixRoundResult EvalExactRadix16Round(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& currentRemainder,
    const Ciphertext<DCRTPoly>& dividend,
    uint32_t round,
    const std::vector<Ciphertext<DCRTPoly>>& multiples,
    DividerMasks& masks,
    const DividerOptions& options,
    DividerTimingReport* timing) {

    const uint32_t sourceLowBit =
        kDividerWordBits - kRadixBits * (round + 1U);

    auto shiftedR = ShiftLeft(cc, currentRemainder, kRadixBits);
    auto nibble = BuildInputNibbleAtLowBits(
        cc, dividend, sourceLowBit, masks);
    auto state = cc->EvalAdd(shiftedR, nibble);
    Synchronize(cc);

    EmitDiagnostic(
        cc, options, "round.input.currentRemainder",
        currentRemainder, DividerDiagnosticLayout::Active512);
    EmitDiagnostic(
        cc, options, "round.input.shiftedR",
        shiftedR, DividerDiagnosticLayout::Active512);
    EmitDiagnostic(
        cc, options, "round.input.nibble",
        nibble, DividerDiagnosticLayout::Active512);
    EmitDiagnostic(
        cc, options, "round.input.T",
        state, DividerDiagnosticLayout::Active512);

    std::array<Ciphertext<DCRTPoly>, kRadixBits> qBits{};

    for (int32_t bit = static_cast<int32_t>(kRadixBits) - 1;
         bit >= 0;
         --bit) {
        const uint32_t bitIndex = static_cast<uint32_t>(bit);
        const uint32_t weight = 1U << bitIndex;
        const std::string trialPrefix =
            "trial." + std::to_string(weight) + "D";

        if (!multiples[weight]) {
            throw std::runtime_error(
                "missing radix-16 divisor multiple " +
                std::to_string(weight) + "D");
        }

        EmitDiagnostic(
            cc, options, trialPrefix + ".input.state",
            state, DividerDiagnosticLayout::Active512);
        EmitDiagnostic(
            cc, options, trialPrefix + ".input.multiple",
            multiples[weight], DividerDiagnosticLayout::Active512);

        // One KS subtraction yields both:
        //   trial.value                 = state - weight*D
        //   trial.carryOrBorrowAtZero   = [state < weight*D]
        auto trial = EvalSubKSWidth(
            cc,
            state,
            multiples[weight],
            kDividerWideBits,
            masks,
            options,
            timing,
            trialPrefix + ".compareKS");

        if (timing) {
            ++timing->candidateSubtractions;
        }

        EmitDiagnostic(
            cc, options, trialPrefix + ".compareKS.savedDifference",
            trial.value, DividerDiagnosticLayout::Active512);

        auto borrow = std::move(trial.carryOrBorrowAtZero);

        // The KS borrow is already precise enough; project it directly
        // instead of introducing another Bootstrap approximation.
        borrow = ProjectBitZeroOnly(
            cc,
            std::move(borrow),
            masks,
            options,
            timing,
            trialPrefix + ".borrow.projectionOnly");

        // qBit = 1 - borrow.
        auto borrowZero = cc->EvalSub(borrow, borrow);
        auto oneAtZero = cc->EvalAdd(
            borrowZero,
            masks.bitZeroMask);
        auto qBit = cc->EvalSub(
            oneAtZero,
            borrow);
        Synchronize(cc);

        EmitDiagnostic(
            cc, options, trialPrefix + ".qBit.borrowZero",
            borrowZero, DividerDiagnosticLayout::BitZero);
        EmitDiagnostic(
            cc, options, trialPrefix + ".qBit.oneAtZero",
            oneAtZero, DividerDiagnosticLayout::BitZero);
        EmitDiagnostic(
            cc, options, trialPrefix + ".qBit.value",
            qBit, DividerDiagnosticLayout::BitZero);

        // borrow=0: accept trial.value
        // borrow=1: retain the old state
        auto selectedState = MuxWideByBitAtZero(
            cc,
            trial.value,
            state,
            borrow,
            options,
            timing,
            trialPrefix + ".mux");

        // Restore the modulus chain after every trial. Removing this Bootstrap
        // lets the next KS reach level 40 before its internal G/P refresh.
        state = BootstrapProjectVectorAlways(
            cc,
            std::move(selectedState),
            kDividerWideBits,
            masks,
            options,
            timing,
            trialPrefix + ".state.bootstrap");

        if (timing) {
            ++timing->roundStateBootstraps;
        }

        EmitDiagnostic(
            cc, options, trialPrefix + ".state.final",
            state, DividerDiagnosticLayout::Active512);

        qBits[bitIndex] = qBit;
    }

    return {std::move(state), std::move(qBits)};
}

inline Ciphertext<DCRTPoly> PlaceRadixDigitInQuotient(
    const CryptoContext<DCRTPoly>& cc,
    Ciphertext<DCRTPoly> quotient,
    const std::array<Ciphertext<DCRTPoly>, kRadixBits>& digitBitsAtZero,
    uint32_t round) {

    const uint32_t destinationLowBit =
        kDividerWordBits - kRadixBits * (round + 1U);
    for (uint32_t bit = 0; bit < kRadixBits; ++bit) {
        const uint32_t destination = destinationLowBit + bit;
        auto placed = destination == 0U
                          ? digitBitsAtZero[bit]
                          : ShiftLeft(cc, digitBitsAtZero[bit], destination);
        auto next = cc->EvalAdd(quotient, placed);
        Synchronize(cc);
        quotient = std::move(next);
    }
    return quotient;
}

inline DividerResult ApplyDivisionByZeroSemantics(
    const CryptoContext<DCRTPoly>& cc,
    Ciphertext<DCRTPoly> quotient,
    Ciphertext<DCRTPoly> remainder,
    const Ciphertext<DCRTPoly>& dividend,
    const Ciphertext<DCRTPoly>& divisorWasZeroAtZero,
    DividerMasks& masks,
    const DividerOptions& options,
    DividerTimingReport* timing) {

    auto zero256 = BroadcastBitZeroRaw(
        cc, divisorWasZeroAtZero, kDividerWordBits);

    // q = 0 when D=0.
    auto zeroTimesQ = cc->EvalMult(quotient, zero256);
    auto selectedQ = cc->EvalSub(quotient, zeroTimesQ);

    // r = A when D=0, otherwise use the computed remainder.
    auto remainder256 = MaskWidth(
        cc, remainder, kDividerWordBits, masks);
    auto aMinusR = cc->EvalSub(dividend, remainder256);
    auto zeroCorrection = cc->EvalMult(zero256, aMinusR);
    auto selectedR = cc->EvalAdd(remainder256, zeroCorrection);
    Synchronize(cc);

    selectedQ = BootstrapRaw(
        cc, std::move(selectedQ), options, timing);
    selectedR = BootstrapRaw(
        cc, std::move(selectedR), options, timing);
    if (timing) {
        timing->finalOutputBootstraps += 2U;
    }
    return {std::move(selectedQ), std::move(selectedR)};
}

} // namespace divider_detail

inline DividerMasks BuildDividerMasksOnCpu(
    const CryptoContext<DCRTPoly>& cc) {
    if (!cc) {
        throw std::invalid_argument("BuildDividerMasksOnCpu: null context");
    }

    DividerMasks masks;
    masks.bitZeroMask = cc->MakeCKKSPackedPlaintext(
        divider_detail::MakeBitZeroMaskValues(), 1, 0, nullptr,
        kDividerTotalSlots);
    masks.active256Mask = cc->MakeCKKSPackedPlaintext(
        divider_detail::MakeActiveMaskValues(kDividerWordBits),
        1, 0, nullptr, kDividerTotalSlots);
    masks.active512Mask = cc->MakeCKKSPackedPlaintext(
        divider_detail::MakeActiveMaskValues(kDividerWideBits),
        1, 0, nullptr, kDividerTotalSlots);

    masks.bitZeroMask->SetLength(kDividerTotalSlots);
    masks.active256Mask->SetLength(kDividerTotalSlots);
    masks.active512Mask->SetLength(kDividerTotalSlots);
    return masks;
}

inline void LoadDividerMasksToGpu(
    const CryptoContext<DCRTPoly>& cc,
    DividerMasks& masks) {
divider_detail::Synchronize(cc);
}

inline DividerResult EvalDivRadix16KS(
    const CryptoContext<DCRTPoly>& cc,
    const Ciphertext<DCRTPoly>& dividend,
    const Ciphertext<DCRTPoly>& divisor,
    DividerMasks& masks,
    const DividerOptions& options = DividerOptions{},
    DividerTimingReport* timing = nullptr) {

    if (!cc || !dividend || !divisor) {
        throw std::invalid_argument("EvalDivRadix16KS received a null input");
    }

    DividerTimingReport local;
    DividerTimingReport& report = timing ? *timing : local;
    report = DividerTimingReport{};
    const auto totalBegin = divider_detail::Clock::now();

    if (options.conciseProgress) {
        std::cout
            << "[divider] preparing encrypted D, 2D, 4D and 8D multiples.\n";
    }
    divider_detail::EmitDiagnostic(
        cc, options, "top.input.dividend",
        dividend, DividerDiagnosticLayout::Active256);
    divider_detail::EmitDiagnostic(
        cc, options, "top.input.divisor",
        divisor, DividerDiagnosticLayout::Active256);

    auto safe = divider_detail::BuildSafeDivisor(
        cc, divisor, masks, options, &report);

    divider_detail::EmitDiagnostic(
        cc, options, "top.safeDivisor",
        safe.valueWide, DividerDiagnosticLayout::Active512);
    divider_detail::EmitDiagnostic(
        cc, options, "top.divisorWasZero",
        safe.divisorWasZeroAtZero, DividerDiagnosticLayout::BitZero);

    auto multiples = divider_detail::BuildRadix16Multiples(
        cc, safe.valueWide, masks, options, &report);

    divider_detail::EmitDiagnostic(
        cc, options, "top.multiple.D",
        multiples[1], DividerDiagnosticLayout::Active512);
    divider_detail::EmitDiagnostic(
        cc, options, "top.multiple.2D",
        multiples[2], DividerDiagnosticLayout::Active512);
    divider_detail::EmitDiagnostic(
        cc, options, "top.multiple.4D",
        multiples[4], DividerDiagnosticLayout::Active512);
    divider_detail::EmitDiagnostic(
        cc, options, "top.multiple.8D",
        multiples[8], DividerDiagnosticLayout::Active512);

    auto remainder = cc->EvalSub(safe.valueWide, safe.valueWide);
    auto quotient = cc->EvalSub(dividend, dividend);
    divider_detail::Synchronize(cc);

    if (options.conciseProgress) {
        std::cout << "[divider] executing " << kRadixRounds
                  << " exact radix-16 rounds "
                  << "(8D -> 4D -> 2D -> D sequential trials).\n";
    }

    const uint32_t roundsToRun =
        std::min(options.debugMaxRounds, kRadixRounds);

    for (uint32_t round = 0; round < roundsToRun; ++round) {
        if (options.conciseProgress &&
            (round == 0U || ((round + 1U) % 4U) == 0U)) {
            std::cout << "[divider][radix16] round "
                      << (round + 1U)
                      << "/" << kRadixRounds
                      << " begin.\n";
        }

        auto roundResult = divider_detail::EvalExactRadix16Round(
            cc,
            remainder,
            dividend,
            round,
            multiples,
            masks,
            options,
            &report);

        divider_detail::EmitDiagnostic(
            cc, options, "round.output.remainder",
            roundResult.remainder,
            DividerDiagnosticLayout::Active512);

        for (uint32_t bit = 0; bit < kRadixBits; ++bit) {
            divider_detail::EmitDiagnostic(
                cc,
                options,
                "round.output.qBit" + std::to_string(bit),
                roundResult.quotientDigitBitsAtZero[bit],
                DividerDiagnosticLayout::BitZero);
        }

        quotient = divider_detail::PlaceRadixDigitInQuotient(
            cc,
            std::move(quotient),
            roundResult.quotientDigitBitsAtZero,
            round);

        divider_detail::EmitDiagnostic(
            cc, options, "round.quotient.placed",
            quotient, DividerDiagnosticLayout::Active256);

        // Every qBit has already been Boolean-projected and is written to a
        // disjoint quotient position, so the accumulated quotient does not
        // require a per-round Bootstrap.
        divider_detail::Synchronize(cc);

        divider_detail::EmitDiagnostic(
            cc, options, "round.quotient.final",
            quotient, DividerDiagnosticLayout::Active256);

        remainder = std::move(roundResult.remainder);

        ++report.radixRounds;

        if (options.diagnosticStopAfterFirstRound) {
            report.totalMs = divider_detail::Milliseconds(
                totalBegin, divider_detail::Clock::now());

            return {
                std::move(quotient),
                std::move(remainder)
            };
        }
    }

    // The safe-divisor construction already defines the division-by-zero
    // behavior. Optimized quotient and remainder outputs are decoded directly.
    divider_detail::Synchronize(cc);

    report.totalMs = divider_detail::Milliseconds(
        totalBegin, divider_detail::Clock::now());

    return {std::move(quotient), std::move(remainder)};
}

} // namespace ks256gpu




