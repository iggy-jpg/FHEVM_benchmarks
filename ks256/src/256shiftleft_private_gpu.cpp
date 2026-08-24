#include <openfhe.h>

#include "shiftleft_private_gpu.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace lbcrypto;

namespace {

using Clock = std::chrono::steady_clock;

struct PhaseTimes {
    double contextMs         = 0.0;
    double keyGenMs          = 0.0;
    double multKeyGenMs      = 0.0;
    double rotateKeyGenMs    = 0.0;
    double bootstrapSetupMs  = 0.0;
    double bootstrapKeyGenMs = 0.0;
    double loadContextMs     = 0.0;
    double maskBuildLoadMs   = 0.0;
    double inputBuildMs      = 0.0;
    double encodingMs        = 0.0;
    double encryptionMs      = 0.0;
    double shiftMs           = 0.0;
    double bootstrapMs       = 0.0;
    double decryptionMs      = 0.0;
    double verificationMs    = 0.0;
    double totalMs           = 0.0;
};

inline double Milliseconds(const Clock::time_point& begin,
                           const Clock::time_point& end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

inline double PrecisionBits(double error) {
    if (error <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }

    return -std::log2(error);
}

uint32_t LevelsRemaining(uint32_t configuredDepth,
                         const Ciphertext<DCRTPoly>& ciphertext) {
    if (!ciphertext) {
        return 0;
    }

    const uint32_t level = ciphertext->GetLevel();

    const uint32_t scalePenalty =
        ciphertext->GetNoiseScaleDeg() > 0
            ? ciphertext->GetNoiseScaleDeg() - 1
            : 0;

    if (configuredDepth <= level + scalePenalty) {
        return 0;
    }

    return configuredDepth - level - scalePenalty;
}

std::string BitsToHex(const std::vector<uint8_t>& bits,
                      uint32_t wordIndex) {
    constexpr char kHex[] = "0123456789abcdef";

    std::string output(
        ks256gpu::kWordBits / 4,
        '0');

    const uint32_t base =
        wordIndex * ks256gpu::kStride;

    for (uint32_t nibble = 0;
         nibble < ks256gpu::kWordBits / 4;
         ++nibble) {
        uint8_t value = 0;

        for (uint32_t offset = 0;
             offset < 4;
             ++offset) {
            const uint32_t bit =
                nibble * 4 + offset;

            value |= static_cast<uint8_t>(
                bits[base + bit] << offset);
        }

        output[output.size() - 1U - nibble] =
            kHex[value];
    }

    return "0x" + output;
}

uint32_t EffectiveShift(const std::vector<uint8_t>& shiftBits,
                        uint32_t wordIndex) {
    const uint32_t base =
        wordIndex * ks256gpu::kStride;

    uint32_t shift = 0;

    for (uint32_t bit = 0;
         bit < ks256gpu::kShiftControlBits;
         ++bit) {
        shift |=
            static_cast<uint32_t>(shiftBits[base + bit])
            << bit;
    }

    return shift;
}

void SetShiftWord(std::vector<uint8_t>* shiftBits,
                  std::vector<double>* shiftSlots,
                  uint32_t wordIndex,
                  uint32_t lowValue,
                  const std::vector<uint32_t>& highOneBits = {}) {
    const uint32_t base =
        wordIndex * ks256gpu::kStride;

    for (uint32_t bit = 0;
         bit < ks256gpu::kWordBits;
         ++bit) {
        (*shiftBits)[base + bit] = 0;
        (*shiftSlots)[base + bit] = 0.0;
    }

    for (uint32_t bit = 0;
         bit < ks256gpu::kShiftControlBits;
         ++bit) {
        const uint8_t value =
            static_cast<uint8_t>(
                (lowValue >> bit) & 1U);

        (*shiftBits)[base + bit] = value;
        (*shiftSlots)[base + bit] =
            static_cast<double>(value);
    }

    for (uint32_t bit : highOneBits) {
        if (bit < ks256gpu::kShiftControlBits ||
            bit >= ks256gpu::kWordBits) {
            throw std::invalid_argument(
                "SetShiftWord received an invalid high shift bit");
        }

        (*shiftBits)[base + bit] = 1;
        (*shiftSlots)[base + bit] = 1.0;
    }
}

void BuildExpectedShiftLeft(
    const std::vector<uint8_t>& input,
    const std::vector<uint8_t>& encryptedShift,
    std::vector<uint8_t>* expected) {

    expected->assign(
        ks256gpu::kTotalSlots,
        0);

    for (uint32_t word = 0;
         word < ks256gpu::kNumWords;
         ++word) {
        const uint32_t base =
            word * ks256gpu::kStride;

        const uint32_t shift =
            EffectiveShift(encryptedShift, word);

        for (uint32_t bit = shift;
             bit < ks256gpu::kWordBits;
             ++bit) {
            (*expected)[base + bit] =
                input[base + bit - shift];
        }
    }
}

void PrintPhase(const std::string& name,
                double milliseconds) {
    std::cout
        << "  "
        << std::left
        << std::setw(40)
        << name
        << std::right
        << std::setw(14)
        << std::fixed
        << std::setprecision(3)
        << milliseconds
        << " ms"
        << "  ("
        << std::setprecision(6)
        << milliseconds / 1000.0
        << " s)\n";
}

void PrintShiftProfile(
    const ks256gpu::TimingReport& profile) {

    std::cout
        << "\n================ CPU 瀵嗘枃绉讳綅閲忔《褰㈠乏绉诲唴閮ㄨ鏃?================\n";

    PrintPhase(
        "绉讳綅閲忎綆 8 浣嶆槑鏂囨帺鐮?,
        profile.lowByteMaskMs);

    double patternTotal = 0.0;

    for (uint32_t layer = 0;
         layer < ks256gpu::kPatternLayers;
         ++layer) {
        const auto& item =
            profile.patternLayers[layer];

        patternTotal += item.wallMs;

        std::cout
            << "\n  Control Pattern Layer "
            << layer
            << "  (step="
            << (8U << layer)
            << ")\n";

        PrintPhase(
            "    Rotate(low-byte pattern)",
            item.rotateMs);

        PrintPhase(
            "    Pattern + Rotate(pattern)",
            item.addMs);

        PrintPhase(
            "    鏈眰鎬诲閽熸椂闂?,
            item.wallMs);
    }

    std::cout << "\n";

    PrintPhase(
        "5 灞備綆 8 浣嶆ā寮忓鍒跺悎璁?,
        patternTotal);

    double barrelTotal = 0.0;

    for (uint32_t layer = 0;
         layer < ks256gpu::kBarrelLayers;
         ++layer) {
        const auto& item =
            profile.barrelLayers[layer];

        barrelTotal += item.wallMs;

        std::cout
            << "\n  Barrel Layer "
            << layer
            << "  (conditional shift="
            << (1U << layer)
            << ")\n";

        PrintPhase(
            "    鎻愬彇褰撳墠鎺у埗浣?lane",
            item.controlLaneMaskMs);

        PrintPhase(
            "    鎺у埗浣嶅榻愬埌 8-slot block",
            item.controlAnchorMs);

        for (uint32_t broadcastLayer = 0;
             broadcastLayer < ks256gpu::kBroadcastLayers;
             ++broadcastLayer) {
            const std::string prefix =
                "    鎺у埗浣嶅箍鎾?step=" +
                std::to_string(1U << broadcastLayer);

            PrintPhase(
                prefix + " Rotate",
                item.controlBroadcastRotateMs[broadcastLayer]);

            PrintPhase(
                prefix + " Add",
                item.controlBroadcastAddMs[broadcastLayer]);
        }

        PrintPhase(
            "    鎺у埗瀵嗘枃鏋勯€犳€绘椂闂?,
            item.controlBuildWallMs);

        PrintPhase(
            "    鏁版嵁鍥哄畾璺濈鏃嬭浆",
            item.dataRotateMs);

        PrintPhase(
            "    宸︾Щ杈圭晫鏄庢枃鎺╃爜",
            item.shiftMaskMs);

        PrintPhase(
            "    candidate-current",
            item.candidateDiffMs);

        PrintPhase(
            "    control*(candidate-current)",
            item.selectMultMs);

        PrintPhase(
            "    current+selected difference",
            item.selectAddMs);

        PrintPhase(
            "    鏈眰鎬诲閽熸椂闂?,
            item.wallMs);
    }

    std::cout << "\n";

    PrintPhase(
        "8 灞傛《褰㈠乏绉诲悎璁?,
        barrelTotal);

    PrintPhase(
        "EvalShiftLeftEncrypted 鎬绘椂闂?,
        profile.totalMs);
}

}  // namespace

int main() {
    try {
        const auto programBegin =
            Clock::now();

        PhaseTimes phase;

        constexpr uint32_t ringDim =
            1U << 17;  // 131072

        constexpr uint32_t numSlots =
            ks256gpu::kTotalSlots;  // 65536

        constexpr uint32_t multiplicativeDepth = 25;

        constexpr uint64_t randomSeed =
            0x53484C454E433235ULL;

        const SecretKeyDist secretKeyDist =
            UNIFORM_TERNARY;

        const std::vector<uint32_t> levelBudget =
            {3, 3};

        const std::vector<uint32_t> bsgsDim =
            {0, 0};

        std::cout
            << "============================================================\n";

        std::cout
            << "OpenFHE CPU: "
            << ks256gpu::kNumWords
            << " x uint256 瀵嗘枃绉讳綅閲忔《褰㈠乏绉?+ 涓€娆?Bootstrap\n";

        std::cout
            << "============================================================\n";

        // --------------------------------------------------------------
        // 1. Create the same CKKS/OpenFHE environment used by the
        //    encrypted comparison operators.
        // --------------------------------------------------------------
        CCParams<CryptoContextCKKSRNS> parameters;

        parameters.SetSecretKeyDist(
            secretKeyDist);

        parameters.SetSecurityLevel(
            HEStd_128_classic);

        parameters.SetRingDim(
            ringDim);

        parameters.SetBatchSize(
            numSlots);

        parameters.SetMultiplicativeDepth(
            multiplicativeDepth);

        parameters.SetScalingTechnique(
            FLEXIBLEAUTO);

        parameters.SetScalingModSize(
            59);

        parameters.SetFirstModSize(
            60);

        parameters.SetKeySwitchTechnique(
            HYBRID);

        parameters.SetNumLargeDigits(
            3);

        const std::vector<int> devices{0, 1};
auto begin = Clock::now();

        CryptoContext<DCRTPoly> cc =
            GenCryptoContext(parameters);

        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);
        cc->Enable(FHE);

        phase.contextMs =
            Milliseconds(
                begin,
                Clock::now());

        std::cout << "鍙傛暟:\n";

        std::cout
            << "  Ring dimension              : "
            << cc->GetRingDimension()
            << "\n";

        std::cout
            << "  CKKS slots                  : "
            << numSlots
            << "\n";

        std::cout
            << "  Parallel 256-bit shifts    : "
            << ks256gpu::kNumWords
            << "\n";

        std::cout
            << "  Slots per shift            : "
            << ks256gpu::kStride
            << " (256 active + 256 guard)\n";

        std::cout
            << "  Data operand width         : uint256\n";

        std::cout
            << "  Encrypted shift width      : uint256\n";

        std::cout
            << "  Effective shift semantics  : low 8 bits (S mod 256)\n";

        std::cout
            << "  Barrel-shifter layers      : "
            << ks256gpu::kBarrelLayers
            << "\n";

        std::cout
            << "  Ciphertext multiplications : 8\n";

        std::cout
            << "  Configured mult depth      : "
            << multiplicativeDepth
            << "\n";

        std::cout
            << "  Bootstrap level budget     : {"
            << levelBudget[0]
            << ", "
            << levelBudget[1]
            << "}\n";

        std::cout
            << "  OpenFHE visible devices   : 0, 1\n";

        // --------------------------------------------------------------
        // 2. Evaluation keys before LoadContext.
        //
        // Negative powers of two are used by data shifts, low-byte pattern
        // replication, and 8-slot control broadcast. Positive 1...7 align
        // each control lane to the start of its 8-slot block.
        // --------------------------------------------------------------
        begin = Clock::now();

        auto keys = cc->KeyGen();

        phase.keyGenMs =
            Milliseconds(
                begin,
                Clock::now());

        begin = Clock::now();

        cc->EvalMultKeyGen(
            keys.secretKey);

        phase.multKeyGenMs =
            Milliseconds(
                begin,
                Clock::now());

        std::vector<int32_t> rotationIndices;
        rotationIndices.reserve(15);

        for (uint32_t layer = 0;
             layer < ks256gpu::kBarrelLayers;
             ++layer) {
            rotationIndices.push_back(
                -static_cast<int32_t>(1U << layer));
        }

        for (uint32_t lane = 1;
             lane < ks256gpu::kShiftControlBits;
             ++lane) {
            rotationIndices.push_back(
                static_cast<int32_t>(lane));
        }

        begin = Clock::now();

        cc->EvalRotateKeyGen(
            keys.secretKey,
            rotationIndices);

        phase.rotateKeyGenMs =
            Milliseconds(
                begin,
                Clock::now());

        begin = Clock::now();

        cc->EvalBootstrapSetup(
            levelBudget,
            bsgsDim,
            numSlots,
            0);

        phase.bootstrapSetupMs =
            Milliseconds(
                begin,
                Clock::now());

        begin = Clock::now();

        cc->EvalBootstrapKeyGen(
            keys.secretKey,
            numSlots);

        phase.bootstrapKeyGenMs =
            Milliseconds(
                begin,
                Clock::now());

        begin = Clock::now();
phase.loadContextMs =
            Milliseconds(
                begin,
                Clock::now());

        begin = Clock::now();

        auto masks =
            ks256gpu::BuildAndLoadMasks(cc);

        phase.maskBuildLoadMs =
            Milliseconds(
                begin,
                Clock::now());

        // --------------------------------------------------------------
        // 3. Build 128 independent uint256 data/shift pairs.
        // --------------------------------------------------------------
        begin = Clock::now();

        std::vector<uint8_t> bitsData(
            ks256gpu::kTotalSlots,
            0);

        std::vector<uint8_t> bitsShift(
            ks256gpu::kTotalSlots,
            0);

        std::vector<double> slotsData(
            ks256gpu::kTotalSlots,
            0.0);

        std::vector<double> slotsShift(
            ks256gpu::kTotalSlots,
            0.0);

        std::mt19937_64 rng(
            randomSeed);

        for (uint32_t word = 0;
             word < ks256gpu::kNumWords;
             ++word) {
            const uint32_t base =
                word * ks256gpu::kStride;

            for (uint32_t bit = 0;
                 bit < ks256gpu::kWordBits;
                 ++bit) {
                const uint8_t dataBit =
                    static_cast<uint8_t>(
                        rng() & 1ULL);

                const uint8_t shiftBit =
                    static_cast<uint8_t>(
                        rng() & 1ULL);

                bitsData[base + bit] = dataBit;
                bitsShift[base + bit] = shiftBit;

                slotsData[base + bit] =
                    static_cast<double>(dataBit);

                slotsShift[base + bit] =
                    static_cast<double>(shiftBit);
            }
        }

        // Deterministic shift-control boundary cases.
        SetShiftWord(&bitsShift, &slotsShift, 0, 0);
        SetShiftWord(&bitsShift, &slotsShift, 1, 1);
        SetShiftWord(&bitsShift, &slotsShift, 2, 2);
        SetShiftWord(&bitsShift, &slotsShift, 3, 7);
        SetShiftWord(&bitsShift, &slotsShift, 4, 8);
        SetShiftWord(&bitsShift, &slotsShift, 5, 31);
        SetShiftWord(&bitsShift, &slotsShift, 6, 63);
        SetShiftWord(&bitsShift, &slotsShift, 7, 127);
        SetShiftWord(&bitsShift, &slotsShift, 8, 128);
        SetShiftWord(&bitsShift, &slotsShift, 9, 255);

        // Full uint256 shift amount 256: low eight bits are zero, so the
        // modulo-256 operator performs shift 0.
        SetShiftWord(
            &bitsShift,
            &slotsShift,
            10,
            0,
            {8});

        // Full uint256 shift amount 257: low eight bits are one.
        SetShiftWord(
            &bitsShift,
            &slotsShift,
            11,
            1,
            {8});

        // All 256 shift bits equal one: effective shift is 255.
        {
            const uint32_t base =
                12U * ks256gpu::kStride;

            for (uint32_t bit = 0;
                 bit < ks256gpu::kWordBits;
                 ++bit) {
                bitsShift[base + bit] = 1;
                slotsShift[base + bit] = 1.0;
            }
        }

        // Deterministic data patterns for the first few cases.
        for (uint32_t word = 0; word < 13; ++word) {
            const uint32_t base =
                word * ks256gpu::kStride;

            for (uint32_t bit = 0;
                 bit < ks256gpu::kWordBits;
                 ++bit) {
                const uint8_t value =
                    static_cast<uint8_t>(
                        ((bit + word) % 5U) == 0U);

                bitsData[base + bit] = value;
                slotsData[base + bit] =
                    static_cast<double>(value);
            }
        }

        std::vector<uint8_t> expected;

        BuildExpectedShiftLeft(
            bitsData,
            bitsShift,
            &expected);

        phase.inputBuildMs =
            Milliseconds(
                begin,
                Clock::now());

        begin = Clock::now();

        Plaintext pData =
            cc->MakeCKKSPackedPlaintext(
                slotsData,
                1,
                0,
                nullptr,
                numSlots);

        Plaintext pShift =
            cc->MakeCKKSPackedPlaintext(
                slotsShift,
                1,
                0,
                nullptr,
                numSlots);

        pData->SetLength(numSlots);
        pShift->SetLength(numSlots);

        phase.encodingMs =
            Milliseconds(
                begin,
                Clock::now());

        begin = Clock::now();

        auto cData =
            cc->Encrypt(
                keys.publicKey,
                pData);

        auto cShift =
            cc->Encrypt(
                keys.publicKey,
                pShift);
phase.encryptionMs =
            Milliseconds(
                begin,
                Clock::now());

        std::cout << "\n杈撳叆鏍蜂緥:\n";

        std::cout
            << "  word 0 data            = "
            << BitsToHex(bitsData, 0)
            << "\n";

        std::cout
            << "  word 0 encrypted shift = "
            << BitsToHex(bitsShift, 0)
            << "\n";

        std::cout
            << "  word 0 effective shift = "
            << EffectiveShift(bitsShift, 0)
            << "\n";

        std::cout
            << "  word 0 expected        = "
            << BitsToHex(expected, 0)
            << "\n";

        std::cout
            << "  word 10 data            = "
            << BitsToHex(bitsData, 10)
            << "\n";

        std::cout
            << "  word 10 encrypted shift = "
            << BitsToHex(bitsShift, 10)
            << "\n";

        std::cout
            << "  word 10 effective shift = "
            << EffectiveShift(bitsShift, 10)
            << " (full shift amount contains bit 8)\n";

        std::cout
            << "  word 10 expected        = "
            << BitsToHex(expected, 10)
            << "\n";

        // --------------------------------------------------------------
        // 4. CPU encrypted-control barrel shifter.
        // --------------------------------------------------------------
        ks256gpu::TimingReport shiftTiming;
begin = Clock::now();

        auto cShifted =
            ks256gpu::EvalShiftLeftEncrypted(
                cc,
                cData,
                cShift,
                masks,
                &shiftTiming);
phase.shiftMs =
            Milliseconds(
                begin,
                Clock::now());

        std::cout
            << "\n瀵嗘枃绉讳綅閲忓乏绉诲悗瀵嗘枃鐘舵€?\n";

        std::cout
            << "  GetLevel()         : "
            << cShifted->GetLevel()
            << "\n";

        std::cout
            << "  GetNoiseScaleDeg() : "
            << cShifted->GetNoiseScaleDeg()
            << "\n";

        std::cout
            << "  浼扮畻鍓╀綑灞傛暟      : "
            << LevelsRemaining(
                   multiplicativeDepth,
                   cShifted)
            << "\n";

        // --------------------------------------------------------------
        // 5. Exactly one CPU bootstrap after the complete barrel shifter.
        // --------------------------------------------------------------
begin = Clock::now();

        auto cBoot =
            cc->EvalBootstrap(
                cShifted);
phase.bootstrapMs =
            Milliseconds(
                begin,
                Clock::now());

        std::cout
            << "\nBootstrap 鍚庡瘑鏂囩姸鎬?\n";

        std::cout
            << "  GetLevel()         : "
            << cBoot->GetLevel()
            << "\n";

        std::cout
            << "  GetNoiseScaleDeg() : "
            << cBoot->GetNoiseScaleDeg()
            << "\n";

        std::cout
            << "  浼扮畻鍓╀綑灞傛暟      : "
            << LevelsRemaining(
                   multiplicativeDepth,
                   cBoot)
            << "\n";

        // --------------------------------------------------------------
        // 6. Copy to CPU and decrypt after bootstrap.
        // --------------------------------------------------------------
        begin = Clock::now();

        Plaintext decrypted;

        const auto decryptResult =
            cc->Decrypt(
                keys.secretKey,
                cBoot,
                &decrypted);

        if (!decryptResult.isValid) {
            throw std::runtime_error(
                "Decryption returned an invalid result");
        }

        decrypted->SetLength(numSlots);

        const std::vector<double> values =
            decrypted->GetRealPackedValue();

        phase.decryptionMs =
            Milliseconds(
                begin,
                Clock::now());

        // --------------------------------------------------------------
        // 7. Recover all active Boolean slots and scan CKKS error.
        // --------------------------------------------------------------
        begin = Clock::now();

        double maxAbsError = 0.0;
        double sumAbsError = 0.0;
        double sumSquaredError = 0.0;
        double maxGuardAbs = 0.0;

        uint64_t bitFailures = 0;
        uint32_t worstWord = 0;
        uint32_t worstBit = 0;
        double worstValue = 0.0;
        double worstExpected = 0.0;

        const uint64_t activeCount =
            static_cast<uint64_t>(
                ks256gpu::kNumWords) *
            ks256gpu::kWordBits;

        for (uint32_t word = 0;
             word < ks256gpu::kNumWords;
             ++word) {
            const uint32_t base =
                word * ks256gpu::kStride;

            for (uint32_t bit = 0;
                 bit < ks256gpu::kWordBits;
                 ++bit) {
                const uint32_t index =
                    base + bit;

                const double target =
                    static_cast<double>(
                        expected[index]);

                const double error =
                    std::abs(
                        values[index] - target);

                sumAbsError += error;
                sumSquaredError += error * error;

                if (error > maxAbsError) {
                    maxAbsError = error;
                    worstWord = word;
                    worstBit = bit;
                    worstValue = values[index];
                    worstExpected = target;
                }

                const uint8_t recovered =
                    values[index] >= 0.5
                        ? 1U
                        : 0U;

                if (recovered != expected[index]) {
                    ++bitFailures;
                }
            }

            for (uint32_t guard = 0;
                 guard < ks256gpu::kGuardBits;
                 ++guard) {
                const uint32_t index =
                    base +
                    ks256gpu::kWordBits +
                    guard;

                maxGuardAbs =
                    std::max(
                        maxGuardAbs,
                        std::abs(values[index]));
            }
        }

        const double meanAbsError =
            sumAbsError /
            static_cast<double>(
                activeCount);

        const double rmse =
            std::sqrt(
                sumSquaredError /
                static_cast<double>(
                    activeCount));

        phase.verificationMs =
            Milliseconds(
                begin,
                Clock::now());

        phase.totalMs =
            Milliseconds(
                programBegin,
                Clock::now());

        std::cout
            << "\n================ 瑙ｅ瘑姝ｇ‘鎬т笌缁忛獙绮惧害 ================\n";

        std::cout
            << std::scientific
            << std::setprecision(10);

        std::cout
            << "  鏈夋晥 bit 妲芥暟                 : "
            << activeCount
            << "\n";

        std::cout
            << "  鏈€澶х粷瀵硅宸?                : "
            << maxAbsError
            << "\n";

        std::cout
            << "  骞冲潎 L1 璇樊                 : "
            << meanAbsError
            << "\n";

        std::cout
            << "  RMSE                         : "
            << rmse
            << "\n";

        std::cout
            << "  Guard 妲芥渶澶х粷瀵瑰€?          : "
            << maxGuardAbs
            << "\n";

        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "  GetLogPrecision              : "
            << decrypted->GetLogPrecision()
            << " bits\n";

        std::cout
            << "  鏈€宸簿搴?-log2(max error)    : "
            << PrecisionBits(maxAbsError)
            << " bits\n";

        std::cout
            << "  骞冲潎绮惧害 -log2(mean L1)      : "
            << PrecisionBits(meanAbsError)
            << " bits\n";

        std::cout
            << "  RMSE 绮惧害 -log2(RMSE)        : "
            << PrecisionBits(rmse)
            << " bits\n";

        std::cout
            << "  0.5 闃堝€兼仮澶嶅け璐ユ暟           : "
            << bitFailures
            << "\n";

        std::cout
            << "  鏈€宸綅缃?                    : word="
            << worstWord
            << ", bit="
            << worstBit
            << "\n";

        std::cout
            << std::scientific
            << std::setprecision(12);

        std::cout
            << "  鏈€宸綅缃湡鏈?瑙ｅ瘑鍊?         : "
            << worstExpected
            << " / "
            << worstValue
            << "\n";

        std::cout << std::fixed;

        std::cout
            << "  鏁翠綋缁撴灉                     : "
            << (bitFailures == 0
                    ? "PASS"
                    : "FAIL")
            << "\n";

        std::cout
            << "璇存槑锛氫笂杩?bit 鏁版槸 CKKS 瑙ｅ瘑璇樊瀵瑰簲鐨勭粡楠岀簿搴︼紝"
               "涓嶆槸 BFV/BGV 鐨勭簿纭櫔澹伴绠椼€俓n";

        PrintShiftProfile(
            shiftTiming);

        std::cout
            << "\n================ 绔埌绔椂闂村垎鏋?================\n";

        PrintPhase(
            "CryptoContext 鐢熸垚涓?Enable",
            phase.contextMs);

        PrintPhase(
            "KeyGen",
            phase.keyGenMs);

        PrintPhase(
            "EvalMultKeyGen",
            phase.multKeyGenMs);

        PrintPhase(
            "妗跺舰绉讳綅鏃嬭浆瀵嗛挜鐢熸垚",
            phase.rotateKeyGenMs);

        PrintPhase(
            "EvalBootstrapSetup",
            phase.bootstrapSetupMs);

        PrintPhase(
            "EvalBootstrapKeyGen",
            phase.bootstrapKeyGenMs);

        PrintPhase(
            "LoadContext 涓庤瘎浼板瘑閽ュ姞杞藉埌 CPU",
            phase.loadContextMs);

        PrintPhase(
            "妗跺舰绉讳綅鎺╃爜缂栫爜骞跺姞杞藉埌 CPU",
            phase.maskBuildLoadMs);

        PrintPhase(
            "杈撳叆鐢熸垚涓?CPU 鍙傝€冪粨鏋?,
            phase.inputBuildMs);

        PrintPhase(
            "CKKS 缂栫爜",
            phase.encodingMs);

        PrintPhase(
            "杈撳叆鍔犲瘑骞跺姞杞藉埌 CPU",
            phase.encryptionMs);

        PrintPhase(
            "CPU 瀵嗘枃绉讳綅閲忔《褰㈠乏绉?,
            phase.shiftMs);

        PrintPhase(
            "CPU 涓€娆?EvalBootstrap",
            phase.bootstrapMs);

        PrintPhase(
            "澶嶅埗鍥?CPU 骞惰В瀵?,
            phase.decryptionMs);

        PrintPhase(
            "姝ｇ‘鎬т笌璇樊鎵弿",
            phase.verificationMs);

        PrintPhase(
            "绋嬪簭鎬诲閽熸椂闂?,
            phase.totalMs);

        const double shiftSeconds =
            phase.shiftMs / 1000.0;

        const double bootSeconds =
            phase.bootstrapMs / 1000.0;

        const double onlineSeconds =
            shiftSeconds + bootSeconds;

        std::cout
            << "\n================ 鍚炲悙涓庢憡閿€ ================\n";

        std::cout
            << std::fixed
            << std::setprecision(6);

        std::cout
            << "  宸︾Щ鎽婇攢鏃堕棿/256-bit word     : "
            << shiftSeconds /
                   ks256gpu::kNumWords
            << " s\n";

        std::cout
            << "  Boot 鎽婇攢鏃堕棿/word            : "
            << bootSeconds /
                   ks256gpu::kNumWords
            << " s\n";

        std::cout
            << "  鍦ㄧ嚎鎬绘憡閿€鏃堕棿/word           : "
            << onlineSeconds /
                   ks256gpu::kNumWords
            << " s\n";

        if (shiftSeconds > 0.0) {
            std::cout
                << "  浠呭瘑鏂囩Щ浣嶉噺宸︾Щ鍚炲悙閲?       : "
                << ks256gpu::kNumWords /
                       shiftSeconds
                << " word/s\n";
        }

        if (onlineSeconds > 0.0) {
            std::cout
                << "  宸︾Щ+Boot 鍦ㄧ嚎鍚炲悙閲?         : "
                << ks256gpu::kNumWords /
                       onlineSeconds
                << " word/s\n";
        }

        return bitFailures == 0 ? 0 : 2;
    }
    catch (const std::exception& error) {
        std::cerr
            << "\n[FATAL] "
            << error.what()
            << "\n";

        return 1;
    }
}




