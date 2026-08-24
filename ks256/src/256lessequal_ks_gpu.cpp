#include <openfhe.h>

#include "lessequal_ks_gpu.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
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
    double comparisonMs      = 0.0;
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

        output[output.size() - 1 - nibble] =
            kHex[value];
    }

    return "0x" + output;
}

void BuildExpectedLessEqual(
    const std::vector<uint8_t>& bitsA,
    const std::vector<uint8_t>& bitsB,
    std::vector<uint8_t>* expected) {

    expected->assign(
        ks256gpu::kNumWords,
        0);

    for (uint32_t word = 0;
         word < ks256gpu::kNumWords;
         ++word) {
        const uint32_t base =
            word * ks256gpu::kStride;

        uint8_t lessEqual = 1;

        for (int32_t bit =
                 static_cast<int32_t>(ks256gpu::kWordBits) - 1;
             bit >= 0;
             --bit) {
            const uint8_t a =
                bitsA[base + static_cast<uint32_t>(bit)];

            const uint8_t b =
                bitsB[base + static_cast<uint32_t>(bit)];

            if (a != b) {
                lessEqual = a < b ? 1U : 0U;
                break;
            }
        }

        (*expected)[word] = lessEqual;
    }
}

void PrintPhase(const std::string& name,
                double milliseconds) {
    std::cout << "  "
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

void PrintComparisonProfile(
    const ks256gpu::TimingReport& profile) {

    std::cout
        << "\n================ CPU 鏈€楂樹笉鍚屼綅 P-only 灏忎簬绛変簬姣旇緝鍐呴儴璁℃椂 ================\n";

    PrintPhase(
        "棰勫鐞?XNOR",
        profile.preprocessMs);

    double prefixTotal = 0.0;

    for (uint32_t layer = 0;
         layer < ks256gpu::kNumLayers;
         ++layer) {
        const auto& item =
            profile.prefixLayers[layer];

        prefixTotal += item.wallMs;

        std::cout << "\n  Suffix Layer "
                  << layer
                  << "  (step="
                  << (1U << layer)
                  << ")\n";

        PrintPhase(
            "    Rotate(P)",
            item.rotatePMs);

        PrintPhase(
            "    P * Rotate(P)",
            item.prefixMultMs);

        PrintPhase(
            "    鏈眰鎬诲閽熸椂闂?,
            item.wallMs);
    }

    std::cout << "\n";
    PrintPhase(
        "8 灞?P-only 鍚庣紑 AND 鍚堣",
        prefixTotal);

    PrintPhase(
        "E[i+1] 鏃嬭浆",
        profile.higherEqualRotateMs);

    PrintPhase(
        "鏈€楂樹笉鍚屼綅 M=E[i+1]-E[i]",
        profile.highestDifferentMs);

    PrintPhase(
        "涓ユ牸澶т簬鏍稿績閫夋嫨 M*A",
        profile.selectAMs);

    double reductionTotal = 0.0;

    for (uint32_t layer = 0;
         layer < ks256gpu::kNumLayers;
         ++layer) {
        const auto& item =
            profile.reductionLayers[layer];

        reductionTotal += item.wallMs;

        std::cout << "\n  Reduction Layer "
                  << layer
                  << "  (step="
                  << (1U << layer)
                  << ")\n";

        PrintPhase(
            "    Rotate(T)",
            item.rotateMs);

        PrintPhase(
            "    T + Rotate(T)",
            item.addMs);

        PrintPhase(
            "    鏈眰鎬诲閽熸椂闂?,
            item.wallMs);
    }

    std::cout << "\n";
    PrintPhase(
        "8 灞傛棆杞姞娉曞綊绾﹀悎璁?,
        reductionTotal);

    PrintPhase(
        "涓ユ牸澶т簬缁撴灉鏄庢枃鎺╃爜涔樻硶",
        profile.resultMaskMs);

    PrintPhase(
        "缁撴灉鍙栧弽 LE=resultMask-GT",
        profile.negateMs);

    PrintPhase(
        "EvalLessEqualPOnly 鎬绘椂闂?,
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
            0x4752454154455232ULL;

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
            << " x 256-bit 鏈€楂樹笉鍚屼綅 P-only 鏃犵鍙峰皬浜庣瓑浜庢瘮杈?+ 涓€娆?Bootstrap\n";

        std::cout
            << "============================================================\n";

        // --------------------------------------------------------------
        // 1. Create a two-CPU OpenFHE CKKS context.
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
            << "  Parallel 256-bit comparisons: "
            << ks256gpu::kNumWords
            << "\n";

        std::cout
            << "  Slots per comparison       : "
            << ks256gpu::kStride
            << " (256 active + 256 guard)\n";

        std::cout
            << "  Result slot per word       : bit 0\n";

        std::cout
            << "  Comparison semantics       : unsigned A <= B\n";

        std::cout
            << "  Ciphertext multiplications : 10\n";

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
        // 2. Generate evaluation keys before LoadContext.
        //
        // +1,+2,...,+128 are used by both:
        //   - P-only suffix-AND;
        //   - rotate-add reduction to bit 0.
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

        rotationIndices.reserve(
            ks256gpu::kNumLayers);

        for (uint32_t layer = 0;
             layer < ks256gpu::kNumLayers;
             ++layer) {
            rotationIndices.push_back(
                static_cast<int32_t>(
                    1U << layer));
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
        // 3. Generate 128 independent 256-bit input pairs.
        // --------------------------------------------------------------
        begin = Clock::now();

        std::vector<uint8_t> bitsA(
            ks256gpu::kTotalSlots,
            0);

        std::vector<uint8_t> bitsB(
            ks256gpu::kTotalSlots,
            0);

        std::vector<double> slotsA(
            ks256gpu::kTotalSlots,
            0.0);

        std::vector<double> slotsB(
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
                const uint8_t a =
                    static_cast<uint8_t>(
                        rng() & 1ULL);

                const uint8_t b =
                    static_cast<uint8_t>(
                        rng() & 1ULL);

                bitsA[base + bit] = a;
                bitsB[base + bit] = b;

                slotsA[base + bit] =
                    static_cast<double>(a);

                slotsB[base + bit] =
                    static_cast<double>(b);
            }

            // Periodically create equal operands.
            if ((word % 8U) == 0U) {
                for (uint32_t bit = 0;
                     bit < ks256gpu::kWordBits;
                     ++bit) {
                    bitsB[base + bit] =
                        bitsA[base + bit];

                    slotsB[base + bit] =
                        slotsA[base + bit];
                }
            }
        }

        // Test 0: 0 <= 0 is true.
        {
            const uint32_t base = 0;

            for (uint32_t bit = 0;
                 bit < ks256gpu::kWordBits;
                 ++bit) {
                bitsA[base + bit] = 0;
                bitsB[base + bit] = 0;
                slotsA[base + bit] = 0.0;
                slotsB[base + bit] = 0.0;
            }
        }

        // Test 1: operands differ only at bit 0, A=1 and B=0; A<=B is false.
        {
            const uint32_t base =
                ks256gpu::kStride;

            for (uint32_t bit = 0;
                 bit < ks256gpu::kWordBits;
                 ++bit) {
                const uint8_t value =
                    static_cast<uint8_t>(
                        (bit & 1U) != 0U);

                bitsA[base + bit] = value;
                bitsB[base + bit] = value;
                slotsA[base + bit] =
                    static_cast<double>(value);
                slotsB[base + bit] =
                    static_cast<double>(value);
            }

            bitsA[base] = 1;
            bitsB[base] = 0;
            slotsA[base] = 1.0;
            slotsB[base] = 0.0;
        }

        // Test 2: 2^256-1 <= 0 is false.
        {
            const uint32_t base =
                2U * ks256gpu::kStride;

            for (uint32_t bit = 0;
                 bit < ks256gpu::kWordBits;
                 ++bit) {
                bitsA[base + bit] = 1;
                bitsB[base + bit] = 0;
                slotsA[base + bit] = 1.0;
                slotsB[base + bit] = 0.0;
            }
        }

        // Test 3: operands differ only at bit 255, A=1 and B=0; A<=B is false.
        {
            const uint32_t base =
                3U * ks256gpu::kStride;

            for (uint32_t bit = 0;
                 bit < ks256gpu::kWordBits;
                 ++bit) {
                const uint8_t value =
                    static_cast<uint8_t>(
                        (bit % 3U) == 0U);

                bitsA[base + bit] = value;
                bitsB[base + bit] = value;
                slotsA[base + bit] =
                    static_cast<double>(value);
                slotsB[base + bit] =
                    static_cast<double>(value);
            }

            const uint32_t msb =
                base + ks256gpu::kWordBits - 1U;

            bitsA[msb] = 1;
            bitsB[msb] = 0;
            slotsA[msb] = 1.0;
            slotsB[msb] = 0.0;
        }

        // Test 4: operands differ only at bit 255, A=0 and B=1; A<=B is true.
        {
            const uint32_t base =
                4U * ks256gpu::kStride;

            for (uint32_t bit = 0;
                 bit < ks256gpu::kWordBits;
                 ++bit) {
                const uint8_t value =
                    static_cast<uint8_t>(
                        (bit % 5U) == 0U);

                bitsA[base + bit] = value;
                bitsB[base + bit] = value;
                slotsA[base + bit] =
                    static_cast<double>(value);
                slotsB[base + bit] =
                    static_cast<double>(value);
            }

            const uint32_t msb =
                base + ks256gpu::kWordBits - 1U;

            bitsA[msb] = 0;
            bitsB[msb] = 1;
            slotsA[msb] = 0.0;
            slotsB[msb] = 1.0;
        }

        // Test 5: operands differ only at bit 0, A=0 and B=1; A<=B is true.
        {
            const uint32_t base =
                5U * ks256gpu::kStride;

            for (uint32_t bit = 0;
                 bit < ks256gpu::kWordBits;
                 ++bit) {
                const uint8_t value =
                    static_cast<uint8_t>(
                        (bit % 7U) == 0U);

                bitsA[base + bit] = value;
                bitsB[base + bit] = value;
                slotsA[base + bit] =
                    static_cast<double>(value);
                slotsB[base + bit] =
                    static_cast<double>(value);
            }

            bitsA[base] = 0;
            bitsB[base] = 1;
            slotsA[base] = 0.0;
            slotsB[base] = 1.0;
        }

        std::vector<uint8_t> expected;

        BuildExpectedLessEqual(
            bitsA,
            bitsB,
            &expected);

        phase.inputBuildMs =
            Milliseconds(
                begin,
                Clock::now());

        begin = Clock::now();

        Plaintext pA =
            cc->MakeCKKSPackedPlaintext(
                slotsA,
                1,
                0,
                nullptr,
                numSlots);

        Plaintext pB =
            cc->MakeCKKSPackedPlaintext(
                slotsB,
                1,
                0,
                nullptr,
                numSlots);

        pA->SetLength(numSlots);
        pB->SetLength(numSlots);

        phase.encodingMs =
            Milliseconds(
                begin,
                Clock::now());

        begin = Clock::now();

        auto cA =
            cc->Encrypt(
                keys.publicKey,
                pA);

        auto cB =
            cc->Encrypt(
                keys.publicKey,
                pB);
phase.encryptionMs =
            Milliseconds(
                begin,
                Clock::now());

        std::cout << "\n杈撳叆鏍蜂緥:\n";

        std::cout
            << "  word 0 A  = "
            << BitsToHex(bitsA, 0)
            << "\n";

        std::cout
            << "  word 0 B  = "
            << BitsToHex(bitsB, 0)
            << "\n";

        std::cout
            << "  word 0 LE = "
            << static_cast<uint32_t>(
                   expected[0])
            << "\n";

        std::cout
            << "  word 1 A  = "
            << BitsToHex(bitsA, 1)
            << "\n";

        std::cout
            << "  word 1 B  = "
            << BitsToHex(bitsB, 1)
            << "\n";

        std::cout
            << "  word 1 LE = "
            << static_cast<uint32_t>(
                   expected[1])
            << "\n";

        // --------------------------------------------------------------
        // 4. CPU highest-different-bit P-only less-than-or-equal circuit.
        // --------------------------------------------------------------
        ks256gpu::TimingReport comparisonTiming;
begin = Clock::now();

        auto cLessEqual =
            ks256gpu::EvalLessEqualPOnly(
                cc,
                cA,
                cB,
                masks,
                &comparisonTiming);
phase.comparisonMs =
            Milliseconds(
                begin,
                Clock::now());

        std::cout
            << "\n灏忎簬绛変簬鍒ゆ柇鍚庡瘑鏂囩姸鎬?\n";

        std::cout
            << "  GetLevel()         : "
            << cLessEqual->GetLevel()
            << "\n";

        std::cout
            << "  GetNoiseScaleDeg() : "
            << cLessEqual->GetNoiseScaleDeg()
            << "\n";

        std::cout
            << "  浼扮畻鍓╀綑灞傛暟      : "
            << LevelsRemaining(
                   multiplicativeDepth,
                   cLessEqual)
            << "\n";

        // --------------------------------------------------------------
        // 5. Exactly one CPU bootstrap after the less-than-or-equal circuit.
        // --------------------------------------------------------------
begin = Clock::now();

        auto cBoot =
            cc->EvalBootstrap(
                cLessEqual);
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
        // 7. Recover the 128 comparison bits and scan CKKS error.
        // --------------------------------------------------------------
        begin = Clock::now();

        double maxAbsError = 0.0;
        double sumAbsError = 0.0;
        double sumSquaredError = 0.0;
        double maxInactiveAbs = 0.0;

        uint64_t bitFailures = 0;
        uint32_t worstWord = 0;
        double worstValue = 0.0;
        double worstExpected = 0.0;

        for (uint32_t word = 0;
             word < ks256gpu::kNumWords;
             ++word) {
            const uint32_t base =
                word * ks256gpu::kStride;

            const double target =
                static_cast<double>(
                    expected[word]);

            const double error =
                std::abs(
                    values[base] - target);

            sumAbsError += error;
            sumSquaredError += error * error;

            if (error > maxAbsError) {
                maxAbsError = error;
                worstWord = word;
                worstValue = values[base];
                worstExpected = target;
            }

            const uint8_t recovered =
                values[base] >= 0.5
                    ? 1U
                    : 0U;

            if (recovered != expected[word]) {
                ++bitFailures;
            }

            for (uint32_t offset = 1;
                 offset < ks256gpu::kStride;
                 ++offset) {
                maxInactiveAbs =
                    std::max(
                        maxInactiveAbs,
                        std::abs(
                            values[base + offset]));
            }
        }

        const double meanAbsError =
            sumAbsError /
            static_cast<double>(
                ks256gpu::kNumWords);

        const double rmse =
            std::sqrt(
                sumSquaredError /
                static_cast<double>(
                    ks256gpu::kNumWords));

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
            << "  鏈夋晥 comparison 妲芥暟          : "
            << ks256gpu::kNumWords
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
            << "  闈炵粨鏋滄Ы鏈€澶х粷瀵瑰€?          : "
            << maxInactiveAbs
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
            << ", result bit=0\n";

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

        PrintComparisonProfile(
            comparisonTiming);

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
            "P-only 姣旇緝鏃嬭浆瀵嗛挜鐢熸垚",
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
            "鎺╃爜缂栫爜骞跺姞杞藉埌 CPU",
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
            "CPU 鏈€楂樹笉鍚屼綅 P-only 灏忎簬绛変簬鍒ゆ柇",
            phase.comparisonMs);

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

        const double comparisonSeconds =
            phase.comparisonMs / 1000.0;

        const double bootSeconds =
            phase.bootstrapMs / 1000.0;

        const double onlineSeconds =
            comparisonSeconds + bootSeconds;

        std::cout
            << "\n================ 鍚炲悙涓庢憡閿€ ================\n";

        std::cout
            << std::fixed
            << std::setprecision(6);

        std::cout
            << "  灏忎簬绛変簬鍒ゆ柇鎽婇攢鏃堕棿/256-bit word : "
            << comparisonSeconds /
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

        if (comparisonSeconds > 0.0) {
            std::cout
                << "  浠呭皬浜庣瓑浜庡垽鏂悶鍚愰噺              : "
                << ks256gpu::kNumWords /
                       comparisonSeconds
                << " word/s\n";
        }

        if (onlineSeconds > 0.0) {
            std::cout
                << "  灏忎簬绛変簬鍒ゆ柇+Boot 鍦ㄧ嚎鍚炲悙閲?     : "
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




