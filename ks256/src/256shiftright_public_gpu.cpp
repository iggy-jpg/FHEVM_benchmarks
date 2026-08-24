#include <openfhe.h>

#include "shiftright_public_gpu.h"

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
    double rotateKeyGenMs    = 0.0;
    double loadContextMs     = 0.0;
    double maskBuildLoadMs   = 0.0;
    double inputBuildMs      = 0.0;
    double encodingMs        = 0.0;
    double encryptionMs      = 0.0;
    double shiftMs           = 0.0;
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

/**
 * Parses a nonnegative decimal Solidity-style public shift amount.
 *
 * The operator only needs to distinguish values 0...255 from values >=256,
 * so any decimal value that reaches 256 is saturated to 256. This also lets
 * the test accept decimal values much larger than uint64_t.
 */
uint32_t ParsePublicShift(int argc, char** argv) {
    constexpr uint32_t kDefaultShift = 13;

    if (argc <= 1) {
        return kDefaultShift;
    }

    const std::string text(argv[1]);

    if (text.empty()) {
        throw std::invalid_argument(
            "public shift must be a nonnegative decimal integer");
    }

    uint32_t value = 0;

    for (char ch : text) {
        if (ch < '0' || ch > '9') {
            throw std::invalid_argument(
                "public shift must contain decimal digits only");
        }

        const uint32_t digit =
            static_cast<uint32_t>(ch - '0');

        if (value >= ks256gpu::kWordBits) {
            return ks256gpu::kWordBits;
        }

        if (value >
            (ks256gpu::kWordBits - digit) / 10U) {
            return ks256gpu::kWordBits;
        }

        value = value * 10U + digit;

        if (value >= ks256gpu::kWordBits) {
            return ks256gpu::kWordBits;
        }
    }

    return value;
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

void BuildExpectedShiftRight(
    const std::vector<uint8_t>& input,
    uint32_t shift,
    std::vector<uint8_t>* expected) {

    expected->assign(
        ks256gpu::kTotalSlots,
        0);

    if (shift >= ks256gpu::kWordBits) {
        return;
    }

    for (uint32_t word = 0;
         word < ks256gpu::kNumWords;
         ++word) {
        const uint32_t base =
            word * ks256gpu::kStride;

        for (uint32_t bit = 0;
             bit < ks256gpu::kWordBits - shift;
             ++bit) {
            (*expected)[base + bit] =
                input[base + bit + shift];
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
    const ks256gpu::TimingReport& profile,
    uint32_t shift) {

    std::cout
        << "\n================ CPU 鏄庢枃绉讳綅閲忓彸绉诲唴閮ㄨ鏃?================\n";

    if (shift == 0U) {
        PrintPhase(
            "shift=0 鐩存帴澶嶇敤杈撳叆瀵嗘枃",
            profile.totalMs);
    }
    else if (shift >= ks256gpu::kWordBits) {
        PrintPhase(
            "shift>=256 鐢熸垚闆跺瘑鏂?,
            profile.zeroMs);
    }
    else {
        PrintPhase(
            "鍥哄畾璺濈 CKKS 鏃嬭浆",
            profile.rotateMs);

        PrintPhase(
            "鍙崇Щ杈圭晫鏄庢枃鎺╃爜涔樻硶",
            profile.maskMs);
    }

    PrintPhase(
        "EvalShiftRightPublic 鎬绘椂闂?,
        profile.totalMs);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const auto programBegin =
            Clock::now();

        PhaseTimes phase;

        const uint32_t publicShift =
            ParsePublicShift(argc, argv);

        constexpr uint32_t ringDim =
            1U << 17;  // 131072

        constexpr uint32_t numSlots =
            ks256gpu::kTotalSlots;  // 65536

        constexpr uint32_t multiplicativeDepth = 2;

        constexpr uint64_t randomSeed =
            0x5348525055423235ULL;

        const SecretKeyDist secretKeyDist =
            UNIFORM_TERNARY;


        std::cout
            << "============================================================\n";

        std::cout
            << "OpenFHE CPU: "
            << ks256gpu::kNumWords
            << " x uint256 Solidity 鏄庢枃绉讳綅閲忓彸绉伙紙鏃?Bootstrap锛塡n";

        std::cout
            << "============================================================\n";

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
            << "  Parallel uint256 shifts     : "
            << ks256gpu::kNumWords
            << "\n";

        std::cout
            << "  Slots per shift             : "
            << ks256gpu::kStride
            << " (256 active + 256 guard)\n";

        std::cout
            << "  Public shift amount         : ";

        if (publicShift >= ks256gpu::kWordBits) {
            std::cout << ">=256 (saturated to zero-result branch)\n";
        }
        else {
            std::cout << publicShift << "\n";
        }

        std::cout
            << "  Solidity semantics          : uint256 x >> shift\n";

        std::cout
            << "  Ciphertext multiplications  : 0\n";

        std::cout
            << "  Required CKKS rotations     : "
            << ((publicShift > 0U &&
                 publicShift < ks256gpu::kWordBits)
                    ? 1
                    : 0)
            << "\n";

        std::cout
            << "  Configured mult depth       : "
            << multiplicativeDepth
            << "\n";

        std::cout
            << "  Bootstrap                   : disabled\n";

        std::cout
            << "  OpenFHE visible devices    : 0\n";

        // --------------------------------------------------------------
        // Generate only the evaluation keys required by the public rotation.
        // For one public shift, generate only its one direct rotation key.
        // --------------------------------------------------------------
        begin = Clock::now();

        auto keys =
            cc->KeyGen();

        phase.keyGenMs =
            Milliseconds(
                begin,
                Clock::now());


        if (publicShift > 0U &&
            publicShift < ks256gpu::kWordBits) {
            const std::vector<int32_t> rotationIndices = {
                static_cast<int32_t>(publicShift)
            };

            begin = Clock::now();

            cc->EvalRotateKeyGen(
                keys.secretKey,
                rotationIndices);

            phase.rotateKeyGenMs =
                Milliseconds(
                    begin,
                    Clock::now());
        }


        begin = Clock::now();
phase.loadContextMs =
            Milliseconds(
                begin,
                Clock::now());

        begin = Clock::now();

        auto masks =
            ks256gpu::BuildAndLoadMasks(
                cc,
                publicShift);

        phase.maskBuildLoadMs =
            Milliseconds(
                begin,
                Clock::now());

        // --------------------------------------------------------------
        // Build 128 independent uint256 Boolean inputs.
        // --------------------------------------------------------------
        begin = Clock::now();

        std::vector<uint8_t> inputBits(
            ks256gpu::kTotalSlots,
            0);

        std::vector<double> inputSlots(
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
                const uint8_t value =
                    static_cast<uint8_t>(
                        rng() & 1ULL);

                inputBits[base + bit] =
                    value;

                inputSlots[base + bit] =
                    static_cast<double>(value);
            }
        }

        // Test 0: all ones, exercises truncation of low bits.
        {
            const uint32_t base = 0;

            for (uint32_t bit = 0;
                 bit < ks256gpu::kWordBits;
                 ++bit) {
                inputBits[base + bit] = 1U;
                inputSlots[base + bit] = 1.0;
            }
        }

        // Test 1: only bit 0 is one.
        {
            const uint32_t base =
                ks256gpu::kStride;

            for (uint32_t bit = 0;
                 bit < ks256gpu::kWordBits;
                 ++bit) {
                inputBits[base + bit] = 0U;
                inputSlots[base + bit] = 0.0;
            }

            inputBits[base] = 1U;
            inputSlots[base] = 1.0;
        }

        // Test 2: only bit 255 is one.
        {
            const uint32_t base =
                2U * ks256gpu::kStride;

            for (uint32_t bit = 0;
                 bit < ks256gpu::kWordBits;
                 ++bit) {
                inputBits[base + bit] = 0U;
                inputSlots[base + bit] = 0.0;
            }

            const uint32_t msb =
                base + ks256gpu::kWordBits - 1U;

            inputBits[msb] = 1U;
            inputSlots[msb] = 1.0;
        }

        // Test 3: alternating bits.
        {
            const uint32_t base =
                3U * ks256gpu::kStride;

            for (uint32_t bit = 0;
                 bit < ks256gpu::kWordBits;
                 ++bit) {
                const uint8_t value =
                    static_cast<uint8_t>(
                        bit & 1U);

                inputBits[base + bit] =
                    value;

                inputSlots[base + bit] =
                    static_cast<double>(value);
            }
        }

        std::vector<uint8_t> expected;

        BuildExpectedShiftRight(
            inputBits,
            publicShift,
            &expected);

        phase.inputBuildMs =
            Milliseconds(
                begin,
                Clock::now());

        begin = Clock::now();

        Plaintext plaintextInput =
            cc->MakeCKKSPackedPlaintext(
                inputSlots,
                1,
                0,
                nullptr,
                numSlots);

        plaintextInput->SetLength(
            numSlots);

        phase.encodingMs =
            Milliseconds(
                begin,
                Clock::now());

        begin = Clock::now();

        auto ciphertextInput =
            cc->Encrypt(
                keys.publicKey,
                plaintextInput);
phase.encryptionMs =
            Milliseconds(
                begin,
                Clock::now());

        std::cout << "\n杈撳叆鏍蜂緥:\n";

        std::cout
            << "  word 0 X       = "
            << BitsToHex(inputBits, 0)
            << "\n";

        std::cout
            << "  word 0 X>>shift = "
            << BitsToHex(expected, 0)
            << "\n";

        std::cout
            << "  word 1 X       = "
            << BitsToHex(inputBits, 1)
            << "\n";

        std::cout
            << "  word 1 X>>shift = "
            << BitsToHex(expected, 1)
            << "\n";

        // --------------------------------------------------------------
        // CPU public-shift logical right-shift operator.
        // --------------------------------------------------------------
        ks256gpu::TimingReport shiftTiming;
begin = Clock::now();

        auto ciphertextShifted =
            ks256gpu::EvalShiftRightPublic(
                cc,
                ciphertextInput,
                publicShift,
                masks,
                &shiftTiming);
phase.shiftMs =
            Milliseconds(
                begin,
                Clock::now());

        std::cout
            << "\n鍙崇Щ鍚庡瘑鏂囩姸鎬?\n";

        std::cout
            << "  GetLevel()         : "
            << ciphertextShifted->GetLevel()
            << "\n";

        std::cout
            << "  GetNoiseScaleDeg() : "
            << ciphertextShifted->GetNoiseScaleDeg()
            << "\n";

        std::cout
            << "  浼扮畻鍓╀綑灞傛暟      : "
            << LevelsRemaining(
                   multiplicativeDepth,
                   ciphertextShifted)
            << "\n";

        // --------------------------------------------------------------
        // No bootstrap: decrypt the right-shift result directly.
        // --------------------------------------------------------------
        // --------------------------------------------------------------
        // Copy to CPU and decrypt the direct shift result.
        // --------------------------------------------------------------
        begin = Clock::now();

        Plaintext decrypted;

        const auto decryptResult =
            cc->Decrypt(
                keys.secretKey,
                ciphertextShifted,
                &decrypted);

        if (!decryptResult.isValid) {
            throw std::runtime_error(
                "Decryption returned an invalid result");
        }

        decrypted->SetLength(
            numSlots);

        const std::vector<double> values =
            decrypted->GetRealPackedValue();

        phase.decryptionMs =
            Milliseconds(
                begin,
                Clock::now());

        // --------------------------------------------------------------
        // Verify all active bits and all guard slots.
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
            shiftTiming,
            publicShift);

        std::cout
            << "\n================ 绔埌绔椂闂村垎鏋?================\n";

        PrintPhase(
            "CryptoContext 鐢熸垚涓?Enable",
            phase.contextMs);

        PrintPhase(
            "KeyGen",
            phase.keyGenMs);


        PrintPhase(
            "鍏紑绉讳綅閲忔棆杞瘑閽ョ敓鎴?,
            phase.rotateKeyGenMs);


        PrintPhase(
            "LoadContext 涓庤瘎浼板瘑閽ュ姞杞藉埌 CPU",
            phase.loadContextMs);

        PrintPhase(
            "鍙崇Щ鎺╃爜缂栫爜骞跺姞杞藉埌 CPU",
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
            "CPU 鏄庢枃绉讳綅閲忓彸绉?,
            phase.shiftMs);


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

        const double onlineSeconds = shiftSeconds;

        std::cout
            << "\n================ 鍚炲悙涓庢憡閿€ ================\n";

        std::cout
            << std::fixed
            << std::setprecision(6);

        std::cout
            << "  鍙崇Щ鎽婇攢鏃堕棿/uint256 word : "
            << shiftSeconds /
                   ks256gpu::kNumWords
            << " s\n";


        std::cout
            << "  鍦ㄧ嚎鍙崇Щ鎽婇攢鏃堕棿/word     : "
            << onlineSeconds /
                   ks256gpu::kNumWords
            << " s\n";

        if (shiftSeconds > 0.0) {
            std::cout
                << "  浠呭彸绉诲悶鍚愰噺              : "
                << ks256gpu::kNumWords /
                       shiftSeconds
                << " word/s\n";
        }

        if (onlineSeconds > 0.0) {
            std::cout
                << "  鍙崇Щ鍦ㄧ嚎鍚炲悙閲?           : "
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




