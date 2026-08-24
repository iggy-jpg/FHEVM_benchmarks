#include <openfhe.h>

#include "bitwise_not_gpu.h"

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
    double operationMs       = 0.0;
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
    std::string output(ks256gpu::kWordBits / 4U, '0');
    const uint32_t base = wordIndex * ks256gpu::kStride;

    for (uint32_t nibble = 0;
         nibble < ks256gpu::kWordBits / 4U;
         ++nibble) {
        uint8_t value = 0;
        for (uint32_t offset = 0; offset < 4U; ++offset) {
            const uint32_t bit = nibble * 4U + offset;
            value |= static_cast<uint8_t>(bits[base + bit] << offset);
        }
        output[output.size() - 1U - nibble] = kHex[value];
    }
    return "0x" + output;
}

void BuildExpectedBitwiseNot(const std::vector<uint8_t>& bitsA,
                             std::vector<uint8_t>* expected) {
    if (!expected) {
        throw std::invalid_argument(
            "BuildExpectedBitwiseNot received a null output vector");
    }
    if (bitsA.size() != ks256gpu::kTotalSlots) {
        throw std::invalid_argument(
            "BuildExpectedBitwiseNot received an invalid input size");
    }

    expected->assign(ks256gpu::kTotalSlots, 0);
    for (uint32_t word = 0; word < ks256gpu::kNumWords; ++word) {
        const uint32_t base = word * ks256gpu::kStride;
        for (uint32_t bit = 0; bit < ks256gpu::kWordBits; ++bit) {
            (*expected)[base + bit] =
                static_cast<uint8_t>(1U - bitsA[base + bit]);
        }
    }
}

void PrintPhase(const std::string& name, double milliseconds) {
    std::cout << "  " << std::left << std::setw(40) << name
              << std::right << std::setw(14)
              << std::fixed << std::setprecision(3)
              << milliseconds << " ms"
              << "  (" << std::setprecision(6)
              << milliseconds / 1000.0 << " s)\n";
}

void PrintBitwiseNotProfile(const ks256gpu::TimingReport& profile) {
    std::cout
        << "\n================ CPU 256-bit 鎸変綅 NOT 鍐呴儴璁℃椂 ================\n";
    PrintPhase("NOT = activeMask - A", profile.subtractMs);
    PrintPhase("EvalBitwiseNot 鎬绘椂闂?, profile.totalMs);
}

}  // namespace

int main() {
    try {
        const auto programBegin = Clock::now();
        PhaseTimes phase;

        constexpr uint32_t ringDim = 1U << 17;  // 131072
        constexpr uint32_t numSlots = ks256gpu::kTotalSlots;  // 65536
        constexpr uint32_t multiplicativeDepth = 25;
        constexpr uint64_t randomSeed = 0x424E4F5432353647ULL;

        const SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
        const std::vector<uint32_t> levelBudget = {3, 3};
        const std::vector<uint32_t> bsgsDim = {0, 0};

        std::cout << "============================================================\n";
        std::cout << "OpenFHE CPU: " << ks256gpu::kNumWords
                  << " x 256-bit 鎸変綅 NOT (1-x) + 涓€娆?Bootstrap\n";
        std::cout << "============================================================\n";

        // --------------------------------------------------------------
        // 1. Create a multi-CPU OpenFHE CKKS context.
        // --------------------------------------------------------------
        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetSecretKeyDist(secretKeyDist);
        parameters.SetSecurityLevel(HEStd_128_classic);
        parameters.SetRingDim(ringDim);
        parameters.SetBatchSize(numSlots);
        parameters.SetMultiplicativeDepth(multiplicativeDepth);
        parameters.SetScalingTechnique(FLEXIBLEAUTO);
        parameters.SetScalingModSize(59);
        parameters.SetFirstModSize(60);
        parameters.SetKeySwitchTechnique(HYBRID);
        parameters.SetNumLargeDigits(3);
        const std::vector<int> devices{0, 1};
auto begin = Clock::now();
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);
        cc->Enable(FHE);
        phase.contextMs = Milliseconds(begin, Clock::now());

        std::cout << "鍙傛暟:\n";
        std::cout << "  Ring dimension              : "
                  << cc->GetRingDimension() << "\n";
        std::cout << "  CKKS slots                  : " << numSlots << "\n";
        std::cout << "  Parallel 256-bit NOTs       : "
                  << ks256gpu::kNumWords << "\n";
        std::cout << "  Slots per NOT               : "
                  << ks256gpu::kStride
                  << " (256 active + 256 guard)\n";
        std::cout << "  Active result slots/word    : "
                  << ks256gpu::kWordBits << "\n";
        std::cout << "  NOT formula                 : 1 - x\n";
        std::cout << "  Ciphertext multiplications  : 0\n";
        std::cout << "  NOT ciphertext mult depth   : 0\n";
        std::cout << "  Configured mult depth       : "
                  << multiplicativeDepth << "\n";
        std::cout << "  Bootstrap level budget      : {"
                  << levelBudget[0] << ", " << levelBudget[1] << "}\n";
        std::cout << "  OpenFHE visible devices    : 0, 1\n";

        // --------------------------------------------------------------
        // 2. Generate evaluation keys before LoadContext.
        // --------------------------------------------------------------
        begin = Clock::now();
        auto keys = cc->KeyGen();
        phase.keyGenMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        cc->EvalMultKeyGen(keys.secretKey);
        phase.multKeyGenMs = Milliseconds(begin, Clock::now());

        phase.rotateKeyGenMs = 0.0;

        begin = Clock::now();
        cc->EvalBootstrapSetup(levelBudget, bsgsDim, numSlots, 0);
        phase.bootstrapSetupMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        cc->EvalBootstrapKeyGen(keys.secretKey, numSlots);
        phase.bootstrapKeyGenMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
phase.loadContextMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        auto masks = ks256gpu::BuildAndLoadBitwiseNotMasks(cc);
        phase.maskBuildLoadMs = Milliseconds(begin, Clock::now());

        // --------------------------------------------------------------
        // 3. Generate 128 independent uint256 inputs.
        // --------------------------------------------------------------
        begin = Clock::now();
        std::vector<uint8_t> bitsA(ks256gpu::kTotalSlots, 0);
        std::vector<double> slotsA(ks256gpu::kTotalSlots, 0.0);

        std::mt19937_64 rng(randomSeed);
        for (uint32_t word = 0; word < ks256gpu::kNumWords; ++word) {
            const uint32_t base = word * ks256gpu::kStride;
            for (uint32_t bit = 0; bit < ks256gpu::kWordBits; ++bit) {
                const uint8_t a = static_cast<uint8_t>(rng() & 1ULL);
                bitsA[base + bit] = a;
                slotsA[base + bit] = static_cast<double>(a);
            }
        }

        // Explicit boundary patterns.
        for (uint32_t bit = 0; bit < ks256gpu::kWordBits; ++bit) {
            bitsA[bit] = 0;
            slotsA[bit] = 0.0;

            const uint32_t base1 = ks256gpu::kStride;
            bitsA[base1 + bit] = 1;
            slotsA[base1 + bit] = 1.0;

            const uint32_t base2 = 2U * ks256gpu::kStride;
            const uint8_t alternating = static_cast<uint8_t>(bit & 1U);
            bitsA[base2 + bit] = alternating;
            slotsA[base2 + bit] = static_cast<double>(alternating);

            const uint32_t base3 = 3U * ks256gpu::kStride;
            const uint8_t inverse = static_cast<uint8_t>((bit & 1U) ^ 1U);
            bitsA[base3 + bit] = inverse;
            slotsA[base3 + bit] = static_cast<double>(inverse);
        }

        std::vector<uint8_t> expected;
        BuildExpectedBitwiseNot(bitsA, &expected);
        phase.inputBuildMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        Plaintext pA = cc->MakeCKKSPackedPlaintext(
            slotsA, 1, 0, nullptr, numSlots);
        pA->SetLength(numSlots);
        phase.encodingMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        auto cA = cc->Encrypt(keys.publicKey, pA);
phase.encryptionMs = Milliseconds(begin, Clock::now());

        std::cout << "\n杈撳叆鏍蜂緥:\n";
        std::cout << "  word 0 A       = " << BitsToHex(bitsA, 0) << "\n";
        std::cout << "  word 0 NOT     = " << BitsToHex(expected, 0) << "\n";
        std::cout << "  word 1 A       = " << BitsToHex(bitsA, 1) << "\n";
        std::cout << "  word 1 NOT     = " << BitsToHex(expected, 1) << "\n";

        // --------------------------------------------------------------
        // 4. CPU uint256 bitwise NOT circuit: activeMask - A.
        // --------------------------------------------------------------
        ks256gpu::TimingReport operationTiming;
begin = Clock::now();
        auto cResult = ks256gpu::EvalBitwiseNot(
            cc, cA, masks, &operationTiming);
phase.operationMs = Milliseconds(begin, Clock::now());

        std::cout << "\n鎸変綅 NOT 鍚庡瘑鏂囩姸鎬?\n";
        std::cout << "  GetLevel()         : " << cResult->GetLevel() << "\n";
        std::cout << "  GetNoiseScaleDeg() : "
                  << cResult->GetNoiseScaleDeg() << "\n";
        std::cout << "  浼扮畻鍓╀綑灞傛暟      : "
                  << LevelsRemaining(multiplicativeDepth, cResult) << "\n";

        // --------------------------------------------------------------
        // 5. Exactly one CPU bootstrap after the operation.
        // --------------------------------------------------------------
begin = Clock::now();
        auto cBoot = cc->EvalBootstrap(cResult);
phase.bootstrapMs = Milliseconds(begin, Clock::now());

        std::cout << "\nBootstrap 鍚庡瘑鏂囩姸鎬?\n";
        std::cout << "  GetLevel()         : " << cBoot->GetLevel() << "\n";
        std::cout << "  GetNoiseScaleDeg() : "
                  << cBoot->GetNoiseScaleDeg() << "\n";
        std::cout << "  浼扮畻鍓╀綑灞傛暟      : "
                  << LevelsRemaining(multiplicativeDepth, cBoot) << "\n";

        // --------------------------------------------------------------
        // 6. Copy to CPU and decrypt after bootstrap.
        // --------------------------------------------------------------
        begin = Clock::now();
        Plaintext decrypted;
        const auto decryptResult =
            cc->Decrypt(keys.secretKey, cBoot, &decrypted);
        if (!decryptResult.isValid) {
            throw std::runtime_error("Decryption returned an invalid result");
        }
        decrypted->SetLength(numSlots);
        const std::vector<double> values = decrypted->GetRealPackedValue();
        phase.decryptionMs = Milliseconds(begin, Clock::now());

        // --------------------------------------------------------------
        // 7. Recover all active bits and scan CKKS error.
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

        for (uint32_t word = 0; word < ks256gpu::kNumWords; ++word) {
            const uint32_t base = word * ks256gpu::kStride;

            for (uint32_t bit = 0; bit < ks256gpu::kWordBits; ++bit) {
                const uint32_t index = base + bit;
                const double target = static_cast<double>(expected[index]);
                const double error = std::abs(values[index] - target);

                sumAbsError += error;
                sumSquaredError += error * error;

                if (error > maxAbsError) {
                    maxAbsError = error;
                    worstWord = word;
                    worstBit = bit;
                    worstValue = values[index];
                    worstExpected = target;
                }

                const uint8_t recovered = values[index] >= 0.5 ? 1U : 0U;
                if (recovered != expected[index]) {
                    ++bitFailures;
                }
            }

            for (uint32_t bit = ks256gpu::kWordBits;
                 bit < ks256gpu::kStride;
                 ++bit) {
                maxGuardAbs = std::max(
                    maxGuardAbs,
                    std::abs(values[base + bit]));
            }
        }

        constexpr uint64_t activeSlots =
            static_cast<uint64_t>(ks256gpu::kNumWords) *
            ks256gpu::kWordBits;

        const double meanAbsError =
            sumAbsError / static_cast<double>(activeSlots);
        const double rmse = std::sqrt(
            sumSquaredError / static_cast<double>(activeSlots));

        phase.verificationMs = Milliseconds(begin, Clock::now());
        phase.totalMs = Milliseconds(programBegin, Clock::now());

        std::cout << "\n================ 瑙ｅ瘑姝ｇ‘鎬т笌缁忛獙绮惧害 ================\n";
        std::cout << std::scientific << std::setprecision(10);
        std::cout << "  鏈夋晥 NOT 妲芥暟                 : "
                  << activeSlots << "\n";
        std::cout << "  鏈€澶х粷瀵硅宸?                : "
                  << maxAbsError << "\n";
        std::cout << "  骞冲潎 L1 璇樊                 : "
                  << meanAbsError << "\n";
        std::cout << "  RMSE                         : " << rmse << "\n";
        std::cout << "  guard 妲芥渶澶х粷瀵瑰€?          : "
                  << maxGuardAbs << "\n";

        std::cout << std::fixed << std::setprecision(3);
        std::cout << "  GetLogPrecision              : "
                  << decrypted->GetLogPrecision() << " bits\n";
        std::cout << "  鏈€宸簿搴?-log2(max error)    : "
                  << PrecisionBits(maxAbsError) << " bits\n";
        std::cout << "  骞冲潎绮惧害 -log2(mean L1)      : "
                  << PrecisionBits(meanAbsError) << " bits\n";
        std::cout << "  RMSE 绮惧害 -log2(RMSE)        : "
                  << PrecisionBits(rmse) << " bits\n";
        std::cout << "  0.5 闃堝€兼仮澶嶅け璐ユ暟           : "
                  << bitFailures << "\n";
        std::cout << "  鏈€宸綅缃?                    : word="
                  << worstWord << ", bit=" << worstBit << "\n";

        std::cout << std::scientific << std::setprecision(12);
        std::cout << "  鏈€宸綅缃湡鏈?瑙ｅ瘑鍊?         : "
                  << worstExpected << " / " << worstValue << "\n";
        std::cout << std::fixed;
        std::cout << "  鏁翠綋缁撴灉                     : "
                  << (bitFailures == 0 ? "PASS" : "FAIL") << "\n";
        std::cout << "璇存槑锛氫笂杩?bit 鏁版槸 CKKS 瑙ｅ瘑璇樊瀵瑰簲鐨勭粡楠岀簿搴︼紝"
                     "涓嶆槸 BFV/BGV 鐨勭簿纭櫔澹伴绠椼€俓n";

        PrintBitwiseNotProfile(operationTiming);

        std::cout << "\n================ 绔埌绔椂闂村垎鏋?================\n";
        PrintPhase("CryptoContext 鐢熸垚涓?Enable", phase.contextMs);
        PrintPhase("KeyGen", phase.keyGenMs);
        PrintPhase("EvalMultKeyGen", phase.multKeyGenMs);
        PrintPhase("EvalRotateKeyGen锛圢OT 鏃犻渶锛?, phase.rotateKeyGenMs);
        PrintPhase("EvalBootstrapSetup", phase.bootstrapSetupMs);
        PrintPhase("EvalBootstrapKeyGen", phase.bootstrapKeyGenMs);
        PrintPhase("LoadContext 涓庤瘎浼板瘑閽ュ姞杞藉埌 CPU", phase.loadContextMs);
        PrintPhase("active mask 缂栫爜骞跺姞杞藉埌 CPU", phase.maskBuildLoadMs);
        PrintPhase("杈撳叆鐢熸垚涓?CPU 鍙傝€冪粨鏋?, phase.inputBuildMs);
        PrintPhase("CKKS 缂栫爜", phase.encodingMs);
        PrintPhase("杈撳叆鍔犲瘑骞跺姞杞藉埌 CPU", phase.encryptionMs);
        PrintPhase("CPU 256-bit 鎸変綅 NOT", phase.operationMs);
        PrintPhase("CPU 涓€娆?EvalBootstrap", phase.bootstrapMs);
        PrintPhase("澶嶅埗鍥?CPU 骞惰В瀵?, phase.decryptionMs);
        PrintPhase("姝ｇ‘鎬т笌璇樊鎵弿", phase.verificationMs);
        PrintPhase("绋嬪簭鎬诲閽熸椂闂?, phase.totalMs);

        const double operationSeconds = phase.operationMs / 1000.0;
        const double bootSeconds = phase.bootstrapMs / 1000.0;
        const double onlineSeconds = operationSeconds + bootSeconds;

        std::cout << "\n================ 鍚炲悙涓庢憡閿€ ================\n";
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "  NOT 鎽婇攢鏃堕棿/256-bit word    : "
                  << operationSeconds / ks256gpu::kNumWords << " s\n";
        std::cout << "  Boot 鎽婇攢鏃堕棿/word           : "
                  << bootSeconds / ks256gpu::kNumWords << " s\n";
        std::cout << "  鍦ㄧ嚎鎬绘憡閿€鏃堕棿/word          : "
                  << onlineSeconds / ks256gpu::kNumWords << " s\n";

        if (operationSeconds > 0.0) {
            std::cout << "  浠?NOT 鍚炲悙閲?               : "
                      << ks256gpu::kNumWords / operationSeconds
                      << " word/s\n";
        }
        if (onlineSeconds > 0.0) {
            std::cout << "  NOT+Boot 鍦ㄧ嚎鍚炲悙閲?         : "
                      << ks256gpu::kNumWords / onlineSeconds
                      << " word/s\n";
        }

        return bitFailures == 0 ? 0 : 2;
    }
    catch (const std::exception& error) {
        std::cerr << "\n[FATAL] " << error.what() << "\n";
        return 1;
    }
}




