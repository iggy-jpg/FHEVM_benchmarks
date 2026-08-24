#include <openfhe.h>

#include "and_gpu.h"

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
    double andMs             = 0.0;
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

void BuildExpectedAnd(const std::vector<uint8_t>& valuesA,
                      const std::vector<uint8_t>& valuesB,
                      std::vector<uint8_t>* expected) {
    if (!expected) {
        throw std::invalid_argument(
            "BuildExpectedAnd received a null output vector");
    }

    if (valuesA.size() != ks256gpu::kNumWords ||
        valuesB.size() != ks256gpu::kNumWords) {
        throw std::invalid_argument(
            "BuildExpectedAnd received an invalid input size");
    }

    expected->assign(ks256gpu::kNumWords, 0);

    for (uint32_t word = 0; word < ks256gpu::kNumWords; ++word) {
        (*expected)[word] =
            static_cast<uint8_t>(valuesA[word] & valuesB[word]);
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

void PrintAndProfile(const ks256gpu::TimingReport& profile) {
    std::cout
        << "\n================ CPU Boolean AND 鍐呴儴璁℃椂 ================\n";
    PrintPhase("AND = A * B", profile.productMs);
    PrintPhase("EvalAnd 鎬绘椂闂?, profile.totalMs);
}

}  // namespace

int main() {
    try {
        const auto programBegin = Clock::now();
        PhaseTimes phase;

        constexpr uint32_t ringDim = 1U << 17;  // 131072
        constexpr uint32_t numSlots = ks256gpu::kTotalSlots;  // 65536
        constexpr uint32_t multiplicativeDepth = 25;
        constexpr uint64_t randomSeed = 0x414E443235364750ULL;

        const SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
        const std::vector<uint32_t> levelBudget = {3, 3};
        const std::vector<uint32_t> bsgsDim = {0, 0};

        std::cout << "============================================================\n";
        std::cout << "OpenFHE CPU: " << ks256gpu::kNumWords
                  << " x Boolean AND (xy) + 涓€娆?Bootstrap\n";
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
        std::cout << "  Parallel Boolean ANDs       : "
                  << ks256gpu::kNumWords << "\n";
        std::cout << "  Slots per Boolean value     : "
                  << ks256gpu::kStride
                  << " (bit 0 active + 511 inactive)\n";
        std::cout << "  Result slot per word        : bit 0\n";
        std::cout << "  AND formula                 : x * y\n";
        std::cout << "  Ciphertext multiplications  : 1\n";
        std::cout << "  AND ciphertext mult depth   : 1\n";
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

        // AND uses no rotations. Keep the phase entry for output alignment
        // with the other logical operators.
        phase.rotateKeyGenMs = 0.0;

        begin = Clock::now();
        cc->EvalBootstrapSetup(levelBudget, bsgsDim, numSlots, 0);
        phase.bootstrapSetupMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        cc->EvalBootstrapKeyGen(keys.secretKey, numSlots);
        phase.bootstrapKeyGenMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
phase.loadContextMs = Milliseconds(begin, Clock::now());

        // AND needs no plaintext masks. Keep the phase entry for output
        // alignment with the other logical operators.
        phase.maskBuildLoadMs = 0.0;

        // --------------------------------------------------------------
        // 3. Generate 128 independent Boolean input pairs. Only bit 0 of
        // every packed word is active, matching comparison-operator output.
        // --------------------------------------------------------------
        begin = Clock::now();

        std::vector<uint8_t> valuesA(ks256gpu::kNumWords, 0);
        std::vector<uint8_t> valuesB(ks256gpu::kNumWords, 0);
        std::vector<double> slotsA(ks256gpu::kTotalSlots, 0.0);
        std::vector<double> slotsB(ks256gpu::kTotalSlots, 0.0);

        std::mt19937_64 rng(randomSeed);
        for (uint32_t word = 0; word < ks256gpu::kNumWords; ++word) {
            valuesA[word] = static_cast<uint8_t>(rng() & 1ULL);
            valuesB[word] = static_cast<uint8_t>(rng() & 1ULL);
        }

        // Explicitly cover all four Boolean input combinations.
        valuesA[0] = 0;
        valuesB[0] = 0;
        valuesA[1] = 0;
        valuesB[1] = 1;
        valuesA[2] = 1;
        valuesB[2] = 0;
        valuesA[3] = 1;
        valuesB[3] = 1;

        for (uint32_t word = 0; word < ks256gpu::kNumWords; ++word) {
            const uint32_t base = word * ks256gpu::kStride;
            slotsA[base] = static_cast<double>(valuesA[word]);
            slotsB[base] = static_cast<double>(valuesB[word]);
        }

        std::vector<uint8_t> expected;
        BuildExpectedAnd(valuesA, valuesB, &expected);
        phase.inputBuildMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        Plaintext pA = cc->MakeCKKSPackedPlaintext(
            slotsA, 1, 0, nullptr, numSlots);
        Plaintext pB = cc->MakeCKKSPackedPlaintext(
            slotsB, 1, 0, nullptr, numSlots);
        pA->SetLength(numSlots);
        pB->SetLength(numSlots);
        phase.encodingMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        auto cA = cc->Encrypt(keys.publicKey, pA);
        auto cB = cc->Encrypt(keys.publicKey, pB);
phase.encryptionMs = Milliseconds(begin, Clock::now());

        std::cout << "\n杈撳叆鏍蜂緥:\n";
        for (uint32_t word = 0; word < 4; ++word) {
            std::cout << "  word " << word
                      << " A/B/AND = "
                      << static_cast<uint32_t>(valuesA[word]) << " / "
                      << static_cast<uint32_t>(valuesB[word]) << " / "
                      << static_cast<uint32_t>(expected[word]) << "\n";
        }

        // --------------------------------------------------------------
        // 4. CPU Boolean AND circuit: A * B.
        // --------------------------------------------------------------
        ks256gpu::TimingReport andTiming;
begin = Clock::now();
        auto cAnd = ks256gpu::EvalAnd(cc, cA, cB, &andTiming);
phase.andMs = Milliseconds(begin, Clock::now());

        std::cout << "\nBoolean AND 鍚庡瘑鏂囩姸鎬?\n";
        std::cout << "  GetLevel()         : " << cAnd->GetLevel() << "\n";
        std::cout << "  GetNoiseScaleDeg() : "
                  << cAnd->GetNoiseScaleDeg() << "\n";
        std::cout << "  浼扮畻鍓╀綑灞傛暟      : "
                  << LevelsRemaining(multiplicativeDepth, cAnd) << "\n";

        // --------------------------------------------------------------
        // 5. Exactly one CPU bootstrap after the AND circuit.
        // --------------------------------------------------------------
begin = Clock::now();
        auto cBoot = cc->EvalBootstrap(cAnd);
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
        // 7. Recover all 128 Boolean result slots and scan CKKS error.
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

        for (uint32_t word = 0; word < ks256gpu::kNumWords; ++word) {
            const uint32_t base = word * ks256gpu::kStride;
            const double target = static_cast<double>(expected[word]);
            const double error = std::abs(values[base] - target);

            sumAbsError += error;
            sumSquaredError += error * error;

            if (error > maxAbsError) {
                maxAbsError = error;
                worstWord = word;
                worstValue = values[base];
                worstExpected = target;
            }

            const uint8_t recovered = values[base] >= 0.5 ? 1U : 0U;
            if (recovered != expected[word]) {
                ++bitFailures;
            }

            for (uint32_t offset = 1; offset < ks256gpu::kStride; ++offset) {
                maxInactiveAbs = std::max(
                    maxInactiveAbs,
                    std::abs(values[base + offset]));
            }
        }

        constexpr uint64_t activeSlots = ks256gpu::kNumWords;

        const double meanAbsError =
            sumAbsError / static_cast<double>(activeSlots);
        const double rmse = std::sqrt(
            sumSquaredError / static_cast<double>(activeSlots));

        phase.verificationMs = Milliseconds(begin, Clock::now());
        phase.totalMs = Milliseconds(programBegin, Clock::now());

        std::cout << "\n================ 瑙ｅ瘑姝ｇ‘鎬т笌缁忛獙绮惧害 ================\n";
        std::cout << std::scientific << std::setprecision(10);
        std::cout << "  鏈夋晥 AND 缁撴灉妲芥暟            : "
                  << activeSlots << "\n";
        std::cout << "  鏈€澶х粷瀵硅宸?                : "
                  << maxAbsError << "\n";
        std::cout << "  骞冲潎 L1 璇樊                 : "
                  << meanAbsError << "\n";
        std::cout << "  RMSE                         : " << rmse << "\n";
        std::cout << "  闈炵粨鏋滄Ы鏈€澶х粷瀵瑰€?          : "
                  << maxInactiveAbs << "\n";

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
                  << worstWord << ", bit=0\n";

        std::cout << std::scientific << std::setprecision(12);
        std::cout << "  鏈€宸綅缃湡鏈?瑙ｅ瘑鍊?         : "
                  << worstExpected << " / " << worstValue << "\n";
        std::cout << std::fixed;
        std::cout << "  鏁翠綋缁撴灉                     : "
                  << (bitFailures == 0 ? "PASS" : "FAIL") << "\n";
        std::cout << "璇存槑锛氫笂杩?bit 鏁版槸 CKKS 瑙ｅ瘑璇樊瀵瑰簲鐨勭粡楠岀簿搴︼紝"
                     "涓嶆槸 BFV/BGV 鐨勭簿纭櫔澹伴绠椼€俓n";

        PrintAndProfile(andTiming);

        std::cout << "\n================ 绔埌绔椂闂村垎鏋?================\n";
        PrintPhase("CryptoContext 鐢熸垚涓?Enable", phase.contextMs);
        PrintPhase("KeyGen", phase.keyGenMs);
        PrintPhase("EvalMultKeyGen", phase.multKeyGenMs);
        PrintPhase("EvalRotateKeyGen锛圓ND 鏃犻渶锛?, phase.rotateKeyGenMs);
        PrintPhase("EvalBootstrapSetup", phase.bootstrapSetupMs);
        PrintPhase("EvalBootstrapKeyGen", phase.bootstrapKeyGenMs);
        PrintPhase("LoadContext 涓庤瘎浼板瘑閽ュ姞杞藉埌 CPU", phase.loadContextMs);
        PrintPhase("鎺╃爜缂栫爜骞跺姞杞藉埌 GPU锛圓ND 鏃犻渶锛?, phase.maskBuildLoadMs);
        PrintPhase("杈撳叆鐢熸垚涓?CPU 鍙傝€冪粨鏋?, phase.inputBuildMs);
        PrintPhase("CKKS 缂栫爜", phase.encodingMs);
        PrintPhase("杈撳叆鍔犲瘑骞跺姞杞藉埌 CPU", phase.encryptionMs);
        PrintPhase("CPU Boolean AND", phase.andMs);
        PrintPhase("CPU 涓€娆?EvalBootstrap", phase.bootstrapMs);
        PrintPhase("澶嶅埗鍥?CPU 骞惰В瀵?, phase.decryptionMs);
        PrintPhase("姝ｇ‘鎬т笌璇樊鎵弿", phase.verificationMs);
        PrintPhase("绋嬪簭鎬诲閽熸椂闂?, phase.totalMs);

        const double andSeconds = phase.andMs / 1000.0;
        const double bootSeconds = phase.bootstrapMs / 1000.0;
        const double onlineSeconds = andSeconds + bootSeconds;

        std::cout << "\n================ 鍚炲悙涓庢憡閿€ ================\n";
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "  AND 鎽婇攢鏃堕棿/Boolean value   : "
                  << andSeconds / ks256gpu::kNumWords << " s\n";
        std::cout << "  Boot 鎽婇攢鏃堕棿/value          : "
                  << bootSeconds / ks256gpu::kNumWords << " s\n";
        std::cout << "  鍦ㄧ嚎鎬绘憡閿€鏃堕棿/value         : "
                  << onlineSeconds / ks256gpu::kNumWords << " s\n";

        if (andSeconds > 0.0) {
            std::cout << "  浠?AND 鍚炲悙閲?               : "
                      << ks256gpu::kNumWords / andSeconds
                      << " value/s\n";
        }
        if (onlineSeconds > 0.0) {
            std::cout << "  AND+Boot 鍦ㄧ嚎鍚炲悙閲?         : "
                      << ks256gpu::kNumWords / onlineSeconds
                      << " value/s\n";
        }

        return bitFailures == 0 ? 0 : 2;
    }
    catch (const std::exception& error) {
        std::cerr << "\n[FATAL] " << error.what() << "\n";
        return 1;
    }
}




