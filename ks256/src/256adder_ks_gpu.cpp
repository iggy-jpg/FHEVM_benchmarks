#include <openfhe.h>

#include "adder_ks_gpu.h"

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
    double additionMs        = 0.0;
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
    std::string output(ks256gpu::kWordBits / 4, '0');
    const uint32_t base = wordIndex * ks256gpu::kStride;

    for (uint32_t nibble = 0; nibble < ks256gpu::kWordBits / 4;
         ++nibble) {
        uint8_t value = 0;
        for (uint32_t offset = 0; offset < 4; ++offset) {
            const uint32_t bit = nibble * 4 + offset;
            value |= static_cast<uint8_t>(bits[base + bit] << offset);
        }
        output[output.size() - 1 - nibble] = kHex[value];
    }
    return "0x" + output;
}

void BuildExpectedSum(const std::vector<uint8_t>& bitsA,
                      const std::vector<uint8_t>& bitsB,
                      std::vector<uint8_t>* expected) {
    expected->assign(ks256gpu::kTotalSlots, 0);

    for (uint32_t word = 0; word < ks256gpu::kNumWords; ++word) {
        const uint32_t base = word * ks256gpu::kStride;
        uint32_t carry = 0;

        for (uint32_t bit = 0; bit < ks256gpu::kWordBits; ++bit) {
            const uint32_t total =
                bitsA[base + bit] + bitsB[base + bit] + carry;
            (*expected)[base + bit] = static_cast<uint8_t>(total & 1U);
            carry = total >> 1U;
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

void PrintAdderProfile(const ks256gpu::TimingReport& profile) {
    std::cout << "\n================ CPU Kogge-Stone 鍐呴儴璁℃椂 ================\n";
    PrintPhase("棰勫鐞?g/p", profile.preprocessMs);

    double layersTotal = 0.0;
    for (uint32_t layer = 0; layer < ks256gpu::kNumLayers; ++layer) {
        const auto& item = profile.layers[layer];
        layersTotal += item.wallMs;

        std::cout << "\n  Layer " << layer
                  << "  (step=" << (1U << layer) << ")\n";
        PrintPhase("    Rotate(G)", item.rotateGMs);
        PrintPhase("    Rotate(P)", item.rotatePMs);
        PrintPhase("    P * Rotate(G)", item.termMultMs);
        if (layer + 1 < ks256gpu::kNumLayers) {
            PrintPhase("    P * Rotate(P)", item.prefixMultMs);
        }
        PrintPhase("    鏇存柊 G", item.updateGMs);
        PrintPhase("    鏈眰鎬诲閽熸椂闂?, item.wallMs);
    }

    std::cout << "\n";
    PrintPhase("8 灞傚墠缂€缃戠粶鍚堣", layersTotal);
    PrintPhase("杩涗綅鍚戦噺鏃嬭浆", profile.carryRotateMs);
    PrintPhase("Carry 鏄庢枃鎺╃爜涔樻硶", profile.carryMaskMs);
    PrintPhase("Group-G 鏄庢枃鎺╃爜涔樻硶", profile.groupMaskMs);
    PrintPhase("S=A+B+C-2G 閲嶆瀯", profile.reconstructMs);
    PrintPhase("EvalAddKS 鎬绘椂闂?, profile.totalMs);
}

}  // namespace

int main() {
    try {
        const auto programBegin = Clock::now();
        PhaseTimes phase;

        constexpr uint32_t ringDim = 1U << 17;   // 131072
        constexpr uint32_t numSlots = ks256gpu::kTotalSlots;  // 65536
        constexpr uint32_t multiplicativeDepth = 25;
        constexpr int directionSign = -1;
        constexpr uint64_t randomSeed = 0x4B53414444323536ULL;

        const SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
        const std::vector<uint32_t> levelBudget = {3, 3};
        const std::vector<uint32_t> bsgsDim = {0, 0};

        std::cout << "============================================================\n";
        std::cout << "OpenFHE CPU: " << ks256gpu::kNumWords
                  << " x 256-bit Kogge-Stone 鍔犳硶 + 涓€娆?Bootstrap\n";
        std::cout << "============================================================\n";

        // --------------------------------------------------------------
        // 1. Create the OpenFHE CPU-enabled CKKS context.
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

        // OpenFHE CPU execution does not require device selection.
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
        std::cout << "  Parallel 256-bit additions : "
                  << ks256gpu::kNumWords << "\n";
        std::cout << "  Slots per addition         : "
                  << ks256gpu::kStride
                  << " (256 active + 256 guard)\n";
        std::cout << "  Configured mult depth      : "
                  << multiplicativeDepth << "\n";
        std::cout << "  Bootstrap level budget     : {"
                  << levelBudget[0] << ", " << levelBudget[1] << "}\n";
        std::cout << "  OpenFHE visible device    : 0\n";

        // --------------------------------------------------------------
        // 2. Generate every evaluation key BEFORE LoadContext.
        // OpenFHE 2.1.x rejects bootstrap setup/key generation after the
        // context has already been loaded to a CPU.
        // --------------------------------------------------------------
        begin = Clock::now();
        auto keys = cc->KeyGen();
        phase.keyGenMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        cc->EvalMultKeyGen(keys.secretKey);
        phase.multKeyGenMs = Milliseconds(begin, Clock::now());

        std::vector<int32_t> rotationIndices;
        rotationIndices.reserve(ks256gpu::kNumLayers);
        for (uint32_t layer = 0; layer < ks256gpu::kNumLayers; ++layer) {
            rotationIndices.push_back(
                directionSign * static_cast<int32_t>(1U << layer));
        }

        begin = Clock::now();
        cc->EvalRotateKeyGen(keys.secretKey, rotationIndices);
        phase.rotateKeyGenMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        cc->EvalBootstrapSetup(levelBudget, bsgsDim, numSlots, 0);
        phase.bootstrapSetupMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        cc->EvalBootstrapKeyGen(keys.secretKey, numSlots);
        phase.bootstrapKeyGenMs = Milliseconds(begin, Clock::now());

        // Transfer context and all generated evaluation keys to CPU 0.
        begin = Clock::now();
phase.loadContextMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        auto masks = ks256gpu::BuildAndLoadMasks(cc, directionSign);
        phase.maskBuildLoadMs = Milliseconds(begin, Clock::now());

        // --------------------------------------------------------------
        // 3. Generate 128 independent 256-bit input pairs.
        // --------------------------------------------------------------
        begin = Clock::now();
        std::vector<uint8_t> bitsA(ks256gpu::kTotalSlots, 0);
        std::vector<uint8_t> bitsB(ks256gpu::kTotalSlots, 0);
        std::vector<double> slotsA(ks256gpu::kTotalSlots, 0.0);
        std::vector<double> slotsB(ks256gpu::kTotalSlots, 0.0);

        std::mt19937_64 rng(randomSeed);
        for (uint32_t word = 0; word < ks256gpu::kNumWords; ++word) {
            const uint32_t base = word * ks256gpu::kStride;
            for (uint32_t bit = 0; bit < ks256gpu::kWordBits; ++bit) {
                bitsA[base + bit] = static_cast<uint8_t>(rng() & 1ULL);
                bitsB[base + bit] = static_cast<uint8_t>(rng() & 1ULL);
                slotsA[base + bit] = static_cast<double>(bitsA[base + bit]);
                slotsB[base + bit] = static_cast<double>(bitsB[base + bit]);
            }
        }

        // Test 0: (2^256 - 1) + 1 = 0 mod 2^256.
        {
            const uint32_t base = 0;
            for (uint32_t bit = 0; bit < ks256gpu::kWordBits; ++bit) {
                bitsA[base + bit] = 1;
                bitsB[base + bit] = 0;
                slotsA[base + bit] = 1.0;
                slotsB[base + bit] = 0.0;
            }
            bitsB[base] = 1;
            slotsB[base] = 1.0;
        }

        // Test 1: (2^255 - 1) + 1 = 2^255.
        {
            const uint32_t base = ks256gpu::kStride;
            for (uint32_t bit = 0; bit < ks256gpu::kWordBits; ++bit) {
                const bool lower255 = bit < 255;
                bitsA[base + bit] = static_cast<uint8_t>(lower255);
                bitsB[base + bit] = 0;
                slotsA[base + bit] = lower255 ? 1.0 : 0.0;
                slotsB[base + bit] = 0.0;
            }
            bitsB[base] = 1;
            slotsB[base] = 1.0;
        }

        std::vector<uint8_t> expected;
        BuildExpectedSum(bitsA, bitsB, &expected);
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
        std::cout << "  word 0 A = " << BitsToHex(bitsA, 0) << "\n";
        std::cout << "  word 0 B = " << BitsToHex(bitsB, 0) << "\n";
        std::cout << "  word 0 S = " << BitsToHex(expected, 0) << "\n";
        std::cout << "  word 1 A = " << BitsToHex(bitsA, 1) << "\n";
        std::cout << "  word 1 B = " << BitsToHex(bitsB, 1) << "\n";
        std::cout << "  word 1 S = " << BitsToHex(expected, 1) << "\n";

        // --------------------------------------------------------------
        // 4. CPU Kogge-Stone addition.
        // --------------------------------------------------------------
        ks256gpu::TimingReport adderTiming;
begin = Clock::now();
        auto cSum = ks256gpu::EvalAddKS(
            cc, cA, cB, directionSign, masks, &adderTiming);
phase.additionMs = Milliseconds(begin, Clock::now());

        std::cout << "\n鍔犳硶鍚庡瘑鏂囩姸鎬?\n";
        std::cout << "  GetLevel()         : " << cSum->GetLevel() << "\n";
        std::cout << "  GetNoiseScaleDeg() : "
                  << cSum->GetNoiseScaleDeg() << "\n";
        std::cout << "  浼扮畻鍓╀綑灞傛暟      : "
                  << LevelsRemaining(multiplicativeDepth, cSum) << "\n";

        // --------------------------------------------------------------
        // 5. Exactly one CPU bootstrap after the addition.
        // --------------------------------------------------------------
begin = Clock::now();
        auto cBoot = cc->EvalBootstrap(cSum);
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
        // 7. Exact bit recovery and empirical CKKS precision analysis.
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
            static_cast<uint64_t>(ks256gpu::kNumWords) *
            ks256gpu::kWordBits;

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

                const uint8_t recovered = values[index] >= 0.5 ? 1 : 0;
                if (recovered != expected[index]) {
                    ++bitFailures;
                }
            }

            for (uint32_t guard = 0; guard < ks256gpu::kGuardBits;
                 ++guard) {
                const uint32_t index =
                    base + ks256gpu::kWordBits + guard;
                maxGuardAbs =
                    std::max(maxGuardAbs, std::abs(values[index]));
            }
        }

        const double meanAbsError =
            sumAbsError / static_cast<double>(activeCount);
        const double rmse = std::sqrt(
            sumSquaredError / static_cast<double>(activeCount));

        phase.verificationMs = Milliseconds(begin, Clock::now());
        phase.totalMs = Milliseconds(programBegin, Clock::now());

        std::cout << "\n================ 瑙ｅ瘑姝ｇ‘鎬т笌缁忛獙绮惧害 ================\n";
        std::cout << std::scientific << std::setprecision(10);
        std::cout << "  鏈夋晥 bit 妲芥暟                 : "
                  << activeCount << "\n";
        std::cout << "  鏈€澶х粷瀵硅宸?                : "
                  << maxAbsError << "\n";
        std::cout << "  骞冲潎 L1 璇樊                 : "
                  << meanAbsError << "\n";
        std::cout << "  RMSE                         : " << rmse << "\n";
        std::cout << "  Guard 妲芥渶澶х粷瀵瑰€?          : "
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

        PrintAdderProfile(adderTiming);

        std::cout << "\n================ 绔埌绔椂闂村垎鏋?================\n";
        PrintPhase("CryptoContext 鐢熸垚涓?Enable", phase.contextMs);
        PrintPhase("KeyGen", phase.keyGenMs);
        PrintPhase("EvalMultKeyGen", phase.multKeyGenMs);
        PrintPhase("Kogge-Stone 鏃嬭浆瀵嗛挜鐢熸垚", phase.rotateKeyGenMs);
        PrintPhase("EvalBootstrapSetup", phase.bootstrapSetupMs);
        PrintPhase("EvalBootstrapKeyGen", phase.bootstrapKeyGenMs);
        PrintPhase("LoadContext 涓庤瘎浼板瘑閽ュ姞杞藉埌 CPU", phase.loadContextMs);
        PrintPhase("鎺╃爜缂栫爜骞跺姞杞藉埌 CPU", phase.maskBuildLoadMs);
        PrintPhase("杈撳叆鐢熸垚涓?CPU 鍙傝€冪粨鏋?, phase.inputBuildMs);
        PrintPhase("CKKS 缂栫爜", phase.encodingMs);
        PrintPhase("杈撳叆鍔犲瘑骞跺姞杞藉埌 CPU", phase.encryptionMs);
        PrintPhase("CPU Kogge-Stone 鍔犳硶", phase.additionMs);
        PrintPhase("CPU 涓€娆?EvalBootstrap", phase.bootstrapMs);
        PrintPhase("澶嶅埗鍥?CPU 骞惰В瀵?, phase.decryptionMs);
        PrintPhase("姝ｇ‘鎬т笌璇樊鎵弿", phase.verificationMs);
        PrintPhase("绋嬪簭鎬诲閽熸椂闂?, phase.totalMs);

        const double addSeconds = phase.additionMs / 1000.0;
        const double bootSeconds = phase.bootstrapMs / 1000.0;
        const double onlineSeconds = addSeconds + bootSeconds;

        std::cout << "\n================ 鍚炲悙涓庢憡閿€ ================\n";
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "  鍔犳硶鎽婇攢鏃堕棿/256-bit word : "
                  << addSeconds / ks256gpu::kNumWords << " s\n";
        std::cout << "  Boot 鎽婇攢鏃堕棿/word        : "
                  << bootSeconds / ks256gpu::kNumWords << " s\n";
        std::cout << "  鍦ㄧ嚎鎬绘憡閿€鏃堕棿/word       : "
                  << onlineSeconds / ks256gpu::kNumWords << " s\n";

        if (addSeconds > 0.0) {
            std::cout << "  浠呭姞娉曞悶鍚愰噺              : "
                      << ks256gpu::kNumWords / addSeconds
                      << " word/s\n";
        }
        if (onlineSeconds > 0.0) {
            std::cout << "  鍔犳硶+Boot 鍦ㄧ嚎鍚炲悙閲?     : "
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




