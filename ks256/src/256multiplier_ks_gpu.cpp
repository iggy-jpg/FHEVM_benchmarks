#include <openfhe.h>

#include "multiplier_ks_gpu.h"

#include <algorithm>
#include <array>
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
    double multiplicationMs  = 0.0;
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

void BuildExpectedProduct(const std::vector<uint8_t>& bitsA,
                          const std::vector<uint8_t>& bitsB,
                          std::vector<uint8_t>* expected) {
    expected->assign(ks256gpu::kMultiplierTotalSlots, 0);

    for (uint32_t word = 0; word < ks256gpu::kMultiplierNumWords; ++word) {
        const uint32_t base = word * ks256gpu::kStride;
        std::array<uint32_t, ks256gpu::kWordBits> columns{};

        // Only terms i+j < 256 contribute to multiplication modulo 2^256.
        for (uint32_t i = 0; i < ks256gpu::kWordBits; ++i) {
            if (bitsA[base + i] == 0) {
                continue;
            }
            for (uint32_t j = 0; j + i < ks256gpu::kWordBits; ++j) {
                columns[i + j] += bitsB[base + j];
            }
        }

        uint32_t carry = 0;
        for (uint32_t bit = 0; bit < ks256gpu::kWordBits; ++bit) {
            const uint32_t total = columns[bit] + carry;
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

void PrintKoggeStoneProfile(const ks256gpu::TimingReport& profile) {
    std::cout << "\n================ 鏈€缁?Kogge-Stone 鍐呴儴璁℃椂 ================\n";
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

void PrintMultiplierProfile(
    const ks256gpu::MultiplierTimingReport& profile) {
    std::cout << "\n================ CPU Dadda 鍐呴儴璁℃椂 ================\n";
    PrintPhase("256 琛岄儴鍒嗙Н娴佸紡鐢熸垚", profile.partialProductsMs);

    double daddaWallTotal = 0.0;
    for (uint32_t stage = 0; stage < ks256gpu::kNumDaddaStages; ++stage) {
        const auto& item = profile.daddaStages[stage];
        daddaWallTotal += item.wallMs;

        std::cout << "\n  Stage " << stage
                  << "  (" << item.inputRows
                  << " -> " << item.targetRows << " rows)";
        std::cout << "\n";
        std::cout << "    Full adders                : "
                  << item.fullAdders << "\n";
        std::cout << "    Passthrough rows           : "
                  << item.passthroughRows << "\n";
        PrintPhase("    鏈眰鎬诲閽熸椂闂?, item.wallMs);
    }

    std::cout << "\n";
    PrintPhase("Dadda 鍚勯樁娈靛閽熸椂闂村悎璁?, daddaWallTotal);
    PrintPhase("鏈€缁堣緭鍑?Active 鎺╃爜", profile.outputMaskMs);
    PrintPhase("EvalMulDaddaKS 鎬绘椂闂?, profile.totalMs);

    std::cout << "\n  杩愮畻璁℃暟:\n";
    std::cout << "    杈撳叆瀵归綈鏃嬭浆              : "
              << profile.inputAlignmentRotations << "\n";
    std::cout << "    bit 骞挎挱鏃嬭浆              : "
              << profile.broadcastRotations << "\n";
    std::cout << "    bit-0 鏄庢枃閫夋嫨涔樻硶        : "
              << profile.selectorPlainMults << "\n";
    std::cout << "    閮ㄥ垎绉瘑鏂囦箻娉?           : "
              << profile.partialCipherMults << "\n";
    std::cout << "    鍘嬬缉鍣ㄥ瘑鏂囦箻娉?           : "
              << profile.compressorCipherMults << "\n";
    std::cout << "    鍘嬬缉鍣?Carry 鏃嬭浆         : "
              << profile.compressorCarryRotations << "\n";
    std::cout << "    娴佸紡鍚屾娆℃暟              : "
              << profile.streamingSynchronizations << "\n";
    std::cout << "    宄板€肩紦瀛樺瘑鏂囪            : "
              << profile.peakBufferedRows << "\n";

    PrintKoggeStoneProfile(profile.finalAdder);
}

}  // namespace

int main() {
    try {
        std::cout << std::unitbuf;
        const auto programBegin = Clock::now();
        PhaseTimes phase;

        constexpr uint32_t ringDim = 1U << 17;  // 131072
        constexpr uint32_t numSlots = ks256gpu::kMultiplierTotalSlots;

        // Security-preserving low-memory configuration.
        // Keep N=131072 and HEStd_128_classic. The current circuit is estimated
        // to consume about 39 levels before the final single bootstrap, so 42
        // gives a small safety margin without the large key-memory cost of 56.
        constexpr bool enableBootstrap = true;
        constexpr uint32_t estimatedCircuitDepth = 39U;
        constexpr uint32_t multiplicativeDepth = 40U;
        constexpr uint32_t compressorBatchSize =
            ks256gpu::kDefaultCompressorBatchSize;
        constexpr int directionSign = -1;
        constexpr uint64_t randomSeed = 0x44414444414D554CULL;

        const SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
        const std::vector<uint32_t> levelBudget = {4, 4};
        const std::vector<uint32_t> bsgsDim = {0, 0};

        std::cout << "============================================================\n";
        std::cout << "OpenFHE CPU: " << ks256gpu::kMultiplierNumWords
                  << " x 256-bit 娴佸紡 Dadda + Kogge-Stone 涔樻硶"
                  << "锛圠OW-MEM v6锛?2璺紝鐙珛澶存枃浠讹紝HEStd-128锛屼竴娆?Bootstrap锛塡n";
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
std::cout << "\n[setup] GenCryptoContext..." << std::endl;
        auto begin = Clock::now();
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        std::cout << "[setup] GenCryptoContext 瀹屾垚" << std::endl;
        std::cout << "[setup] Enable(PKE/KEYSWITCH/LEVELEDSHE/ADVANCEDSHE)..." << std::endl;
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);
        if constexpr (enableBootstrap) {
            cc->Enable(FHE);
        }
        phase.contextMs = Milliseconds(begin, Clock::now());
        std::cout << "[setup] Enable 瀹屾垚" << std::endl;

        std::cout << "鍙傛暟:\n";
        std::cout << "  Build tag                      : lowmem-v6-standalone-r131072-d42-128bit-boot-32words\n";
        std::cout << "  Ring dimension                 : "
                  << cc->GetRingDimension() << "\n";
        std::cout << "  CKKS slots                     : " << numSlots << "\n";
        std::cout << "  Parallel 256-bit multiplications: "
                  << ks256gpu::kMultiplierNumWords << "\n";
        std::cout << "  Slots per multiplication      : "
                  << ks256gpu::kStride
                  << " (256 active + 256 guard)\n";
        std::cout << "  Product semantics              : A * B mod 2^256\n";
        std::cout << "  Dadda target heights           : ";
        for (uint32_t i = 0; i < ks256gpu::kNumDaddaStages; ++i) {
            if (i != 0) {
                std::cout << ",";
            }
            std::cout << ks256gpu::kDaddaTargets[i];
        }
        std::cout << "\n";
        std::cout << "  Estimated circuit depth       : "
                  << estimatedCircuitDepth << "\n";
        std::cout << "  Configured mult depth         : "
                  << multiplicativeDepth << "\n";
        std::cout << "  Streaming compressor batch    : "
                  << compressorBatchSize << "\n";
        std::cout << "  Bootstrap enabled             : yes\n";
        std::cout << "  HYBRID large digits           : 2\n";
        std::cout << "  OpenFHE visible device       : 0\n";

        // --------------------------------------------------------------
        // 2. Generate every evaluation key BEFORE LoadContext.
        // --------------------------------------------------------------
        std::cout << "\n[setup] KeyGen..." << std::endl;
        begin = Clock::now();
        auto keys = cc->KeyGen();
        phase.keyGenMs = Milliseconds(begin, Clock::now());
        std::cout << "[setup] KeyGen 瀹屾垚" << std::endl;

        std::cout << "[setup] EvalMultKeyGen..." << std::endl;
        begin = Clock::now();
        cc->EvalMultKeyGen(keys.secretKey);
        phase.multKeyGenMs = Milliseconds(begin, Clock::now());
        std::cout << "[setup] EvalMultKeyGen 瀹屾垚" << std::endl;

        // Required rotations:
        //   -1,-2,...,-128 for broadcast and Kogge-Stone;
        //   +1 only for advancing B one source bit at a time.
        // Positive +2,...,+128 keys are unnecessary and would consume memory.
        std::vector<int32_t> rotationIndices;
        rotationIndices.reserve(ks256gpu::kNumLayers + 1);
        rotationIndices.push_back(1);
        for (uint32_t layer = 0; layer < ks256gpu::kNumLayers; ++layer) {
            rotationIndices.push_back(
                -static_cast<int32_t>(1U << layer));
        }

        std::cout << "[setup] EvalRotateKeyGen锛? 涓棆杞敭锛?.." << std::endl;
        begin = Clock::now();
        cc->EvalRotateKeyGen(keys.secretKey, rotationIndices);
        phase.rotateKeyGenMs = Milliseconds(begin, Clock::now());
        std::cout << "[setup] EvalRotateKeyGen 瀹屾垚" << std::endl;

        if constexpr (enableBootstrap) {
            std::cout << "[setup] EvalBootstrapSetup...\n";
            begin = Clock::now();
            cc->EvalBootstrapSetup(levelBudget, bsgsDim, numSlots, 0);
            phase.bootstrapSetupMs = Milliseconds(begin, Clock::now());
            std::cout << "[setup] EvalBootstrapSetup 瀹屾垚\n";

            std::cout << "[setup] EvalBootstrapKeyGen...\n";
            begin = Clock::now();
            cc->EvalBootstrapKeyGen(keys.secretKey, numSlots);
            phase.bootstrapKeyGenMs = Milliseconds(begin, Clock::now());
            std::cout << "[setup] EvalBootstrapKeyGen 瀹屾垚\n";
        }

        std::cout << "[setup] LoadContext..." << std::endl;
        begin = Clock::now();
phase.loadContextMs = Milliseconds(begin, Clock::now());
        std::cout << "[setup] LoadContext 瀹屾垚" << std::endl;

        std::cout << "[setup] 鏋勫缓骞跺姞杞芥帺鐮?.." << std::endl;
        begin = Clock::now();
        auto masks =
            ks256gpu::BuildAndLoadMultiplierMasks(cc, directionSign);
        phase.maskBuildLoadMs = Milliseconds(begin, Clock::now());
        std::cout << "[setup] 鎺╃爜鍔犺浇瀹屾垚" << std::endl;

        // --------------------------------------------------------------
        // 3. Generate independent 256-bit input pairs.
        // --------------------------------------------------------------
        begin = Clock::now();
        std::vector<uint8_t> bitsA(ks256gpu::kMultiplierTotalSlots, 0);
        std::vector<uint8_t> bitsB(ks256gpu::kMultiplierTotalSlots, 0);
        std::vector<double> slotsA(ks256gpu::kMultiplierTotalSlots, 0.0);
        std::vector<double> slotsB(ks256gpu::kMultiplierTotalSlots, 0.0);

        std::mt19937_64 rng(randomSeed);
        for (uint32_t word = 0; word < ks256gpu::kMultiplierNumWords; ++word) {
            const uint32_t base = word * ks256gpu::kStride;
            for (uint32_t bit = 0; bit < ks256gpu::kWordBits; ++bit) {
                bitsA[base + bit] = static_cast<uint8_t>(rng() & 1ULL);
                bitsB[base + bit] = static_cast<uint8_t>(rng() & 1ULL);
                slotsA[base + bit] = static_cast<double>(bitsA[base + bit]);
                slotsB[base + bit] = static_cast<double>(bitsB[base + bit]);
            }
        }

        // Test 0: (2^256 - 1) * 1 = 2^256 - 1.
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

        // Test 1: 2^128 * 2^127 = 2^255.
        {
            const uint32_t base = ks256gpu::kStride;
            for (uint32_t bit = 0; bit < ks256gpu::kWordBits; ++bit) {
                bitsA[base + bit] = 0;
                bitsB[base + bit] = 0;
                slotsA[base + bit] = 0.0;
                slotsB[base + bit] = 0.0;
            }
            bitsA[base + 128] = 1;
            bitsB[base + 127] = 1;
            slotsA[base + 128] = 1.0;
            slotsB[base + 127] = 1.0;
        }

        // Test 2: (2^256 - 1)^2 = 1 mod 2^256.
        {
            const uint32_t base = 2 * ks256gpu::kStride;
            for (uint32_t bit = 0; bit < ks256gpu::kWordBits; ++bit) {
                bitsA[base + bit] = 1;
                bitsB[base + bit] = 1;
                slotsA[base + bit] = 1.0;
                slotsB[base + bit] = 1.0;
            }
        }

        std::vector<uint8_t> expected;
        BuildExpectedProduct(bitsA, bitsB, &expected);
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
        std::cout << "  word 0 P = " << BitsToHex(expected, 0) << "\n";
        std::cout << "  word 1 A = " << BitsToHex(bitsA, 1) << "\n";
        std::cout << "  word 1 B = " << BitsToHex(bitsB, 1) << "\n";
        std::cout << "  word 1 P = " << BitsToHex(expected, 1) << "\n";
        std::cout << "  word 2 A = " << BitsToHex(bitsA, 2) << "\n";
        std::cout << "  word 2 B = " << BitsToHex(bitsB, 2) << "\n";
        std::cout << "  word 2 P = " << BitsToHex(expected, 2) << "\n";

        // --------------------------------------------------------------
        // 4. CPU Dadda + Kogge-Stone multiplication.
        // --------------------------------------------------------------
        ks256gpu::MultiplierTimingReport multiplierTiming;
begin = Clock::now();
        auto cProduct = ks256gpu::EvalMulDaddaKS(
            cc, cA, cB, directionSign, masks,
            &multiplierTiming, compressorBatchSize);
phase.multiplicationMs = Milliseconds(begin, Clock::now());

        std::cout << "\n涔樻硶鍚庡瘑鏂囩姸鎬?\n";
        std::cout << "  GetLevel()         : " << cProduct->GetLevel() << "\n";
        std::cout << "  GetNoiseScaleDeg() : "
                  << cProduct->GetNoiseScaleDeg() << "\n";
        std::cout << "  浼扮畻鍓╀綑灞傛暟      : "
                  << LevelsRemaining(multiplicativeDepth, cProduct) << "\n";

        // --------------------------------------------------------------
        // 5. Optional bootstrap. It is disabled for the first memory-bounded
        //    correctness run; decrypt the multiplication result directly.
        // --------------------------------------------------------------
        auto cResult = cProduct;
        if constexpr (enableBootstrap) {
begin = Clock::now();
            cResult = cc->EvalBootstrap(cProduct);
phase.bootstrapMs = Milliseconds(begin, Clock::now());

            std::cout << "\nBootstrap 鍚庡瘑鏂囩姸鎬?\n";
            std::cout << "  GetLevel()         : "
                      << cResult->GetLevel() << "\n";
            std::cout << "  GetNoiseScaleDeg() : "
                      << cResult->GetNoiseScaleDeg() << "\n";
            std::cout << "  浼扮畻鍓╀綑灞傛暟      : "
                      << LevelsRemaining(
                             multiplicativeDepth, cResult)
                      << "\n";
        }
        else {
            std::cout << "\nBootstrap 宸插叧闂細鐩存帴瑙ｅ瘑涔樻硶缁撴灉銆俓n";
        }

        // --------------------------------------------------------------
        // 6. Copy to CPU and decrypt.
        // --------------------------------------------------------------
        begin = Clock::now();
        Plaintext decrypted;
        const auto decryptResult =
            cc->Decrypt(keys.secretKey, cResult, &decrypted);
        if (!decryptResult.isValid) {
            throw std::runtime_error(
                "Decryption returned an invalid result");
        }
        decrypted->SetLength(numSlots);
        const std::vector<double> values =
            decrypted->GetRealPackedValue();
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
            static_cast<uint64_t>(ks256gpu::kMultiplierNumWords) *
            ks256gpu::kWordBits;

        for (uint32_t word = 0; word < ks256gpu::kMultiplierNumWords; ++word) {
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

        PrintMultiplierProfile(multiplierTiming);

        std::cout << "\n================ 绔埌绔椂闂村垎鏋?================\n";
        PrintPhase("CryptoContext 鐢熸垚涓?Enable", phase.contextMs);
        PrintPhase("KeyGen", phase.keyGenMs);
        PrintPhase("EvalMultKeyGen", phase.multKeyGenMs);
        PrintPhase("涔樻硶鍣ㄦ棆杞瘑閽ョ敓鎴?, phase.rotateKeyGenMs);
        if constexpr (enableBootstrap) {
            PrintPhase("EvalBootstrapSetup", phase.bootstrapSetupMs);
            PrintPhase("EvalBootstrapKeyGen", phase.bootstrapKeyGenMs);
        }
        PrintPhase("LoadContext 涓庤瘎浼板瘑閽ュ姞杞藉埌 CPU", phase.loadContextMs);
        PrintPhase("鎺╃爜缂栫爜骞跺姞杞藉埌 CPU", phase.maskBuildLoadMs);
        PrintPhase("杈撳叆鐢熸垚涓?CPU 鍙傝€冪粨鏋?, phase.inputBuildMs);
        PrintPhase("CKKS 缂栫爜", phase.encodingMs);
        PrintPhase("杈撳叆鍔犲瘑骞跺姞杞藉埌 CPU", phase.encryptionMs);
        PrintPhase("CPU Dadda + Kogge-Stone 涔樻硶", phase.multiplicationMs);
        if constexpr (enableBootstrap) {
            PrintPhase("CPU 涓€娆?EvalBootstrap", phase.bootstrapMs);
        }
        PrintPhase("澶嶅埗鍥?CPU 骞惰В瀵?, phase.decryptionMs);
        PrintPhase("姝ｇ‘鎬т笌璇樊鎵弿", phase.verificationMs);
        PrintPhase("绋嬪簭鎬诲閽熸椂闂?, phase.totalMs);

        const double mulSeconds = phase.multiplicationMs / 1000.0;
        const double bootSeconds = phase.bootstrapMs / 1000.0;
        const double onlineSeconds = mulSeconds + bootSeconds;

        std::cout << "\n================ 鍚炲悙涓庢憡閿€ ================\n";
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "  涔樻硶鎽婇攢鏃堕棿/256-bit word : "
                  << mulSeconds / ks256gpu::kMultiplierNumWords << " s\n";
        if constexpr (enableBootstrap) {
            std::cout << "  Boot 鎽婇攢鏃堕棿/word        : "
                      << bootSeconds / ks256gpu::kMultiplierNumWords << " s\n";
        }
        std::cout << "  鍦ㄧ嚎鎬绘憡閿€鏃堕棿/word       : "
                  << onlineSeconds / ks256gpu::kMultiplierNumWords << " s\n";

        if (mulSeconds > 0.0) {
            std::cout << "  浠呬箻娉曞悶鍚愰噺              : "
                      << ks256gpu::kMultiplierNumWords / mulSeconds
                      << " word/s\n";
        }
        if (onlineSeconds > 0.0) {
            std::cout << "  "
                      << (enableBootstrap
                              ? "涔樻硶+Boot 鍦ㄧ嚎鍚炲悙閲?     : "
                              : "涔樻硶鍦ㄧ嚎鍚炲悙閲?          : ")
                      << ks256gpu::kMultiplierNumWords / onlineSeconds
                      << " word/s\n";
        }

        return bitFailures == 0 ? 0 : 2;
    }
    catch (const std::exception& error) {
        std::cerr << "\n[FATAL] " << error.what() << "\n";
        return 1;
    }
}




