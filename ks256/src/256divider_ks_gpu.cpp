#include "divider_ks_gpu.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace lbcrypto;
using boost::multiprecision::cpp_int;

namespace {

using Clock = std::chrono::steady_clock;

struct PhaseTimes {
    double contextMs = 0.0;
    double keyGenMs = 0.0;
    double multKeyGenMs = 0.0;
    double rotateKeyGenMs = 0.0;
    double bootstrapSetupMs = 0.0;
    double bootstrapKeyGenMs = 0.0;
    double loadContextMs = 0.0;
    double masksMs = 0.0;
    double inputMs = 0.0;
    double encodingMs = 0.0;
    double encryptionMs = 0.0;
    double divisionMs = 0.0;
    double decryptionMs = 0.0;
    double verificationMs = 0.0;
    double totalMs = 0.0;
};

inline double Milliseconds(const Clock::time_point& begin,
                           const Clock::time_point& end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

void PrintPhase(const std::string& name, double milliseconds) {
    std::cout << "  " << std::left << std::setw(42) << name
              << std::right << std::setw(14)
              << std::fixed << std::setprecision(3)
              << milliseconds << " ms"
              << "  (" << std::setprecision(6)
              << milliseconds / 1000.0 << " s)\n";
}

cpp_int RandomUint256(std::mt19937_64& rng) {
    cpp_int value = 0;
    for (uint32_t limb = 0; limb < 4; ++limb) {
        value |= cpp_int(rng()) << (64U * limb);
    }
    return value;
}

void EncodeWord(const cpp_int& value,
                uint32_t word,
                std::vector<double>* slots) {
    const uint32_t base = word * ks256gpu::kDividerStride;
    for (uint32_t bit = 0; bit < ks256gpu::kDividerWordBits; ++bit) {
        const unsigned digit =
            ((value >> bit) & 1).convert_to<unsigned>();
        (*slots)[base + bit] = static_cast<double>(digit);
    }
}

cpp_int DecodeWord(const std::vector<uint8_t>& bits, uint32_t word) {
    cpp_int value = 0;
    const uint32_t base = word * ks256gpu::kDividerStride;
    for (uint32_t bit = 0; bit < ks256gpu::kDividerWordBits; ++bit) {
        if (bits[base + bit] != 0U) {
            value |= cpp_int(1) << bit;
        }
    }
    return value;
}

std::string ToHex256(cpp_int value) {
    constexpr char digits[] = "0123456789abcdef";
    std::string text(64, '0');
    for (int i = 63; i >= 0; --i) {
        text[static_cast<size_t>(i)] =
            digits[(value & 0xF).convert_to<unsigned>()];
        value >>= 4;
    }
    return "0x" + text;
}

std::vector<int32_t> BuildRotationIndices() {
    std::vector<int32_t> indices;
    for (int32_t step = 1; step <= 512; step <<= 1) {
        indices.push_back(step);
        indices.push_back(-step);
    }
    return indices;
}

void PrintDividerReport(const ks256gpu::DividerTimingReport& report) {
    std::cout << "\n================ Radix-16 divider counters ================\n";
    std::cout << "  Radix-16 rounds                    : "
              << report.radixRounds << "\n";
    std::cout << "  Parallel candidate KS subtractions : "
              << report.candidateSubtractions << "\n";
    std::cout << "  Total KS additions                 : "
              << report.ksAdditions << "\n";
    std::cout << "  Total KS subtractions              : "
              << report.ksSubtractions << "\n";
    std::cout << "  Balanced candidate mux gates       : "
              << report.candidateMuxGates << "\n";
    std::cout << "  Ordinary CKKS bootstraps           : "
              << report.ordinaryBootstraps << "\n";
    std::cout << "  Boolean projection calls           : "
              << report.booleanProjections << "\n";
    std::cout << "  Per-trial state bootstraps          : "
              << report.roundStateBootstraps << "\n";
    std::cout << "  Multiple-preparation bootstraps     : "
              << report.multiplePreparationBootstraps << "\n";
    std::cout << "  Final-output bootstraps             : "
              << report.finalOutputBootstraps << "\n";
    PrintPhase("EvalDivRadix16KS total", report.totalMs);
}

} // namespace

int main() {
    try {
        std::cout << std::unitbuf;
        const auto programBegin = Clock::now();
        PhaseTimes phase;

        constexpr uint32_t ringDim = 1U << 17;
        constexpr uint32_t numSlots = ks256gpu::kDividerTotalSlots;
        constexpr uint32_t multiplicativeDepth = 41;
        constexpr uint32_t numLargeDigits = 3;
        constexpr uint64_t randomSeed = 0x524144495831364BULL;

        const std::vector<uint32_t> levelBudget = {4, 4};
        const std::vector<uint32_t> bsgsDim = {0, 0};

        ks256gpu::DividerOptions dividerOptions{
            .bootstrapInputMaxLevel = 39,
            .conciseProgress = true,
            .debugProgress = false,
        };
        if (const char* value = std::getenv("KS256_MAX_ROUNDS")) {
            dividerOptions.debugMaxRounds =
                static_cast<uint32_t>(std::stoul(value));
        }

        const bool diagnosticMode =
            std::getenv("KS256_DIAG_ALL") != nullptr;
        if (diagnosticMode) {
            dividerOptions.debugMaxRounds = 1U;
            dividerOptions.diagnosticStopAfterFirstRound = true;
            dividerOptions.conciseProgress = true;
        }

        constexpr double diagnosticWarningThreshold = 1.0e-3;
        constexpr double diagnosticFailureThreshold = 1.0e-2;

        uint32_t diagnosticPassed = 0;
        uint32_t diagnosticWarnings = 0;
        uint32_t diagnosticFailed = 0;
        std::vector<std::string> warnedDiagnosticNames;
        std::vector<std::string> failedDiagnosticNames;
        std::cout << "============================================================\n";
        std::cout << "OpenFHE CPU: " << ks256gpu::kDividerNumWords
                  << " x uint256 exact radix-16 Kogge-Stone divider"
                  << " (final optimized 768-BTS diagnostic build)\n";
        std::cout << "Boolean input/output; 4 Boolean bits per radix-16 digit\n";
        std::cout << "Single KS + projected borrow + refreshed state + serial G/P bootstrap\n";
        std::cout << "============================================================\n";

        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetSecretKeyDist(SPARSE_ENCAPSULATED);
        parameters.SetSecurityLevel(HEStd_128_classic);
        parameters.SetRingDim(ringDim);
        parameters.SetBatchSize(numSlots);
        parameters.SetMultiplicativeDepth(multiplicativeDepth);
        parameters.SetScalingTechnique(FLEXIBLEAUTO);
        parameters.SetScalingModSize(59);
        parameters.SetFirstModSize(60);
        parameters.SetKeySwitchTechnique(HYBRID);
        parameters.SetNumLargeDigits(numLargeDigits);
        const std::vector<int> devices{0, 1, 2};
auto begin = Clock::now();
        lbcrypto::CryptoContext<lbcrypto::DCRTPoly> cc =
            GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);
        cc->Enable(FHE);
        phase.contextMs = Milliseconds(begin, Clock::now());

        std::cout << "Parameters:\n";
        std::cout << "  Ring dimension                 : "
                  << cc->GetRingDimension() << "\n";
        std::cout << "  CKKS slots                     : " << numSlots << "\n";
        std::cout << "  Parallel uint256 divisions     : "
                  << ks256gpu::kDividerNumWords << "\n";
        std::cout << "  Radix                          : 16\n";
        std::cout << "  Boolean bits per digit         : 4\n";
        std::cout << "  Radix rounds                   : "
                  << ks256gpu::kRadixRounds << "\n";
        std::cout << "  Configured multiplicative depth: "
                  << multiplicativeDepth << "\n";
        std::cout << "  Secret-key/bootstrap mode       : SPARSE_ENCAPSULATED\n";
        std::cout << "  HYBRID large digits            : "
                  << numLargeDigits << "\n";
        std::cout << "  Logical OpenFHE devices       : ";
        for (size_t i = 0; i < devices.size(); ++i) {
            if (i != 0U) {
                std::cout << ", ";
            }
            std::cout << devices[i];
        }
        std::cout << "\n";
        std::cout << "  Trials per radix round          : 4 (8D,4D,2D,D)\n";
        std::cout << "  KS split                        : after stride 8\n";
        std::cout << "  Borrow handling                 : projection only\n";
        std::cout << "  State handling per trial        : Bootstrap + projection\n";
        std::cout << "  Quotient Bootstrap per round    : disabled\n";
        std::cout << "  G/P Bootstrap scheduling        : serial\n";
        std::cout << "  Diagnostic mode                 : "
                  << (diagnosticMode
                          ? "all first-round checkpoints"
                          : "off")
                  << "\n";


        // Build CPU-side test vectors and every plaintext before LoadContext.
        begin = Clock::now();
        std::vector<double> slotsA(numSlots, 0.0);
        std::vector<double> slotsD(numSlots, 0.0);
        std::vector<cpp_int> expectedQ(ks256gpu::kDividerNumWords);
        std::vector<cpp_int> expectedR(ks256gpu::kDividerNumWords);

        std::mt19937_64 rng(randomSeed);
        const cpp_int max256 = (cpp_int(1) << 256) - 1;
        for (uint32_t word = 0; word < ks256gpu::kDividerNumWords; ++word) {
            cpp_int a = RandomUint256(rng);
            cpp_int d = RandomUint256(rng);
            if (word == 0) { a = 100; d = 7; }
            else if (word == 1) { a = max256; d = 1; }
            else if (word == 2) { a = 0; d = 17; }
            else if (word == 3) { a = 12345; d = 67890; }
            else if (word == 4) { a = max256; d = max256; }
            else if (word == 5) { a = max256; d = 0; }
            else if (d == 0) { d = 1; }

            EncodeWord(a, word, &slotsA);
            EncodeWord(d, word, &slotsD);
            if (d == 0) {
                expectedQ[word] = 0;
                expectedR[word] = a;
            }
            else {
                expectedQ[word] = a / d;
                expectedR[word] = a % d;
            }
        }
        phase.inputMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        auto pA = cc->MakeCKKSPackedPlaintext(
            slotsA, 1, 0, nullptr, numSlots);
        auto pD = cc->MakeCKKSPackedPlaintext(
            slotsD, 1, 0, nullptr, numSlots);
        pA->SetLength(numSlots);
        pD->SetLength(numSlots);
        auto masks = ks256gpu::BuildDividerMasksOnCpu(cc);
        phase.encodingMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        auto keys = cc->KeyGen();
        phase.keyGenMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        cc->EvalMultKeyGen(keys.secretKey);
        phase.multKeyGenMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        cc->EvalRotateKeyGen(keys.secretKey, BuildRotationIndices());
        phase.rotateKeyGenMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        cc->EvalBootstrapSetup(levelBudget, bsgsDim, numSlots, 0);
        phase.bootstrapSetupMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        cc->EvalBootstrapKeyGen(keys.secretKey, numSlots);
        phase.bootstrapKeyGenMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        std::cout << "[init] loading CPU context and bootstrap keys...\n";
phase.loadContextMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        ks256gpu::LoadDividerMasksToGpu(cc, masks);
        phase.masksMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        auto cA = cc->Encrypt(keys.publicKey, pA);
        auto cD = cc->Encrypt(keys.publicKey, pD);
phase.encryptionMs = Milliseconds(begin, Clock::now());
        std::cout << "[init] inputs encrypted; starting exact-candidate radix-16 divider.\n";

        if (diagnosticMode) {
            dividerOptions.diagnosticSink =
                [&](const std::string& name,
                    const ks256gpu::Ciphertext<ks256gpu::DCRTPoly>& ciphertext,
                    ks256gpu::DividerDiagnosticLayout layout) {
                    std::cout << "\n[diag] " << name
                              << " | level=" << ciphertext->GetLevel()
                              << ", noiseScaleDeg="
                              << ciphertext->GetNoiseScaleDeg()
                              << "\n";

                    try {
                        auto diagnosticCiphertext = ciphertext;
                        Plaintext plaintext;
                        cc->Decrypt(
                            keys.secretKey,
                            diagnosticCiphertext,
                            &plaintext);
                        plaintext->SetLength(numSlots);

                        const auto& values =
                            plaintext->GetRealPackedValue();

                        double maxBooleanError = 0.0;
                        double maxAbsoluteValue = 0.0;
                        uint64_t checked = 0;
                        uint64_t nonFinite = 0;

                        auto inspectIndex = [&](uint32_t index) {
                            const double value = values[index];
                            ++checked;
                            if (!std::isfinite(value)) {
                                ++nonFinite;
                                return;
                            }

                            maxAbsoluteValue = std::max(
                                maxAbsoluteValue,
                                std::abs(value));
                            maxBooleanError = std::max(
                                maxBooleanError,
                                std::abs(value - std::round(value)));
                        };

                        for (uint32_t word = 0;
                             word < ks256gpu::kDividerNumWords;
                             ++word) {
                            const uint32_t base =
                                word * ks256gpu::kDividerStride;

                            if (layout ==
                                ks256gpu::DividerDiagnosticLayout::BitZero) {
                                inspectIndex(base);
                                continue;
                            }

                            const uint32_t width =
                                layout ==
                                        ks256gpu::DividerDiagnosticLayout::Active256
                                    ? ks256gpu::kDividerWordBits
                                    : ks256gpu::kDividerWideBits;

                            for (uint32_t bit = 0; bit < width; ++bit) {
                                inspectIndex(base + bit);
                            }
                        }

                        const bool hardFailure =
                            nonFinite != 0U ||
                            maxBooleanError >
                                diagnosticFailureThreshold;
                        const bool warning =
                            !hardFailure &&
                            maxBooleanError >
                                diagnosticWarningThreshold;

                        const char* status = "[diag][PASS]";
                        if (hardFailure) {
                            status = "[diag][FAIL]";
                            ++diagnosticFailed;
                            failedDiagnosticNames.push_back(name);
                        }
                        else if (warning) {
                            status = "[diag][WARN]";
                            ++diagnosticWarnings;
                            warnedDiagnosticNames.push_back(name);
                        }
                        else {
                            ++diagnosticPassed;
                        }

                        std::cout << status << " " << name
                                  << " | logPrecision="
                                  << plaintext->GetLogPrecision()
                                  << ", checked=" << checked
                                  << ", nonFinite=" << nonFinite
                                  << ", maxIntegerResidual="
                                  << std::scientific
                                  << maxBooleanError
                                  << ", maxAbs="
                                  << maxAbsoluteValue
                                  << std::fixed
                                  << "\n";
                    }
                    catch (const std::exception& error) {
                        ++diagnosticFailed;
                        failedDiagnosticNames.push_back(name);
                        std::cout << "[diag][FAIL] " << name
                                  << " | " << error.what()
                                  << "\n";
                    }
                };
        }

        ks256gpu::DividerTimingReport dividerReport;
        begin = Clock::now();
        auto result = ks256gpu::EvalDivRadix16KS(
            cc, cA, cD, masks, dividerOptions, &dividerReport);
phase.divisionMs = Milliseconds(begin, Clock::now());

        std::cout << "\n[pre-decrypt] quotient level="
                  << result.quotient->GetLevel()
                  << ", remainder level="
                  << result.remainder->GetLevel() << "\n";
        PrintDividerReport(dividerReport);

        if (diagnosticMode) {
            phase.totalMs = Milliseconds(programBegin, Clock::now());

            std::cout
                << "\n================ All-points diagnostic summary ================\n";
            std::cout << "  Checkpoints passed : "
                      << diagnosticPassed << "\n";
            std::cout << "  Checkpoint warnings: "
                      << diagnosticWarnings << "\n";
            std::cout << "  Checkpoints failed : "
                      << diagnosticFailed << "\n";
            std::cout << "  Warning threshold  : "
                      << std::scientific
                      << diagnosticWarningThreshold << "\n";
            std::cout << "  Failure threshold  : "
                      << diagnosticFailureThreshold
                      << std::fixed << "\n";

            if (!warnedDiagnosticNames.empty()) {
                std::cout << "  Warning checkpoints:\n";
                for (const auto& name : warnedDiagnosticNames) {
                    std::cout << "    - " << name << "\n";
                }
                std::cout
                    << "  The first warning in the log is the earliest "
                       "observed precision-risk point.\n";
            }

            if (!failedDiagnosticNames.empty()) {
                std::cout << "  Failed checkpoints:\n";
                for (const auto& name : failedDiagnosticNames) {
                    std::cout << "    - " << name << "\n";
                }
                std::cout
                    << "  The first failure in the log is the earliest "
                       "confirmed excessive-error point.\n";
            }
            else {
                std::cout
                    << "  No checkpoint exceeded the hard residual threshold "
                       "during the complete first radix-16 round.\n";
            }

            std::cout
                << "\n================ Diagnostic timing ================\n";
            PrintPhase("Encrypted diagnostic execution", phase.divisionMs);
            PrintPhase("Whole diagnostic program", phase.totalMs);

            return diagnosticFailed == 0U ? 0 : 3;
        }

        begin = Clock::now();
        std::cout << "[decrypt] remainder...\n";
        Plaintext pR;
        auto rStatus = cc->Decrypt(
            keys.secretKey, result.remainder, &pR);
        pR->SetLength(numSlots);
        std::cout << "[decrypt] remainder passed, log precision="
                << pR->GetLogPrecision() << "\n";

        std::cout << "[decrypt] quotient...\n";
        Plaintext pQ;
        auto qStatus = cc->Decrypt(
            keys.secretKey, result.quotient, &pQ);
        pQ->SetLength(numSlots);
        std::cout << "[decrypt] quotient passed, log precision="
                << pQ->GetLogPrecision() << "\n";
        phase.decryptionMs = Milliseconds(begin, Clock::now());

        begin = Clock::now();
        const auto& qValues = pQ->GetRealPackedValue();
        const auto& rValues = pR->GetRealPackedValue();
        std::vector<uint8_t> gotQ(numSlots, 0);
        std::vector<uint8_t> gotR(numSlots, 0);
        double maxQError = 0.0;
        double maxRError = 0.0;

        for (uint32_t word = 0; word < ks256gpu::kDividerNumWords; ++word) {
            const uint32_t base = word * ks256gpu::kDividerStride;
            for (uint32_t bit = 0; bit < ks256gpu::kDividerWordBits; ++bit) {
                const uint32_t index = base + bit;
                const double qr = std::round(qValues[index]);
                const double rr = std::round(rValues[index]);
                maxQError = std::max(maxQError, std::abs(qValues[index] - qr));
                maxRError = std::max(maxRError, std::abs(rValues[index] - rr));
                gotQ[index] = static_cast<uint8_t>(qValues[index] >= 0.5);
                gotR[index] = static_cast<uint8_t>(rValues[index] >= 0.5);
            }
        }

        uint32_t failed = 0;
        for (uint32_t word = 0; word < ks256gpu::kDividerNumWords; ++word) {
            const cpp_int q = DecodeWord(gotQ, word);
            const cpp_int r = DecodeWord(gotR, word);
            if (q != expectedQ[word] || r != expectedR[word]) {
                ++failed;
                if (failed <= 8U) {
                    std::cout << "[FAIL] word " << word << "\n"
                              << "  expected q: " << ToHex256(expectedQ[word]) << "\n"
                              << "  actual   q: " << ToHex256(q) << "\n"
                              << "  expected r: " << ToHex256(expectedR[word]) << "\n"
                              << "  actual   r: " << ToHex256(r) << "\n";
                }
            }
        }
        phase.verificationMs = Milliseconds(begin, Clock::now());
        phase.totalMs = Milliseconds(programBegin, Clock::now());

        std::cout << "\n================ Verification ================\n";
        std::cout << "  Failed packed words             : "
                  << failed << " / " << ks256gpu::kDividerNumWords << "\n";
        std::cout << "  Max quotient Boolean error      : "
                  << std::scientific << maxQError << "\n";
        std::cout << "  Max remainder Boolean error     : "
                  << std::scientific << maxRError << "\n";
        std::cout << std::fixed;

        std::cout << "\n================ End-to-end timing ================\n";
        PrintPhase("GenCryptoContext", phase.contextMs);
        PrintPhase("KeyGen", phase.keyGenMs);
        PrintPhase("EvalMultKeyGen", phase.multKeyGenMs);
        PrintPhase("EvalRotateKeyGen", phase.rotateKeyGenMs);
        PrintPhase("EvalBootstrapSetup", phase.bootstrapSetupMs);
        PrintPhase("EvalBootstrapKeyGen", phase.bootstrapKeyGenMs);
        PrintPhase("LoadContext", phase.loadContextMs);
        PrintPhase("Load divider masks", phase.masksMs);
        PrintPhase("Build CPU inputs", phase.inputMs);
        PrintPhase("CKKS encoding/masks", phase.encodingMs);
        PrintPhase("Encryption", phase.encryptionMs);
        PrintPhase("Encrypted radix-16 division", phase.divisionMs);
        PrintPhase("Decryption", phase.decryptionMs);
        PrintPhase("Exact verification", phase.verificationMs);
        PrintPhase("Whole program", phase.totalMs);

        if (failed != 0U) {
            return 2;
        }
        std::cout << "\nAll packed uint256 radix-16 divisions passed.\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "[FATAL] " << error.what() << "\n";
        return 1;
    }
}




