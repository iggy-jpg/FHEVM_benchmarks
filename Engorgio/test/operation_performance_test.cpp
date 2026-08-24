#include "ENGORGIO/comp.h"
#include "ENGORGIO/quant.h"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace lbcrypto;
using namespace openfhe;

namespace {

constexpr usint kRingDim = 131072;
constexpr usint kSlots = 65536;
constexpr usint kLevels = 32;
// OpenFHE adds the first modulus to the requested multiplicative-depth chain.
constexpr usint kMultiplicativeDepth = kLevels - 1;
constexpr usint kScalingModSize = 50;
constexpr usint kFirstModSize = 60;
constexpr std::uint32_t kLimbBits = 4;
constexpr std::uint32_t kPolynomialDegree = 119;
constexpr std::size_t kBasicOpRepetitions = 10;

using Clock = std::chrono::steady_clock;

double elapsedMilliseconds(const Clock::time_point& begin, const Clock::time_point& end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

void printTiming(const std::string& operation, std::uint32_t bits, double milliseconds) {
    const double millisecondsPerSlot = milliseconds / static_cast<double>(kSlots);
    const double nanosecondsPerSlot = millisecondsPerSlot * 1.0e6;

    std::cout << std::left << std::setw(6) << bits
              << std::setw(8) << operation
              << std::right << std::fixed << std::setprecision(6)
              << std::setw(16) << milliseconds << " ms"
              << std::setw(18) << millisecondsPerSlot << " ms/slot"
              << std::setw(16) << nanosecondsPerSlot << " ns/slot\n";
}

void printModulusChain(const CryptoContext<DCRTPoly>& cc) {
    const auto cryptoParams =
        std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    const auto& moduli = cryptoParams->GetElementParams()->GetParams();

    std::uint64_t sumOfPrimeBits = 0;
    long double log2Q = 0;
    std::cout << "Generated modulus chain (" << moduli.size() << " primes): ";
    for (std::size_t i = 0; i < moduli.size(); ++i) {
        const auto bits = moduli[i]->GetModulus().GetMSB();
        sumOfPrimeBits += bits;
        log2Q += std::log2(
            static_cast<long double>(moduli[i]->GetModulus().ConvertToDouble()));
        if (i != 0)
            std::cout << '/';
        std::cout << bits;
    }
    std::cout << " bits\n"
              << "Sum of modulus-prime bit sizes: " << sumOfPrimeBits << " bits\n"
              << "Q bit length: " << (static_cast<std::uint64_t>(std::floor(log2Q)) + 1)
              << " bits\n";
}

std::vector<Ciphertext<DCRTPoly>> encryptLimbs(
    const CryptoContext<DCRTPoly>& cc,
    const PublicKey<DCRTPoly>& publicKey,
    const std::vector<std::vector<double>>& limbs) {
    std::vector<Ciphertext<DCRTPoly>> ciphertexts;
    ciphertexts.reserve(limbs.size());
    for (const auto& limb : limbs)
        ciphertexts.push_back(cc->Encrypt(publicKey, cc->MakeCKKSPackedPlaintext(limb)));
    return ciphertexts;
}

void benchmarkBitWidth(const CryptoContext<DCRTPoly>& cc,
                       KeyPair<DCRTPoly>& keyPair,
                       std::uint32_t bitWidth) {
    std::mt19937_64 generator(0x454E474F5247494FULL + bitWidth);
    constexpr std::uint32_t limbMaximum = (std::uint32_t{1} << (kLimbBits - 1)) - 1;
    std::uniform_int_distribution<std::uint32_t> distribution(0, limbMaximum);

    std::vector<std::uint32_t> valuesA(kSlots);
    std::vector<std::uint32_t> valuesB(kSlots);
    std::vector<double> packedA(kSlots);
    std::vector<double> packedB(kSlots);
    for (std::size_t i = 0; i < kSlots; ++i) {
        for (std::uint32_t offset = 0; offset < bitWidth; offset += kLimbBits) {
            valuesA[i] |= distribution(generator) << offset;
            valuesB[i] |= distribution(generator) << offset;
        }
        packedA[i] = static_cast<double>(valuesA[i]);
        packedB[i] = static_cast<double>(valuesB[i]);
    }

    // Split through ENGORGIO and encrypt every 4-bit group independently.
    auto ciphertextA = encryptLimbs(cc, keyPair.publicKey,
                                    quantization(packedA, bitWidth, kLimbBits));
    auto ciphertextB = encryptLimbs(cc, keyPair.publicKey,
                                    quantization(packedB, bitWidth, kLimbBits));

    std::vector<Ciphertext<DCRTPoly>> results(ciphertextA.size());
    auto begin = Clock::now();
    for (std::size_t repetition = 0; repetition < kBasicOpRepetitions; ++repetition) {
        for (std::size_t limb = 0; limb < ciphertextA.size(); ++limb)
            results[limb] = cc->EvalAdd(ciphertextA[limb], ciphertextB[limb]);
    }
    auto end = Clock::now();
    printTiming("add", bitWidth,
                elapsedMilliseconds(begin, end) / static_cast<double>(kBasicOpRepetitions));

    begin = Clock::now();
    for (std::size_t repetition = 0; repetition < kBasicOpRepetitions; ++repetition) {
        for (std::size_t limb = 0; limb < ciphertextA.size(); ++limb)
            results[limb] = cc->EvalMult(ciphertextA[limb], ciphertextB[limb]);
    }
    end = Clock::now();
    printTiming("mult", bitWidth,
                elapsedMilliseconds(begin, end) / static_cast<double>(kBasicOpRepetitions));

    Ciphertext<DCRTPoly> greaterResult;
    constexpr double limbSignBound = (std::uint32_t{1} << (kLimbBits - 1)) - 1;

    begin = Clock::now();
    comp_greater_than_modular(ciphertextA, ciphertextB, limbSignBound, kPolynomialDegree,
                              greaterResult, keyPair.secretKey);
    end = Clock::now();
    printTiming("gt", bitWidth, elapsedMilliseconds(begin, end));
}

}  // namespace

int main() {
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecurityLevel(HEStd_128_classic);
    parameters.SetSecretKeyDist(SPARSE_TERNARY);
    parameters.SetRingDim(kRingDim);
    parameters.SetBatchSize(kSlots);
    parameters.SetMultiplicativeDepth(kMultiplicativeDepth);
    parameters.SetScalingModSize(kScalingModSize);
    parameters.SetFirstModSize(kFirstModSize);

    auto cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);

    auto keyPair = cc->KeyGen();
    if (!keyPair.good()) {
        std::cerr << "Key generation failed\n";
        return 1;
    }
    cc->EvalMultKeyGen(keyPair.secretKey);

    std::cout << "Operation performance test\n"
              << "N=" << cc->GetRingDimension()
              << ", slots=" << kSlots
              << ", levels=" << kLevels
              << ", multiplicative depth=" << kMultiplicativeDepth
              << ", first modulus=" << kFirstModSize
              << " bits, scaling modulus=" << kScalingModSize << " bits\n";
    printModulusChain(cc);
    std::cout << "\n"
              << std::left << std::setw(6) << "bits"
              << std::setw(8) << "op"
              << std::right << std::setw(19) << "latency"
              << std::setw(26) << "amortized"
              << std::setw(24) << "amortized\n";

    for (const auto bitWidth : std::array<std::uint32_t, 4>{4, 8, 16, 32})
        benchmarkBitWidth(cc, keyPair, bitWidth);

    return 0;
}
