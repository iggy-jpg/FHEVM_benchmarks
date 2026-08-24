#include "compare/compsigneval.hpp"
#include "openfhe.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace lbcrypto;
using Clock = std::chrono::steady_clock;

constexpr std::uint32_t kRingDimension = 65536;
constexpr std::uint32_t kSlots = 32768;

template <class Function>
double Measure(Function&& function)
{
    const auto begin = Clock::now();
    function();
    return std::chrono::duration<double, std::milli>(Clock::now() - begin).count();
}

double DecryptFirst(const CryptoContext<DCRTPoly>& cc,
                    const PrivateKey<DCRTPoly>& sk,
                    const Ciphertext<DCRTPoly>& value)
{
    Plaintext plaintext;
    cc->Decrypt(sk, value, &plaintext);
    plaintext->SetLength(1);
    return plaintext->GetRealPackedValue()[0];
}

int main()
{
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(20);
    parameters.SetScalingModSize(59);
    parameters.SetFirstModSize(60);
    parameters.SetSecurityLevel(HEStd_128_classic);
    parameters.SetRingDim(kRingDimension);
    parameters.SetScalingTechnique(FLEXIBLEAUTO);
    parameters.SetBatchSize(kSlots);
    auto cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    const auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    int failures = 0;
    auto encrypted = [&](double value) {
        return cc->Encrypt(
            keys.publicKey,
            cc->MakeCKKSPackedPlaintext(std::vector<double>(kSlots, value)));
    };
    auto report = [&](const char* op, int width, double ms, double actual, double expected) {
        const bool correct = std::fabs(actual - expected) < 0.2;
        failures += !correct;
        std::cout << std::fixed << std::setprecision(3)
                  << "PERF op=" << op << " W=" << width << " ms=" << ms
                  << " actual=" << actual << " expected=" << expected
                  << " correct=" << (correct ? "yes" : "no") << '\n';
    };

    for (int width : {4, 8, 16, 32}) {
        const double a = std::ldexp(1.0, width - 2) + 1.0;
        const double b = std::ldexp(1.0, width - 3) + 1.0;
        const auto ca = encrypted(a), cb = encrypted(b);
        Ciphertext<DCRTPoly> result;
        const double ms = Measure([&] { result = cc->EvalAdd(ca, cb); });
        report("add", width, ms, DecryptFirst(cc, keys.secretKey, result), a + b);
    }
    for (int width : {4, 8, 16}) {
        const double a = std::ldexp(1.0, width / 2 - 1) + 1.0;
        const double b = std::ldexp(1.0, width / 2 - 2) + 1.0;
        const auto ca = encrypted(a), cb = encrypted(b);
        Ciphertext<DCRTPoly> result;
        const double ms = Measure([&] { result = cc->EvalMult(ca, cb); });
        report("mult", width, ms, DecryptFirst(cc, keys.secretKey, result), a * b);
    }
    for (int width : {4, 8}) {
        const double equalValue = std::ldexp(1.0, width - 2) + 1.0;
        const auto equalA = encrypted(equalValue), equalB = encrypted(equalValue);
        Ciphertext<DCRTPoly> result;
        double ms = Measure([&] { result = compare_minimax::Equal(equalA, equalB); });
        report("eq", width, ms, DecryptFirst(cc, keys.secretKey, result), 1.0);

        const double a = std::ldexp(1.0, width - 2) + 1.0;
        const double b = std::ldexp(1.0, width - 3);
        const auto ca = encrypted(a), cb = encrypted(b);
        ms = Measure([&] { result = compare_minimax::GreaterThan(ca, cb); });
        report("gt", width, ms, DecryptFirst(cc, keys.secretKey, result), 1.0);
    }
    return failures == 0 ? 0 : 2;
}
