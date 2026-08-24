#include <openfhe.h>

#include <algorithm>
#include <chrono>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <vector>

namespace {

using namespace lbcrypto;
using Ct = Ciphertext<DCRTPoly>;
using Cc = CryptoContext<DCRTPoly>;
using Clock = std::chrono::steady_clock;

constexpr std::uint32_t kRingDimension = 131072;
constexpr std::uint32_t kSlots = 65536;
constexpr std::uint32_t kWordBits = 64;
constexpr std::uint32_t kGuardBits = 64;
constexpr std::uint32_t kStride = kWordBits + kGuardBits;
constexpr std::uint32_t kPackedCapacity = kSlots / kStride;  // 512
constexpr std::uint32_t kTransactions = 3;
constexpr std::uint32_t kLayers = 6;
constexpr std::uint32_t kDepth = 25;
constexpr std::uint32_t kScalingBits = 59;
constexpr std::uint32_t kFirstBits = 60;

class TeeBuffer : public std::streambuf {
public:
    TeeBuffer(std::streambuf* console, std::streambuf* log)
        : console_(console), log_(log) {}

protected:
    int overflow(int character) override {
        if (character == traits_type::eof()) {
            return traits_type::not_eof(character);
        }
        const char value = static_cast<char>(character);
        const bool a = console_->sputc(value) != traits_type::eof();
        const bool b = log_->sputc(value) != traits_type::eof();
        return a && b ? character : traits_type::eof();
    }

    std::streamsize xsputn(const char* data,
                           std::streamsize size) override {
        return std::min(console_->sputn(data, size),
                        log_->sputn(data, size));
    }

    int sync() override {
        return console_->pubsync() == 0 && log_->pubsync() == 0
                   ? 0
                   : -1;
    }

private:
    std::streambuf* console_;
    std::streambuf* log_;
};

class ScopedLogTee {
public:
    explicit ScopedLogTee(const std::string& path)
        : log_(path, std::ios::out | std::ios::trunc),
          original_(std::cout.rdbuf()),
          tee_(original_, log_.rdbuf()) {
        if (!log_) {
            throw std::runtime_error("cannot create log file: " + path);
        }
        std::cout.rdbuf(&tee_);
    }

    ~ScopedLogTee() {
        std::cout.flush();
        std::cout.rdbuf(original_);
    }

private:
    std::ofstream log_;
    std::streambuf* original_;
    TeeBuffer tee_;
};

struct Masks {
    Plaintext active;
    Plaintext carry;
    Plaintext result;
};

struct InputRow {
    std::uint64_t sender;
    std::uint64_t receiver;
    std::uint64_t amount;
};

double Milliseconds(Clock::time_point begin) {
    return std::chrono::duration<double, std::milli>(
               Clock::now() - begin)
        .count();
}

Plaintext MakePlaintext(const Cc& cc,
                        const std::vector<double>& values) {
    auto plaintext =
        cc->MakeCKKSPackedPlaintext(values, 1, 0, nullptr, kSlots);
    plaintext->SetLength(kSlots);
    return plaintext;
}

Masks BuildMasks(const Cc& cc) {
    std::vector<double> active(kSlots, 0.0);
    std::vector<double> carry(kSlots, 0.0);
    std::vector<double> result(kSlots, 0.0);
    for (std::uint32_t word = 0; word < kPackedCapacity; ++word) {
        const std::uint32_t base = word * kStride;
        result[base] = 1.0;
        for (std::uint32_t bit = 0; bit < kWordBits; ++bit) {
            active[base + bit] = 1.0;
            if (bit != 0) {
                carry[base + bit] = 1.0;
            }
        }
    }
    return {MakePlaintext(cc, active),
            MakePlaintext(cc, carry),
            MakePlaintext(cc, result)};
}

void SetWord(std::vector<double>& slots,
             std::uint32_t word,
             std::uint64_t value) {
    const std::uint32_t base = word * kStride;
    for (std::uint32_t bit = 0; bit < kWordBits; ++bit) {
        slots[base + bit] = static_cast<double>((value >> bit) & 1ULL);
    }
}

Ct Encrypt(const Cc& cc,
           const PublicKey<DCRTPoly>& publicKey,
           const std::vector<double>& values) {
    return cc->Encrypt(publicKey, MakePlaintext(cc, values));
}

Ct EvalAddKS(const Cc& cc, const Ct& a, const Ct& b,
             const Masks& masks) {
    auto g = cc->EvalMult(a, b);
    auto p = cc->EvalSub(cc->EvalAdd(a, b), cc->EvalAdd(g, g));
    auto G = g;
    auto P = p;
    for (std::uint32_t layer = 0; layer < kLayers; ++layer) {
        const int32_t rotation = -static_cast<int32_t>(1U << layer);
        const auto oldG = G;
        const auto oldP = P;
        const auto rotatedG = cc->EvalRotate(oldG, rotation);
        const auto rotatedP = cc->EvalRotate(oldP, rotation);
        G = cc->EvalAdd(oldG, cc->EvalMult(rotatedG, oldP));
        if (layer + 1 < kLayers) {
            P = cc->EvalMult(rotatedP, oldP);
        }
    }
    auto carry = cc->EvalMult(cc->EvalRotate(G, -1), masks.carry);
    G = cc->EvalMult(G, masks.active);
    auto sum = cc->EvalAdd(cc->EvalAdd(a, b), carry);
    return cc->EvalSub(sum, cc->EvalAdd(G, G));
}

Ct EvalSubKS(const Cc& cc, const Ct& a, const Ct& b,
             const Masks& masks) {
    auto ab = cc->EvalMult(a, b);
    auto xorAB =
        cc->EvalSub(cc->EvalAdd(a, b), cc->EvalAdd(ab, ab));
    auto G = cc->EvalSub(b, ab);
    auto P = cc->EvalSub(masks.active, xorAB);
    for (std::uint32_t layer = 0; layer < kLayers; ++layer) {
        const int32_t rotation = -static_cast<int32_t>(1U << layer);
        const auto oldG = G;
        const auto oldP = P;
        const auto rotatedG = cc->EvalRotate(oldG, rotation);
        const auto rotatedP = cc->EvalRotate(oldP, rotation);
        G = cc->EvalAdd(oldG, cc->EvalMult(rotatedG, oldP));
        if (layer + 1 < kLayers) {
            P = cc->EvalMult(rotatedP, oldP);
        }
    }
    auto borrow = cc->EvalMult(cc->EvalRotate(G, -1), masks.carry);
    G = cc->EvalMult(G, masks.active);
    auto difference = cc->EvalSub(cc->EvalSub(a, b), borrow);
    return cc->EvalAdd(difference, cc->EvalAdd(G, G));
}

// Returns one result bit per word: 1 exactly when a < b.
Ct EvalLessThan(const Cc& cc, const Ct& a, const Ct& b,
                const Masks& masks) {
    auto ab = cc->EvalMult(a, b);
    auto xorAB =
        cc->EvalSub(cc->EvalAdd(a, b), cc->EvalAdd(ab, ab));
    auto G = cc->EvalSub(b, ab);  // (~a)&b
    auto P = cc->EvalSub(masks.active, xorAB);
    for (std::uint32_t layer = 0; layer < kLayers; ++layer) {
        const int32_t step = static_cast<int32_t>(1U << layer);
        const auto oldG = G;
        const auto oldP = P;
        const auto highG = cc->EvalRotate(oldG, step);
        const auto highP = cc->EvalRotate(oldP, step);
        G = cc->EvalAdd(highG, cc->EvalMult(highP, oldG));
        if (layer + 1 < kLayers) {
            P = cc->EvalMult(highP, oldP);
        }
    }
    return cc->EvalMult(G, masks.result);
}

Ct BroadcastSign(const Cc& cc, Ct sign) {
    for (std::uint32_t layer = 0; layer < kLayers; ++layer) {
        const int32_t step = -static_cast<int32_t>(1U << layer);
        sign = cc->EvalAdd(sign, cc->EvalRotate(sign, step));
    }
    return sign;
}

// Implements sign*a + (1-sign)*b with two ciphertext multiplications.
Ct EvalSelect(const Cc& cc,
              const Ct& signWord,
              const Ct& a,
              const Ct& b,
              const Masks& masks) {
    auto signA = cc->EvalMult(signWord, a);
    auto oneMinusSign = cc->EvalSub(masks.active, signWord);
    auto notSignB = cc->EvalMult(oneMinusSign, b);
    return cc->EvalAdd(signA, notSignB);
}

std::uint64_t DecodeWord(const std::vector<std::complex<double>>& slots,
                         std::uint32_t word) {
    std::uint64_t value = 0;
    const std::uint32_t base = word * kStride;
    for (std::uint32_t bit = 0; bit < kWordBits; ++bit) {
        if (slots[base + bit].real() > 0.5) {
            value |= std::uint64_t{1} << bit;
        }
    }
    return value;
}

std::vector<std::complex<double>> DecryptSlots(
    const Cc& cc,
    const PrivateKey<DCRTPoly>& secretKey,
    const Ct& ciphertext) {
    Plaintext plaintext;
    cc->Decrypt(secretKey, ciphertext, &plaintext);
    plaintext->SetLength(kSlots);
    return plaintext->GetCKKSPackedValue();
}

}  // namespace

int main() {
    try {
        ScopedLogTee log("transaction_64bit_cpu.log");
        std::cout
            << "64-bit Boolean confidential transaction (CPU)\n"
            << "--------------------------------------------------------\n"
            << "inputs              = sender, receiver, amount\n"
            << "circuit             = FHE.lt + FHE.sub + 2xFHE.select + FHE.add\n"
            << "select              = sign*a + (1-sign)*b\n"
            << "layout capacity     = [64 data][64 guard] x 512\n"
            << "active transactions = " << kTransactions << '\n'
            << "ring dimension N    = " << kRingDimension << '\n'
            << "slots               = " << kSlots << '\n'
            << "multiplicative depth= " << kDepth << '\n'
            << "requested Q         = 60 + 25*59 = 1535 bits\n"
            << "log file            = transaction_64bit_cpu.log\n";

        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetSecurityLevel(HEStd_128_classic);
        parameters.SetSecretKeyDist(SPARSE_TERNARY);
        parameters.SetRingDim(kRingDimension);
        parameters.SetBatchSize(kSlots);
        parameters.SetMultiplicativeDepth(kDepth);
        parameters.SetScalingTechnique(FIXEDAUTO);
        parameters.SetScalingModSize(kScalingBits);
        parameters.SetFirstModSize(kFirstBits);
        parameters.SetKeySwitchTechnique(HYBRID);
        parameters.SetNumLargeDigits(3);

        auto cc = GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);

        const auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);
        std::vector<int32_t> rotations;
        for (std::uint32_t step = 1; step < kWordBits; step <<= 1U) {
            rotations.push_back(static_cast<int32_t>(step));
            rotations.push_back(-static_cast<int32_t>(step));
        }
        cc->EvalRotateKeyGen(keys.secretKey, rotations);
        const Masks masks = BuildMasks(cc);

        std::mt19937_64 rng(0x5452414E53414354ULL);
        std::uniform_int_distribution<std::uint64_t> balance(1000, 1000000);
        std::uniform_int_distribution<std::uint64_t> amount(1, 1500000);
        std::vector<InputRow> rows(kTransactions);
        std::vector<double> senderSlots(kSlots, 0.0);
        std::vector<double> receiverSlots(kSlots, 0.0);
        std::vector<double> amountSlots(kSlots, 0.0);
        std::vector<double> zeroSlots(kSlots, 0.0);
        for (std::uint32_t word = 0; word < kTransactions; ++word) {
            rows[word] = {balance(rng), balance(rng), amount(rng)};
            SetWord(senderSlots, word, rows[word].sender);
            SetWord(receiverSlots, word, rows[word].receiver);
            SetWord(amountSlots, word, rows[word].amount);
        }

        const auto cSender = Encrypt(cc, keys.publicKey, senderSlots);
        const auto cReceiver = Encrypt(cc, keys.publicKey, receiverSlots);
        const auto cAmount = Encrypt(cc, keys.publicKey, amountSlots);
        const auto cZero = Encrypt(cc, keys.publicKey, zeroSlots);

        const auto transactionBegin = Clock::now();

        auto phaseBegin = Clock::now();
        auto sign = EvalLessThan(cc, cSender, cAmount, masks);
        const double lessThanMs = Milliseconds(phaseBegin);

        phaseBegin = Clock::now();
        auto debited = EvalSubKS(cc, cSender, cAmount, masks);
        const double subtractMs = Milliseconds(phaseBegin);

        phaseBegin = Clock::now();
        auto signWord = BroadcastSign(cc, sign);
        const double signBroadcastMs = Milliseconds(phaseBegin);

        phaseBegin = Clock::now();
        auto newSender =
            EvalSelect(cc, signWord, cSender, debited, masks);
        const double selectSenderMs = Milliseconds(phaseBegin);

        phaseBegin = Clock::now();
        auto transferred =
            EvalSelect(cc, signWord, cZero, cAmount, masks);
        const double selectAmountMs = Milliseconds(phaseBegin);

        phaseBegin = Clock::now();
        auto newReceiver =
            EvalAddKS(cc, cReceiver, transferred, masks);
        const double addMs = Milliseconds(phaseBegin);

        const double transactionMs = Milliseconds(transactionBegin);
        const double amortizedMs =
            transactionMs / static_cast<double>(kTransactions);
        const double throughput = 1000.0 / amortizedMs;

        std::cout << std::fixed << std::setprecision(3)
                  << "\nbatch 1/1 (" << kTransactions
                  << " transactions) lt-to-final-add: "
                  << transactionMs << " ms\n"
                  << "\nTiming (ms, encryption excluded)\n"
                  << "  FHE.lt                  : " << lessThanMs << '\n'
                  << "  FHE.sub                 : " << subtractMs << '\n'
                  << "  sign broadcast          : " << signBroadcastMs << '\n'
                  << "  FHE.select sender       : " << selectSenderMs << '\n'
                  << "  FHE.select amount       : " << selectAmountMs << '\n'
                  << "  FHE.add                 : " << addMs << '\n'
                  << "  transaction total       : " << transactionMs << '\n'
                  << "  amortized per transaction: " << amortizedMs << '\n'
                  << "  throughput (ops/s)      : " << throughput << '\n'
                  << "\nOutput levels\n"
                  << "  sign=" << sign->GetLevel()
                  << ", debited=" << debited->GetLevel()
                  << ", newSender=" << newSender->GetLevel()
                  << ", transferred=" << transferred->GetLevel()
                  << ", newReceiver=" << newReceiver->GetLevel()
                  << '\n';

        const auto senderResult =
            DecryptSlots(cc, keys.secretKey, newSender);
        const auto receiverResult =
            DecryptSlots(cc, keys.secretKey, newReceiver);
        const auto signResult =
            DecryptSlots(cc, keys.secretKey, sign);
        bool pass = true;
        for (std::uint32_t word = 0; word < kTransactions; ++word) {
            const bool insufficient =
                rows[word].sender < rows[word].amount;
            const std::uint64_t expectedSender =
                insufficient
                    ? rows[word].sender
                    : rows[word].sender - rows[word].amount;
            const std::uint64_t expectedReceiver =
                insufficient
                    ? rows[word].receiver
                    : rows[word].receiver + rows[word].amount;
            pass = pass &&
                   (signResult[word * kStride].real() > 0.5) ==
                       insufficient &&
                   DecodeWord(senderResult, word) == expectedSender &&
                   DecodeWord(receiverResult, word) == expectedReceiver;
        }

        const auto& first = rows.front();
        const bool firstInsufficient = first.sender < first.amount;
        std::cout << "\nCorrectness\n"
                  << "  first input sender   : " << first.sender << '\n'
                  << "  first input receiver : " << first.receiver << '\n'
                  << "  first input amount   : " << first.amount << '\n'
                  << "  first sign (lt)      : "
                  << (firstInsufficient ? 1 : 0) << '\n'
                  << "  decrypted sign       : "
                  << (signResult[0].real() > 0.5 ? 1 : 0) << '\n'
                  << "  decrypted sender     : "
                  << DecodeWord(senderResult, 0) << '\n'
                  << "  decrypted receiver   : "
                  << DecodeWord(receiverResult, 0) << '\n'
                  << "  all 3 transactions   : "
                  << (pass ? "PASS" : "FAIL") << '\n';

        return pass ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
