#include <openfhe.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <streambuf>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace lbcrypto;
using Ct = Ciphertext<DCRTPoly>;
using Cc = CryptoContext<DCRTPoly>;
using Clock = std::chrono::steady_clock;

constexpr std::uint32_t kRingDimension = 131072;
constexpr std::uint32_t kSlots = 65536;
constexpr std::uint32_t kStandardMultiplicativeDepth = 25;
constexpr std::uint32_t kMultiplierMultiplicativeDepth = 40;
constexpr std::uint32_t kScalingModBits = 59;
constexpr std::uint32_t kFirstModBits = 60;
constexpr std::uint32_t kRepetitions = 1;

class TeeBuffer : public std::streambuf {
public:
    TeeBuffer(std::streambuf* console, std::streambuf* log)
        : console_(console), log_(log) {}

protected:
    int overflow(int character) override {
        if (character == traits_type::eof()) {
            return traits_type::not_eof(character);
        }
        const auto value = static_cast<char>(character);
        const bool consoleOk =
            console_->sputc(value) != traits_type::eof();
        const bool logOk = log_->sputc(value) != traits_type::eof();
        return consoleOk && logOk ? character : traits_type::eof();
    }

    std::streamsize xsputn(const char* data,
                           std::streamsize size) override {
        const auto consoleWritten = console_->sputn(data, size);
        const auto logWritten = log_->sputn(data, size);
        return std::min(consoleWritten, logWritten);
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

    ScopedLogTee(const ScopedLogTee&) = delete;
    ScopedLogTee& operator=(const ScopedLogTee&) = delete;

    std::ostream& Console() {
        return console_;
    }

private:
    std::ofstream log_;
    std::streambuf* original_;
    TeeBuffer tee_;
    std::ostream console_{original_};
};

struct Layout {
    std::uint32_t bits;
    std::uint32_t guard;
    std::uint32_t stride;
    std::uint32_t words;
    std::uint32_t layers;
};

struct Masks {
    Plaintext active;
    Plaintext carry;
    Plaintext result;
    Plaintext bitZero;
};

struct BenchResult {
    std::string op;
    std::uint32_t bits;
    std::uint32_t batch;
    double repeatedTotalMs;
    double averageMs;
};

std::vector<int> InputWord(const std::vector<double>& input,
                           const Layout& layout) {
    std::vector<int> bits(layout.bits);
    for (std::uint32_t bit = 0; bit < layout.bits; ++bit) {
        bits[bit] = input[bit] > 0.5 ? 1 : 0;
    }
    return bits;
}

std::string BitString(const std::vector<int>& bits) {
    std::string text = "0b";
    text.reserve(bits.size() + 2);
    for (auto it = bits.rbegin(); it != bits.rend(); ++it) {
        text.push_back(*it ? '1' : '0');
    }
    return text;
}

std::vector<int> ExpectedArithmetic(const std::vector<int>& a,
                                    const std::vector<int>& b,
                                    const std::string& op) {
    const std::size_t width = a.size();
    std::vector<int> result(width, 0);
    if (op == "add") {
        int carry = 0;
        for (std::size_t bit = 0; bit < width; ++bit) {
            const int value = a[bit] + b[bit] + carry;
            result[bit] = value & 1;
            carry = value >> 1;
        }
    } else if (op == "sub") {
        int borrow = 0;
        for (std::size_t bit = 0; bit < width; ++bit) {
            int value = a[bit] - b[bit] - borrow;
            borrow = value < 0 ? 1 : 0;
            if (value < 0) {
                value += 2;
            }
            result[bit] = value;
        }
    } else if (op == "mult") {
        for (std::size_t left = 0; left < width; ++left) {
            if (!a[left]) {
                continue;
            }
            for (std::size_t right = 0;
                 right + left < width;
                 ++right) {
                if (!b[right]) {
                    continue;
                }
                std::size_t bit = left + right;
                int carry = 1;
                while (carry && bit < width) {
                    const int value = result[bit] + carry;
                    result[bit] = value & 1;
                    carry = value >> 1;
                    ++bit;
                }
            }
        }
    } else {
        throw std::invalid_argument("unsupported arithmetic validation op");
    }
    return result;
}

int CompareWords(const std::vector<int>& a,
                 const std::vector<int>& b) {
    for (std::size_t bit = a.size(); bit-- > 0;) {
        if (a[bit] != b[bit]) {
            return a[bit] > b[bit] ? 1 : -1;
        }
    }
    return 0;
}

void VerifyResult(std::ostream& terminal,
                  const Cc& cc,
                  const PrivateKey<DCRTPoly>& secretKey,
                  const Ct& result,
                  const Layout& layout,
                  const std::vector<double>& inputA,
                  const std::vector<double>& inputB,
                  const std::string& op) {
    Plaintext decrypted;
    cc->Decrypt(secretKey, result, &decrypted);
    decrypted->SetLength(kSlots);
    const auto values = decrypted->GetCKKSPackedValue();
    const auto a = InputWord(inputA, layout);
    const auto b = InputWord(inputB, layout);

    terminal << "  correctness (first word, terminal only)\n";
    if (op == "gt" || op == "lt" || op == "eq" ||
        op == "ge" || op == "le") {
        const int comparison = CompareWords(a, b);
        int expected = 0;
        if (op == "gt") expected = comparison > 0;
        else if (op == "lt") expected = comparison < 0;
        else if (op == "eq") expected = comparison == 0;
        else if (op == "ge") expected = comparison >= 0;
        else expected = comparison <= 0;
        const int actual = values[0].real() > 0.5 ? 1 : 0;
        terminal << "    plaintext result = " << expected << '\n'
                 << "    decrypted result = " << actual << '\n'
                 << "    status           = "
                 << (expected == actual ? "PASS" : "FAIL")
                 << std::endl;
        return;
    }

    const auto expected = ExpectedArithmetic(a, b, op);
    std::vector<int> actual(layout.bits);
    bool pass = true;
    for (std::uint32_t bit = 0; bit < layout.bits; ++bit) {
        actual[bit] = values[bit].real() > 0.5 ? 1 : 0;
        pass = pass && actual[bit] == expected[bit];
    }
    terminal << "    plaintext result = " << BitString(expected) << '\n'
             << "    decrypted result = " << BitString(actual) << '\n'
             << "    status           = " << (pass ? "PASS" : "FAIL")
             << std::endl;
}

double ElapsedMs(Clock::time_point begin) {
    return std::chrono::duration<double, std::milli>(
               Clock::now() - begin)
        .count();
}

Layout MakeLayout(std::uint32_t bits) {
    if (bits < 4 || bits > 256 || (bits & (bits - 1U)) != 0) {
        throw std::invalid_argument(
            "bit width must be one of 4,8,16,32,64,128,256");
    }
    Layout layout{bits, bits, 2U * bits, kSlots / (2U * bits), 0};
    for (std::uint32_t value = bits; value > 1; value >>= 1U) {
        ++layout.layers;
    }
    return layout;
}

Plaintext MakePlaintext(const Cc& cc,
                        const std::vector<double>& values) {
    auto plaintext =
        cc->MakeCKKSPackedPlaintext(values, 1, 0, nullptr, kSlots);
    plaintext->SetLength(kSlots);
    return plaintext;
}

Masks BuildMasks(const Cc& cc, const Layout& layout) {
    std::vector<double> active(kSlots, 0.0);
    std::vector<double> carry(kSlots, 0.0);
    std::vector<double> result(kSlots, 0.0);
    std::vector<double> bitZero(kSlots, 0.0);

    for (std::uint32_t word = 0; word < layout.words; ++word) {
        const std::uint32_t base = word * layout.stride;
        result[base] = 1.0;
        bitZero[base] = 1.0;
        for (std::uint32_t bit = 0; bit < layout.bits; ++bit) {
            active[base + bit] = 1.0;
            if (bit != 0) {
                carry[base + bit] = 1.0;
            }
        }
    }
    return {MakePlaintext(cc, active), MakePlaintext(cc, carry),
            MakePlaintext(cc, result), MakePlaintext(cc, bitZero)};
}

std::vector<double> MakeInput(const Layout& layout,
                              std::uint64_t salt) {
    std::vector<double> values(kSlots, 0.0);
    for (std::uint32_t word = 0; word < layout.words; ++word) {
        std::uint64_t state =
            salt ^ (0x9e3779b97f4a7c15ULL * (word + 1ULL));
        const std::uint32_t base = word * layout.stride;
        for (std::uint32_t bit = 0; bit < layout.bits; ++bit) {
            // SplitMix-style deterministic bit stream; unlike integer-only
            // input, this exercises all bit positions for 128/256-bit tests.
            state += 0x9e3779b97f4a7c15ULL;
            std::uint64_t z = state;
            z = (z ^ (z >> 30U)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27U)) * 0x94d049bb133111ebULL;
            z ^= z >> 31U;
            values[base + bit] = static_cast<double>(z & 1ULL);
        }
    }
    return values;
}

Ct Encrypt(const Cc& cc,
           const PublicKey<DCRTPoly>& publicKey,
           const std::vector<double>& values) {
    return cc->Encrypt(publicKey, MakePlaintext(cc, values));
}

Ct EvalAddKS(const Cc& cc, const Ct& a, const Ct& b,
             const Layout& layout, const Masks& masks) {
    auto g = cc->EvalMult(a, b);
    auto p = cc->EvalSub(cc->EvalAdd(a, b), cc->EvalAdd(g, g));
    auto G = g;
    auto P = p;

    for (std::uint32_t layer = 0; layer < layout.layers; ++layer) {
        const int32_t rotation =
            -static_cast<int32_t>(1U << layer);
        const auto oldG = G;
        const auto oldP = P;
        auto rotatedG = cc->EvalRotate(oldG, rotation);
        auto rotatedP = cc->EvalRotate(oldP, rotation);
        G = cc->EvalAdd(oldG, cc->EvalMult(rotatedG, oldP));
        if (layer + 1 < layout.layers) {
            P = cc->EvalMult(rotatedP, oldP);
        }
    }

    auto carry = cc->EvalRotate(G, -1);
    carry = cc->EvalMult(carry, masks.carry);
    G = cc->EvalMult(G, masks.active);
    auto sum = cc->EvalAdd(a, b);
    sum = cc->EvalAdd(sum, carry);
    return cc->EvalSub(sum, cc->EvalAdd(G, G));
}

Ct EvalSubKS(const Cc& cc, const Ct& a, const Ct& b,
             const Layout& layout, const Masks& masks) {
    auto ab = cc->EvalMult(a, b);
    auto xorAB =
        cc->EvalSub(cc->EvalAdd(a, b), cc->EvalAdd(ab, ab));
    auto G = cc->EvalSub(b, ab);
    auto P = cc->EvalSub(masks.active, xorAB);

    for (std::uint32_t layer = 0; layer < layout.layers; ++layer) {
        const int32_t rotation =
            -static_cast<int32_t>(1U << layer);
        const auto oldG = G;
        const auto oldP = P;
        auto rotatedG = cc->EvalRotate(oldG, rotation);
        auto rotatedP = cc->EvalRotate(oldP, rotation);
        G = cc->EvalAdd(oldG, cc->EvalMult(rotatedG, oldP));
        if (layer + 1 < layout.layers) {
            P = cc->EvalMult(rotatedP, oldP);
        }
    }

    auto borrow = cc->EvalRotate(G, -1);
    borrow = cc->EvalMult(borrow, masks.carry);
    G = cc->EvalMult(G, masks.active);
    auto difference = cc->EvalSub(a, b);
    difference = cc->EvalSub(difference, borrow);
    return cc->EvalAdd(difference, cc->EvalAdd(G, G));
}

// Computes A>B when greater=true, and A>=B when greater=false.
Ct EvalOrderedCompare(const Cc& cc, const Ct& a, const Ct& b,
                      const Layout& layout, const Masks& masks,
                      bool greater) {
    auto ab = cc->EvalMult(a, b);
    auto xorAB =
        cc->EvalSub(cc->EvalAdd(a, b), cc->EvalAdd(ab, ab));
    // A>B uses A&~B. A>=B is the complement of A<B, whose generate
    // state is ~A&B.
    auto G = greater ? cc->EvalSub(a, ab) : cc->EvalSub(b, ab);
    auto P = cc->EvalSub(masks.active, xorAB);

    for (std::uint32_t layer = 0; layer < layout.layers; ++layer) {
        const int32_t step = static_cast<int32_t>(1U << layer);
        const auto oldG = G;
        const auto oldP = P;
        auto highG = cc->EvalRotate(oldG, step);
        auto highP = cc->EvalRotate(oldP, step);
        G = cc->EvalAdd(highG, cc->EvalMult(highP, oldG));
        if (layer + 1 < layout.layers) {
            P = cc->EvalMult(highP, oldP);
        }
    }

    auto result = cc->EvalMult(G, masks.result);
    if (!greater) {
        result = cc->EvalSub(masks.result, result);
    }
    return result;
}

Ct EvalEqualKS(const Cc& cc, const Ct& a, const Ct& b,
               const Layout& layout, const Masks& masks) {
    auto ab = cc->EvalMult(a, b);
    auto xorAB =
        cc->EvalSub(cc->EvalAdd(a, b), cc->EvalAdd(ab, ab));
    auto G = xorAB;
    auto P = cc->EvalSub(masks.active, xorAB);

    for (std::uint32_t layer = 0; layer < layout.layers; ++layer) {
        const int32_t rotation =
            -static_cast<int32_t>(1U << layer);
        const auto oldG = G;
        const auto oldP = P;
        auto lowerG = cc->EvalRotate(oldG, rotation);
        auto lowerP = cc->EvalRotate(oldP, rotation);
        G = cc->EvalAdd(oldG, cc->EvalMult(oldP, lowerG));
        if (layer + 1 < layout.layers) {
            P = cc->EvalMult(oldP, lowerP);
        }
    }

    auto notEqual = cc->EvalRotate(
        G, static_cast<int32_t>(layout.bits - 1U));
    notEqual = cc->EvalMult(notEqual, masks.result);
    return cc->EvalSub(masks.result, notEqual);
}

struct CompressedRows {
    Ct sum;
    Ct carry;
};

CompressedRows Compress3To2(const Cc& cc,
                            const Ct& x, const Ct& y, const Ct& z) {
    auto xy = cc->EvalMult(x, y);
    auto p = cc->EvalSub(cc->EvalAdd(x, y), cc->EvalAdd(xy, xy));
    auto pz = cc->EvalMult(p, z);
    auto sum =
        cc->EvalSub(cc->EvalAdd(p, z), cc->EvalAdd(pz, pz));
    auto carry = cc->EvalRotate(cc->EvalAdd(xy, pz), -1);
    return {std::move(sum), std::move(carry)};
}

std::vector<std::uint32_t> DaddaTargets(std::uint32_t bits) {
    std::vector<std::uint32_t> ascending{2};
    while (ascending.back() < bits) {
        ascending.push_back((3U * ascending.back()) / 2U);
    }
    std::vector<std::uint32_t> targets;
    for (auto it = ascending.rbegin(); it != ascending.rend(); ++it) {
        if (*it < bits) {
            targets.push_back(*it);
        }
    }
    if (targets.empty() || targets.back() != 2) {
        targets.push_back(2);
    }
    return targets;
}

struct DaddaStage {
    std::uint32_t inputRows = 0;
    std::uint32_t targetRows = 0;
    std::uint32_t fullAdders = 0;
    std::uint32_t received = 0;
    std::uint32_t emitted = 0;
    std::uint32_t completed = 0;
    std::array<Ct, 3> pending{};
    std::uint32_t pendingCount = 0;
};

class StreamingDadda {
public:
    StreamingDadda(const Cc& cc, std::uint32_t bits)
        : cc_(cc) {
        std::uint32_t rows = bits;
        for (std::uint32_t target : DaddaTargets(bits)) {
            DaddaStage stage;
            stage.inputRows = rows;
            stage.targetRows = target;
            stage.fullAdders = rows - target;
            stages_.push_back(std::move(stage));
            rows = target;
        }
    }

    void Push(Ct row) {
        PushAt(0, std::move(row));
    }

    std::array<Ct, 2> Finish() {
        for (const auto& stage : stages_) {
            if (stage.received != stage.inputRows ||
                stage.emitted != stage.targetRows ||
                stage.completed != stage.fullAdders ||
                stage.pendingCount != 0) {
                throw std::runtime_error("invalid Dadda pipeline state");
            }
        }
        if (finalRows_.size() != 2) {
            throw std::runtime_error("Dadda tree did not produce two rows");
        }
        return {std::move(finalRows_[0]), std::move(finalRows_[1])};
    }

private:
    void PushAt(std::size_t index, Ct row) {
        if (index == stages_.size()) {
            finalRows_.push_back(std::move(row));
            return;
        }
        auto& stage = stages_[index];
        ++stage.received;
        const std::uint32_t compressedInputs = 3U * stage.fullAdders;
        if (stage.received <= compressedInputs) {
            stage.pending[stage.pendingCount++] = std::move(row);
            if (stage.pendingCount != 3) {
                return;
            }
            auto compressed = Compress3To2(
                cc_, stage.pending[0], stage.pending[1], stage.pending[2]);
            stage.pending = {};
            stage.pendingCount = 0;
            ++stage.completed;
            stage.emitted += 2;
            PushAt(index + 1, std::move(compressed.sum));
            PushAt(index + 1, std::move(compressed.carry));
        } else {
            ++stage.emitted;
            PushAt(index + 1, std::move(row));
        }
    }

    Cc cc_;
    std::vector<DaddaStage> stages_;
    std::vector<Ct> finalRows_;
};

Ct EvalMultiplyDadda(const Cc& cc, const Ct& a, const Ct& b,
                     const Layout& layout, const Masks& masks) {
    StreamingDadda pipeline(cc, layout.bits);
    auto shiftedA = a;
    auto alignedB = b;

    for (std::uint32_t row = 0; row < layout.bits; ++row) {
        auto selected = cc->EvalMult(alignedB, masks.bitZero);
        auto broadcast = selected;
        for (std::uint32_t layer = 0; layer < layout.layers; ++layer) {
            const int32_t step =
                -static_cast<int32_t>(1U << layer);
            broadcast = cc->EvalAdd(
                broadcast, cc->EvalRotate(broadcast, step));
        }
        pipeline.Push(cc->EvalMult(shiftedA, broadcast));
        if (row + 1 < layout.bits) {
            shiftedA = cc->EvalRotate(shiftedA, -1);
            alignedB = cc->EvalRotate(alignedB, 1);
        }
    }

    auto rows = pipeline.Finish();
    auto product =
        EvalAddKS(cc, rows[0], rows[1], layout, masks);
    return cc->EvalMult(product, masks.active);
}

template <class Fn, class VerifyFn>
BenchResult RunOne(const std::string& name,
                   const Layout& layout,
                   Fn&& fn,
                   VerifyFn&& verify) {
    const auto begin = Clock::now();
    Ct result;
    for (std::uint32_t repetition = 0;
         repetition < kRepetitions;
         ++repetition) {
        result = fn();
        if (!result) {
            throw std::runtime_error(name + " returned a null ciphertext");
        }
    }
    const double repeatedTotalMs = ElapsedMs(begin);
    const double averageMs =
        repeatedTotalMs / static_cast<double>(kRepetitions);

    std::cout << "completed " << layout.bits << "-bit " << name
              << " x " << kRepetitions
              << ": repeated total = " << std::fixed
              << std::setprecision(3) << repeatedTotalMs
              << " ms, average = " << averageMs << " ms"
              << std::endl;

    verify(result);
    return {name, layout.bits, layout.words,
            repeatedTotalMs, averageMs};
}

void PrintTableRow(const BenchResult& result) {
    const double amortized =
        result.averageMs / static_cast<double>(result.batch);
    std::cout << std::left
              << std::setw(6) << result.bits
              << std::setw(8) << result.op
              << std::right
              << std::setw(8) << result.batch
              << std::setw(18) << std::fixed << std::setprecision(3)
              << result.averageMs
              << std::setw(20) << amortized
              << std::endl;
}

void PrintActualParameters(const std::string& label,
                           const Cc& cc,
                           std::uint32_t requestedDepth) {
    const auto params =
        std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(
            cc->GetCryptoParameters());
    const auto moduli = params->GetElementParams()->GetParams();
    std::uint64_t totalBits = 0;
    for (const auto& modulus : moduli) {
        totalBits += modulus->GetModulus().GetMSB();
    }
    std::cout << label << " context: requested depth = "
              << requestedDepth
              << ", actual Q towers = " << moduli.size()
              << ", available levels = "
              << (moduli.empty() ? 0 : moduli.size() - 1)
              << ", actual total Q bits = " << totalBits << '\n';
    const std::uint64_t requestedBits =
        kFirstModBits +
        static_cast<std::uint64_t>(requestedDepth) * kScalingModBits;
    if (totalBits != requestedBits) {
        std::cout
            << "WARNING: OpenFHE generated a modulus chain whose actual "
               "prime-bit sum differs from requested "
            << requestedBits << " bits.\n";
    }
    std::cout << std::endl;
}

}  // namespace

int main() {
    try {
        ScopedLogTee logTee("boolean_operator_benchmark_cpu.log");
        auto& terminal = logTee.Console();
        std::cout
            << "CKKS Boolean operator benchmark (CPU)\n"
            << "layout: [W data bits][W guard slots], total slots=65536\n"
            << "N=131072, slots=65536, first modulus=60, "
               "scaling modulus=59\n"
            << "add/sub/gt/eq/ge: depth=25, requested Q=1535 bits\n"
            << "mult: depth=40, requested Q=2420 bits\n\n";
        std::cout << "repetitions per operator = " << kRepetitions
                  << "\nlog file = boolean_operator_benchmark_cpu.log"
                  << "\n\n";

        std::vector<int32_t> rotations;
        for (std::uint32_t step = 1; step <= 128; step <<= 1U) {
            rotations.push_back(static_cast<int32_t>(step));
            rotations.push_back(-static_cast<int32_t>(step));
        }
        for (std::uint32_t bits : {4U, 8U, 16U, 32U, 64U, 128U, 256U}) {
            rotations.push_back(static_cast<int32_t>(bits - 1U));
        }
        std::sort(rotations.begin(), rotations.end());
        rotations.erase(
            std::unique(rotations.begin(), rotations.end()),
            rotations.end());

        auto makeContext = [&](std::uint32_t depth,
                               const std::string& label) {
            CCParams<CryptoContextCKKSRNS> parameters;
            parameters.SetSecurityLevel(HEStd_128_classic);
            parameters.SetSecretKeyDist(SPARSE_TERNARY);
            parameters.SetRingDim(kRingDimension);
            parameters.SetBatchSize(kSlots);
            parameters.SetMultiplicativeDepth(depth);
            parameters.SetScalingTechnique(FLEXIBLEAUTO);
            parameters.SetScalingModSize(kScalingModBits);
            parameters.SetFirstModSize(kFirstModBits);
            parameters.SetKeySwitchTechnique(HYBRID);
            parameters.SetNumLargeDigits(3);

            auto context = GenCryptoContext(parameters);
            context->Enable(PKE);
            context->Enable(KEYSWITCH);
            context->Enable(LEVELEDSHE);
            context->Enable(ADVANCEDSHE);
            PrintActualParameters(label, context, depth);

            auto keyPair = context->KeyGen();
            context->EvalMultKeyGen(keyPair.secretKey);
            context->EvalRotateKeyGen(keyPair.secretKey, rotations);
            return std::make_pair(context, keyPair);
        };

        auto standard = makeContext(
            kStandardMultiplicativeDepth, "standard operators");
        auto multiplier = makeContext(
            kMultiplierMultiplicativeDepth, "multiplier");
        const auto& standardCc = standard.first;
        const auto& standardKeys = standard.second;
        const auto& multiplierCc = multiplier.first;
        const auto& multiplierKeys = multiplier.second;

        std::vector<BenchResult> results;
        results.reserve(7U * 6U);

        for (std::uint32_t bits :
            {4U, 8U, 16U, 32U, 64U, 128U, 256U}) {
            const Layout layout = MakeLayout(bits);
            const Masks standardMasks = BuildMasks(standardCc, layout);
            const Masks multiplierMasks = BuildMasks(multiplierCc, layout);
            const auto inputA = MakeInput(layout, 0x123456789abcdef0ULL);
            const auto inputB = MakeInput(layout, 0xfedcba9876543210ULL);
            const auto standardA =
                Encrypt(standardCc, standardKeys.publicKey, inputA);
            const auto standardB =
                Encrypt(standardCc, standardKeys.publicKey, inputB);
            const auto multiplierA =
                Encrypt(multiplierCc, multiplierKeys.publicKey, inputA);
            const auto multiplierB =
                Encrypt(multiplierCc, multiplierKeys.publicKey, inputB);

            const auto validateStandard = [&](const std::string& op) {
                return [&, op](const Ct& result) {
                    VerifyResult(
                        terminal, standardCc, standardKeys.secretKey,
                        result, layout, inputA, inputB, op);
                };
            };
            const auto validateMultiplier = [&](const Ct& result) {
                    VerifyResult(
                        terminal, multiplierCc,
                        multiplierKeys.secretKey, result,
                        layout, inputA, inputB, "mult");
                };

            results.push_back(RunOne(
                "add", layout,
                [&] {
                    return EvalAddKS(
                        standardCc, standardA, standardB,
                        layout, standardMasks);
                },
                validateStandard("add")));
            results.push_back(RunOne(
                "sub", layout,
                [&] {
                    return EvalSubKS(
                        standardCc, standardA, standardB,
                        layout, standardMasks);
                },
                validateStandard("sub")));
            results.push_back(RunOne(
                "mult", layout,
                [&] {
                    return EvalMultiplyDadda(
                        multiplierCc, multiplierA, multiplierB,
                        layout, multiplierMasks);
                },
                validateMultiplier));
            results.push_back(RunOne(
                "gt", layout,
                [&] {
                    return EvalOrderedCompare(
                        standardCc, standardA, standardB,
                        layout, standardMasks, true);
                },
                validateStandard("gt")));
            results.push_back(RunOne(
                "lt", layout,
                [&] {
                    return EvalOrderedCompare(
                        standardCc, standardB, standardA,
                        layout, standardMasks, true);
                },
                validateStandard("lt")));
            results.push_back(RunOne(
                "eq", layout,
                [&] {
                    return EvalEqualKS(
                        standardCc, standardA, standardB,
                        layout, standardMasks);
                },
                validateStandard("eq")));
            results.push_back(RunOne(
                "ge", layout,
                [&] {
                    return EvalOrderedCompare(
                        standardCc, standardA, standardB,
                        layout, standardMasks, false);
                },
                validateStandard("ge")));
            results.push_back(RunOne(
                "le", layout,
                [&] {
                    return EvalOrderedCompare(
                        standardCc, standardB, standardA,
                        layout, standardMasks, false);
                },
                validateStandard("le")));
        }

        std::cout << "\nFinal averages\n";
        std::cout << std::left
                  << std::setw(6) << "bits"
                  << std::setw(8) << "op"
                  << std::right
                  << std::setw(8) << "batch"
                  << std::setw(18) << "latency(ms)"
                  << std::setw(20) << "amortized(ms)"
                  << std::endl;
        std::cout << std::string(60, '-') << std::endl;
        for (const auto& result : results) {
            PrintTableRow(result);
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
