#include "tpch_boolean_ops.h"

#include <openfhe.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <vector>

namespace {

using namespace lbcrypto;
using Clock = std::chrono::steady_clock;
using CkksCiphertext = Ciphertext<DCRTPoly>;

constexpr std::uint32_t kRingDimension = 131072;
constexpr std::uint32_t kSlots = 65536;
constexpr std::uint32_t kMultiplicativeDepth = 25;
constexpr std::uint32_t kScalingModBits = 59;
constexpr std::uint32_t kFirstModBits = 60;
constexpr std::uint32_t kExpectedQBits = 1535;

constexpr std::uint16_t kShipDateBegin = 10592;
constexpr std::uint16_t kShipDateEnd = 10957;
constexpr std::uint16_t kDiscountLow = 5;
constexpr std::uint16_t kDiscountHigh = 7;
constexpr std::uint16_t kQuantityLimit = 24;

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

private:
    std::ofstream log_;
    std::streambuf* original_;
    TeeBuffer tee_;
};

struct LineItem {
    std::uint16_t extendedPrice;
    std::uint16_t discountPercent;
    std::uint16_t quantity;
    std::uint16_t shipDate;
};

struct Timings {
    double contextMs = 0.0;
    double keyGenMs = 0.0;
    double evalKeyGenMs = 0.0;
    double inputBuildMs = 0.0;
    double encryptionMs = 0.0;
    double comparisonMs = 0.0;
    double predicateAndMs = 0.0;
    double multiplicationMs = 0.0;
    double selectionMs = 0.0;
    double aggregationMs = 0.0;
    double decryptionMs = 0.0;
};

double Milliseconds(Clock::time_point begin) {
    return std::chrono::duration<double, std::milli>(
               Clock::now() - begin)
        .count();
}

std::size_t ParseRecordCount(int argc, char** argv) {
    if (argc == 1) {
        return kSlots;
    }
    if (argc != 2) {
        throw std::invalid_argument(
            "usage: tpch_q6_16bit_cpu [records, 1..65536]");
    }
    const unsigned long parsed = std::stoul(argv[1]);
    if (parsed == 0 || parsed > kSlots) {
        throw std::invalid_argument("records must be in [1, 65536]");
    }
    return static_cast<std::size_t>(parsed);
}

std::vector<LineItem> MakeInput(std::size_t records) {
    std::mt19937 engine(0x51364231U);
    std::uniform_int_distribution<std::uint16_t> price(100, 7000);
    std::uniform_int_distribution<std::uint16_t> discount(3, 9);
    std::uniform_int_distribution<std::uint16_t> quantity(1, 50);
    std::uniform_int_distribution<std::uint16_t> shipDate(10000, 11500);

    std::vector<LineItem> rows(records);
    for (auto& row : rows) {
        row = {price(engine), discount(engine), quantity(engine),
               shipDate(engine)};
    }
    return rows;
}

bool MatchesQ6(const LineItem& row) {
    return row.shipDate >= kShipDateBegin &&
           row.shipDate < kShipDateEnd &&
           row.discountPercent >= kDiscountLow &&
           row.discountPercent <= kDiscountHigh &&
           row.quantity < kQuantityLimit;
}

std::uint64_t PlainQ6(const std::vector<LineItem>& rows) {
    std::uint64_t revenue = 0;
    for (const auto& row : rows) {
        if (MatchesQ6(row)) {
            revenue += static_cast<std::uint64_t>(row.extendedPrice) *
                       row.discountPercent;
        }
    }
    return revenue;
}

void SetWord(std::vector<double>& slots,
             std::uint32_t word,
             std::uint16_t value) {
    const std::uint32_t base = word * tpch16::kStride;
    for (std::uint32_t bit = 0; bit < tpch16::kWordBits; ++bit) {
        slots[base + bit] = static_cast<double>((value >> bit) & 1U);
    }
}

std::vector<double> ConstantWords(std::uint16_t value) {
    std::vector<double> slots(kSlots, 0.0);
    for (std::uint32_t word = 0;
         word < tpch16::kRecordsPerBatch;
         ++word) {
        SetWord(slots, word, value);
    }
    return slots;
}

std::vector<double> ActiveResultMask(std::size_t activeRecords) {
    std::vector<double> mask(kSlots, 0.0);
    for (std::size_t word = 0; word < activeRecords; ++word) {
        mask[word * tpch16::kStride] = 1.0;
    }
    return mask;
}

std::vector<double> ActiveWordMask(std::size_t activeRecords) {
    std::vector<double> mask(kSlots, 0.0);
    for (std::size_t word = 0; word < activeRecords; ++word) {
        const std::size_t base = word * tpch16::kStride;
        for (std::size_t bit = 0; bit < tpch16::kWordBits; ++bit) {
            mask[base + bit] = 1.0;
        }
    }
    return mask;
}

std::vector<double> BitWeights() {
    std::vector<double> weights(kSlots, 0.0);
    for (std::uint32_t word = 0;
         word < tpch16::kRecordsPerBatch;
         ++word) {
        const std::uint32_t base = word * tpch16::kStride;
        for (std::uint32_t bit = 0; bit < tpch16::kWordBits; ++bit) {
            weights[base + bit] =
                static_cast<double>(std::uint32_t{1} << bit);
        }
    }
    return weights;
}

std::vector<double> ResultBitMask() {
    return ActiveResultMask(tpch16::kRecordsPerBatch);
}

CkksCiphertext EncryptSlots(
    const CryptoContext<DCRTPoly>& cc,
    const PublicKey<DCRTPoly>& publicKey,
    const std::vector<double>& values) {
    auto plaintext =
        cc->MakeCKKSPackedPlaintext(values, 1, 0, nullptr, kSlots);
    plaintext->SetLength(kSlots);
    return cc->Encrypt(publicKey, plaintext);
}

CkksCiphertext BroadcastPredicate(
    const CryptoContext<DCRTPoly>& cc,
    CkksCiphertext predicate) {
    for (std::uint32_t step = 1;
         step < tpch16::kWordBits;
         step <<= 1U) {
        predicate =
            cc->EvalAdd(predicate,
                        cc->EvalRotate(predicate,
                                       -static_cast<int32_t>(step)));
    }
    return predicate;
}

CkksCiphertext SumBatchRevenue(
    const CryptoContext<DCRTPoly>& cc,
    const CkksCiphertext& selectedRevenue,
    const Plaintext& bitWeights,
    const Plaintext& resultBitMask) {
    auto weighted = cc->EvalMult(selectedRevenue, bitWeights);

    // Collapse the 16 bit positions of every word into its bit-0 slot.
    auto wordValues = weighted;
    for (std::uint32_t bit = 1; bit < tpch16::kWordBits; ++bit) {
        wordValues =
            cc->EvalAdd(wordValues,
                        cc->EvalRotate(weighted,
                                       static_cast<int32_t>(bit)));
    }
    wordValues = cc->EvalMult(wordValues, resultBitMask);

    // Sum all 2,048 words. Only rotations/additions are used, so this part
    // consumes no additional ciphertext-multiplication depth.
    for (std::uint32_t words = 1;
         words < tpch16::kRecordsPerBatch;
         words <<= 1U) {
        const int32_t rotation =
            static_cast<int32_t>(words * tpch16::kStride);
        wordValues =
            cc->EvalAdd(wordValues, cc->EvalRotate(wordValues, rotation));
    }
    // Slot 0 already contains the sum of all packed words. The program only
    // decrypts slot 0, so a final first-slot plaintext multiplication would
    // be redundant and would consume one additional CKKS level.
    return wordValues;
}

void PrintParameters(const CryptoContext<DCRTPoly>& cc) {
    const auto params =
        std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(
            cc->GetCryptoParameters());
    const auto moduli = params->GetElementParams()->GetParams();
    std::uint64_t totalBits = 0;
    for (const auto& modulus : moduli) {
        totalBits += modulus->GetModulus().GetMSB();
    }
    std::cout << "actual Q towers     = " << moduli.size() << '\n'
              << "available levels    = "
              << (moduli.empty() ? 0 : moduli.size() - 1) << '\n'
              << "actual total Q bits = " << totalBits << '\n';
    if (totalBits != kExpectedQBits) {
        std::cout << "WARNING: expected Q=1535 bits, installed OpenFHE "
                     "generated Q="
                  << totalBits << " bits\n";
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        ScopedLogTee logTee("tpch_q6_16bit_cpu.log");
        const std::size_t records = ParseRecordCount(argc, argv);
        const std::size_t batches =
            (records + tpch16::kRecordsPerBatch - 1) /
            tpch16::kRecordsPerBatch;
        Timings timing;

        std::cout
            << "TPC-H Q6 - repository Boolean operators, 16-bit CPU\n"
            << "--------------------------------------------------------\n"
            << "records             = " << records << '\n'
            << "records per batch   = " << tpch16::kRecordsPerBatch << '\n'
            << "ciphertext batches  = " << batches << '\n'
            << "layout              = [16 data][16 guard] x 2048\n"
            << "ring dimension N    = " << kRingDimension << '\n'
            << "slots               = " << kSlots << '\n'
            << "multiplicative depth= " << kMultiplicativeDepth << '\n'
            << "requested Q         = 60 + 25*59 = 1535 bits\n"
            << "log file            = tpch_q6_16bit_cpu.log\n";

        auto begin = Clock::now();
        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetSecurityLevel(HEStd_128_classic);
        parameters.SetSecretKeyDist(SPARSE_TERNARY);
        parameters.SetRingDim(kRingDimension);
        parameters.SetBatchSize(kSlots);
        parameters.SetMultiplicativeDepth(kMultiplicativeDepth);
        parameters.SetScalingTechnique(FIXEDAUTO);
        parameters.SetScalingModSize(kScalingModBits);
        parameters.SetFirstModSize(kFirstModBits);
        parameters.SetKeySwitchTechnique(HYBRID);
        parameters.SetNumLargeDigits(3);

        auto cc = GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);
        timing.contextMs = Milliseconds(begin);

        begin = Clock::now();
        const auto keys = cc->KeyGen();
        timing.keyGenMs = Milliseconds(begin);

        begin = Clock::now();
        cc->EvalMultKeyGen(keys.secretKey);
        std::vector<int32_t> rotations;
        for (int32_t rotation = -15; rotation <= 15; ++rotation) {
            if (rotation != 0) {
                rotations.push_back(rotation);
            }
        }
        for (std::uint32_t words = 1;
             words < tpch16::kRecordsPerBatch;
             words <<= 1U) {
            rotations.push_back(
                static_cast<int32_t>(words * tpch16::kStride));
        }
        std::sort(rotations.begin(), rotations.end());
        rotations.erase(
            std::unique(rotations.begin(), rotations.end()),
            rotations.end());
        cc->EvalRotateKeyGen(keys.secretKey, rotations);
        timing.evalKeyGenMs = Milliseconds(begin);

        PrintParameters(cc);

        begin = Clock::now();
        const auto rows = MakeInput(records);
        const std::uint64_t expectedRevenue = PlainQ6(rows);
        timing.inputBuildMs = Milliseconds(begin);

        begin = Clock::now();
        const auto cShipBegin = EncryptSlots(
            cc, keys.publicKey, ConstantWords(kShipDateBegin - 1));
        const auto cShipEnd = EncryptSlots(
            cc, keys.publicKey, ConstantWords(kShipDateEnd));
        const auto cDiscountLow = EncryptSlots(
            cc, keys.publicKey, ConstantWords(kDiscountLow - 1));
        const auto cDiscountHigh = EncryptSlots(
            cc, keys.publicKey, ConstantWords(kDiscountHigh + 1));
        const auto cQuantityLimit = EncryptSlots(
            cc, keys.publicKey, ConstantWords(kQuantityLimit));
        const auto ptBitWeights =
            cc->MakeCKKSPackedPlaintext(BitWeights(), 1, 0, nullptr, kSlots);
        const auto ptResultBitMask =
            cc->MakeCKKSPackedPlaintext(ResultBitMask(), 1, 0, nullptr, kSlots);
        timing.encryptionMs += Milliseconds(begin);

        CkksCiphertext encryptedGrandTotal;
        std::vector<double> batchEvaluationMs;
        batchEvaluationMs.reserve(batches);

        for (std::size_t batch = 0; batch < batches; ++batch) {
            const std::size_t rowBegin =
                batch * tpch16::kRecordsPerBatch;
            const std::size_t activeRecords =
                std::min<std::size_t>(
                    tpch16::kRecordsPerBatch, records - rowBegin);

            begin = Clock::now();
            std::vector<double> price(kSlots, 0.0);
            std::vector<double> discount(kSlots, 0.0);
            std::vector<double> quantity(kSlots, 0.0);
            std::vector<double> shipDate(kSlots, 0.0);
            for (std::size_t word = 0; word < activeRecords; ++word) {
                const auto& row = rows[rowBegin + word];
                SetWord(price, word, row.extendedPrice);
                SetWord(discount, word, row.discountPercent);
                SetWord(quantity, word, row.quantity);
                SetWord(shipDate, word, row.shipDate);
            }
            const auto cPrice = EncryptSlots(cc, keys.publicKey, price);
            const auto cDiscount = EncryptSlots(cc, keys.publicKey, discount);
            const auto cQuantity = EncryptSlots(cc, keys.publicKey, quantity);
            const auto cShipDate = EncryptSlots(cc, keys.publicKey, shipDate);
            const auto ptActiveResults = cc->MakeCKKSPackedPlaintext(
                ActiveResultMask(activeRecords), 1, 0, nullptr, kSlots);
            const auto ptActiveWords = cc->MakeCKKSPackedPlaintext(
                ActiveWordMask(activeRecords), 1, 0, nullptr, kSlots);
            timing.encryptionMs += Milliseconds(begin);

            // Per-batch query timing deliberately starts after input
            // construction/encryption and ends after aggregation.
            const auto batchEvaluationBegin = Clock::now();
            double elapsed = 0.0;
            auto pShipBegin = tpch16::EvalGreaterThan(
                cc, cShipDate, cShipBegin, &elapsed);
            timing.comparisonMs += elapsed;
            auto pShipEnd = tpch16::EvalGreaterThan(
                cc, cShipEnd, cShipDate, &elapsed);
            timing.comparisonMs += elapsed;
            auto pDiscountLow = tpch16::EvalGreaterThan(
                cc, cDiscount, cDiscountLow, &elapsed);
            timing.comparisonMs += elapsed;
            auto pDiscountHigh = tpch16::EvalGreaterThan(
                cc, cDiscountHigh, cDiscount, &elapsed);
            timing.comparisonMs += elapsed;
            auto pQuantity = tpch16::EvalGreaterThan(
                cc, cQuantityLimit, cQuantity, &elapsed);
            timing.comparisonMs += elapsed;

            auto pDate = tpch16::EvalAnd(
                cc, pShipBegin, pShipEnd, &elapsed);
            timing.predicateAndMs += elapsed;
            auto pDiscount = tpch16::EvalAnd(
                cc, pDiscountLow, pDiscountHigh, &elapsed);
            timing.predicateAndMs += elapsed;
            auto pFour = tpch16::EvalAnd(
                cc, pDate, pDiscount, &elapsed);
            timing.predicateAndMs += elapsed;
            auto predicate = tpch16::EvalAnd(
                cc, pFour, pQuantity, &elapsed);
            timing.predicateAndMs += elapsed;
            predicate = cc->EvalMult(predicate, ptActiveResults);

            auto product = tpch16::EvalMultiply(
                cc, cPrice, cDiscount, &elapsed);
            timing.multiplicationMs += elapsed;

            begin = Clock::now();
            auto broadcast = BroadcastPredicate(cc, predicate);
            broadcast = cc->EvalMult(broadcast, ptActiveWords);
            auto selected =
                tpch16::EvalAnd(cc, product, broadcast, &elapsed);
            timing.selectionMs += Milliseconds(begin);

            begin = Clock::now();
            auto batchTotal = SumBatchRevenue(
                cc, selected, ptBitWeights,
                ptResultBitMask);
            if (!encryptedGrandTotal) {
                encryptedGrandTotal = batchTotal;
            } else {
                encryptedGrandTotal =
                    cc->EvalAdd(encryptedGrandTotal, batchTotal);
            }
            timing.aggregationMs += Milliseconds(begin);

            const double currentBatchMs =
                std::chrono::duration<double, std::milli>(
                    Clock::now() - batchEvaluationBegin)
                    .count();
            batchEvaluationMs.push_back(currentBatchMs);
            std::cout << std::fixed << std::setprecision(3)
                      << "batch " << (batch + 1) << "/" << batches
                      << " (" << activeRecords << " records)"
                      << " comparison-to-aggregation: "
                      << currentBatchMs << " ms\n"
                      << "  levels: product=" << product->GetLevel()
                      << ", selected=" << selected->GetLevel()
                      << ", batchTotal=" << batchTotal->GetLevel()
                      << '\n'
                      << std::flush;
        }

        begin = Clock::now();
        Plaintext decrypted;
        cc->Decrypt(keys.secretKey, encryptedGrandTotal, &decrypted);
        decrypted->SetLength(1);
        const double actualRevenue =
            decrypted->GetCKKSPackedValue()[0].real();
        timing.decryptionMs = Milliseconds(begin);

        const double absoluteError =
            std::abs(actualRevenue -
                     static_cast<double>(expectedRevenue));
        const double relativeError =
            absoluteError /
            std::max(1.0, static_cast<double>(expectedRevenue));
        const double queryMs =
            std::accumulate(batchEvaluationMs.begin(),
                            batchEvaluationMs.end(), 0.0);

        std::cout << std::fixed << std::setprecision(3)
                  << "\nTiming (ms)\n"
                  << "  GenCryptoContext       : " << timing.contextMs << '\n'
                  << "  KeyGen                 : " << timing.keyGenMs << '\n'
                  << "  evaluation-key Gen     : " << timing.evalKeyGenMs << '\n'
                  << "  input construction     : " << timing.inputBuildMs << '\n'
                  << "  encryption             : " << timing.encryptionMs << '\n'
                  << "  Boolean comparisons    : " << timing.comparisonMs << '\n'
                  << "  Boolean AND            : " << timing.predicateAndMs << '\n'
                  << "  Boolean multiplier     : " << timing.multiplicationMs << '\n'
                  << "  predicate selection    : " << timing.selectionMs << '\n'
                  << "  Rotate/Add aggregation : " << timing.aggregationMs << '\n'
                  << "  query evaluation       : " << queryMs << '\n'
                  << "  decrypt                : " << timing.decryptionMs << '\n'
                  << "\nBatch computation total\n"
                  << "  batches                : " << batches << '\n'
                  << "  total (no encryption)  : " << queryMs << " ms\n"
                  << "\nQuery result\n"
                  << "  encrypted revenue : " << actualRevenue << '\n'
                  << "  plaintext revenue : " << expectedRevenue << '\n'
                  << "  absolute error    : " << absoluteError << '\n'
                  << "  relative error    : " << relativeError << '\n';

        return relativeError < 0.01 ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
