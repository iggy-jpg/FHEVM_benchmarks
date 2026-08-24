#include "fhe_integer_bits.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>

using namespace arcedb;
using namespace arcedb_integer;

template <class Function>
double Measure(Function&& function)
{
    const auto begin = std::chrono::steady_clock::now();
    function();
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - begin).count();
}

int main()
{
    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    std::mt19937_64 rng(0x415243454442ULL);
    int failures = 0;

    for (int width : {4, 8, 16, 32, 64, 128, 256}) {
        cpp_int a = 0, b = 0;
        for (int offset = 0; offset < width; offset += 64) {
            a |= cpp_int(rng()) << offset;
            b |= cpp_int(rng()) << offset;
        }
        const cpp_int mask = (cpp_int(1) << width) - 1;
        a &= mask; b &= mask;
        const auto ca = Encrypt(a, width, sk), cb = Encrypt(b, width, sk);
        auto reportWord = [&](const char* op, double ms, const Word& actual, cpp_int expected) {
            const bool ok = Decrypt(actual, sk) == (expected & mask);
            failures += !ok;
            std::cout << std::fixed << std::setprecision(3) << "PERF op=" << op
                      << " W=" << width << " ms=" << ms
                      << " correct=" << (ok ? "yes" : "no") << '\n';
        };
        auto reportBit = [&](const char* op, double ms, const Bit& actual, bool expected) {
            const bool ok = Decrypt(actual, sk) == expected;
            failures += !ok;
            std::cout << std::fixed << std::setprecision(3) << "PERF op=" << op
                      << " W=" << width << " ms=" << ms
                      << " correct=" << (ok ? "yes" : "no") << '\n';
        };

        Word word;
        double ms = Measure([&] { word = Add(ca, cb, ek); });
        reportWord("add", ms, word, a + b);
        ms = Measure([&] { word = Sub(ca, cb, ek); });
        reportWord("sub", ms, word, a - b);
        ms = Measure([&] { word = Multiply(ca, cb, ek); });
        reportWord("mult", ms, word, a * b);

        Comparison comparison;
        ms = Measure([&] { comparison = Compare(ca, cb, ek); });
        // One comparison pass computes all five predicates; report the shared
        // latency for each requested comparison operation.
        reportBit("gt", ms, comparison.gt, a > b);
        reportBit("lt", ms, comparison.lt, a < b);
        reportBit("eq", ms, comparison.eq, a == b);
        reportBit("ge", ms, comparison.ge, a >= b);
        reportBit("le", ms, comparison.le, a <= b);
    }
    return failures == 0 ? 0 : 2;
}
