// Core unsigned integer performance matrix.
// Usage: test_PERF_CORE [epochs] [comma-separated-widths]
// Defaults: one measured run at 4,8,16,32,64,128,256 bits.
#include "HomOP.h"
#include "workload_common.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace TFHEpp;
using namespace autohog_workload;

static std::vector<int> parse_widths(const std::string &text)
{
    std::vector<int> widths;
    std::stringstream input(text);
    std::string item;
    while (std::getline(input, item, ',')) {
        const int width = std::stoi(item);
        if (width <= 0) throw std::runtime_error("bit width must be positive");
        widths.push_back(width);
    }
    if (widths.empty()) throw std::runtime_error("no bit widths supplied");
    return widths;
}

template <class Function>
static double measure(Function &&function)
{
    const auto begin = std::chrono::steady_clock::now();
    function();
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

static cpp_int random_uint(std::mt19937_64 &rng, int width)
{
    cpp_int value = 0;
    for (int generated = 0; generated < width; generated += 64) {
        value <<= std::min(64, width - generated);
        value |= cpp_int(rng()) & bit_mask(std::min(64, width - generated));
    }
    return value & bit_mask(width);
}

int main(int argc, char **argv)
{
    try {
        const int epochs = argc > 1 ? std::stoi(argv[1]) : 1;
        const auto widths = parse_widths(
            argc > 2 ? argv[2] : "4,8,16,32,64,128,256");
        if (epochs <= 0) throw std::runtime_error("epochs must be positive");

        const int threads = thread_count();
        std::cout << "Initializing TFHE keys (threads=" << threads << ")...\n";
        SecretKey sk;
        EvalKey ek;
        ek.emplacebkfft<lvl02param>(sk);
        ek.emplaceiksk<lvl20param>(sk);
        init_constants(sk);
        std::mt19937_64 rng(0x4155544f484f47ULL);

        auto print = [](const char *op, int width, double total_ms, int ok, int n) {
            std::cout << std::fixed << std::setprecision(3)
                      << "PERF op=" << op << " W=" << width
                      << " ms=" << total_ms / n
                      << " correct=" << ok << "/" << n << '\n';
        };

        int failures = 0;
        for (const int width : widths) {
            double add_ms = 0, sub_ms = 0, mul_ms = 0, gt_ms = 0;
            double lt_ms = 0, eq_ms = 0, ge_ms = 0, le_ms = 0;
            int add_ok = 0, sub_ok = 0, mul_ok = 0, gt_ok = 0;
            int lt_ok = 0, eq_ok = 0, ge_ok = 0, le_ok = 0;

            for (int epoch = 0; epoch < epochs; ++epoch) {
                const cpp_int a = random_uint(rng, width);
                const cpp_int b = random_uint(rng, width);
                const auto encrypted_a = encrypt_uint(a, width, sk);
                const auto encrypted_b = encrypt_uint(b, width, sk);

                HomAddResult add_result;
                add_ms += measure([&] {
                    add_result = HomAdd_unsigned_pal(
                        encrypted_a, encrypted_b, width, ek, false, threads);
                });
                add_ok += decrypt_uint(add_result.SUM, sk) ==
                          ((a + b) & bit_mask(width));

                HomAddResult sub_result;
                sub_ms += measure([&] {
                    sub_result = HomSub_unsigned_pal(
                        encrypted_a, encrypted_b, width, ek, false, threads);
                });
                sub_ok += decrypt_uint(sub_result.SUM, sk) ==
                          ((a - b) & bit_mask(width));

                std::vector<TLWE<lvl2param>> mul_result;
                mul_ms += measure([&] {
                    mul_result = HomMul_pal(
                        encrypted_a, encrypted_b, width, ek, threads);
                });
                mul_ok += decrypt_uint(mul_result, sk) == a * b;

                TLWE<lvl2param> result;
                gt_ms += measure([&] { result = HomGT_unsigned_pal(
                    encrypted_a, encrypted_b, width, ek, threads); });
                gt_ok += decrypt_bit(result, sk) == (a > b);
                lt_ms += measure([&] { result = HomLT_unsigned_pal(
                    encrypted_a, encrypted_b, width, ek, threads); });
                lt_ok += decrypt_bit(result, sk) == (a < b);
                eq_ms += measure([&] { result = HomEQ_pal(
                    encrypted_a, encrypted_b, width, ek, threads); });
                eq_ok += decrypt_bit(result, sk) == (a == b);
                ge_ms += measure([&] { result = HomGE_unsigned_pal(
                    encrypted_a, encrypted_b, width, ek, threads); });
                ge_ok += decrypt_bit(result, sk) == (a >= b);
                le_ms += measure([&] { result = HomLE_unsigned_pal(
                    encrypted_a, encrypted_b, width, ek, threads); });
                le_ok += decrypt_bit(result, sk) == (a <= b);
            }

            print("add", width, add_ms, add_ok, epochs);
            print("sub", width, sub_ms, sub_ok, epochs);
            print("mult", width, mul_ms, mul_ok, epochs);
            print("gt", width, gt_ms, gt_ok, epochs);
            print("lt", width, lt_ms, lt_ok, epochs);
            print("eq", width, eq_ms, eq_ok, epochs);
            print("ge", width, ge_ms, ge_ok, epochs);
            print("le", width, le_ms, le_ok, epochs);
            failures += add_ok != epochs || sub_ok != epochs || mul_ok != epochs ||
                        gt_ok != epochs || lt_ok != epochs || eq_ok != epochs ||
                        ge_ok != epochs || le_ok != epochs;
        }
        return failures == 0 ? 0 : 2;
    } catch (const std::exception &error) {
        std::cerr << "test_PERF_CORE: " << error.what() << '\n';
        return 1;
    }
}
