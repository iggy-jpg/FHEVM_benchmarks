#include "Autohog_Cop.hpp"
#include "HomOP.h"
#include "timer.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace TFHEpp;
using boost::multiprecision::cpp_int;

namespace {

struct EpochStats {
    int serial_ok = 0;
    int parallel_ok = 0;
    int serial_parallel_same = 0;
    int total = 0;
    std::vector<double> serial_ms;
    std::vector<double> parallel_ms;
};

static void print_usage(const char *prog)
{
    std::cerr
        << "Usage: " << prog << " <mul|div|both> <bit-width> [epochs] "
        << "[--no-epoch-detail|--show-epochs]\n\n"
        << "Examples:\n"
        << "  " << prog << " mul 4 1 --show-epochs\n"
        << "  " << prog << " div 4 2 --show-epochs\n"
        << "  " << prog << " both 4 1 --no-epoch-detail\n";
}

static int max_threads()
{
#ifdef _OPENMP
    return std::max(1, omp_get_max_threads());
#else
    return 1;
#endif
}

static cpp_int mask_bits(int width)
{
    if (width <= 0) return 0;
    return (cpp_int(1) << width) - 1;
}

static std::vector<uint8_t> uint_to_bits_msb(cpp_int value, int width)
{
    value &= mask_bits(width);

    std::vector<uint8_t> bits(width);
    for (int i = 0; i < width; ++i) {
        const int shift = width - 1 - i;
        bits[i] = static_cast<uint8_t>(((value >> shift) & 1) != 0);
    }
    return bits;
}

static cpp_int bits_to_uint_msb(const std::vector<uint8_t> &bits)
{
    cpp_int value = 0;
    for (uint8_t bit : bits) {
        value <<= 1;
        value |= (bit & 1u);
    }
    return value;
}

static std::string bits_to_string(const std::vector<uint8_t> &bits)
{
    std::string out;
    out.reserve(bits.size());
    for (uint8_t bit : bits) {
        out.push_back((bit & 1u) ? '1' : '0');
    }
    return out;
}

static cpp_int random_uint(int width, std::default_random_engine &rng)
{
    cpp_int value = 0;
    int produced = 0;

    while (produced < width) {
        value <<= 32;
        value |= static_cast<unsigned int>(rng());
        produced += 32;
    }

    return value & mask_bits(width);
}

static std::vector<TLWE<lvl2param>> encrypt_uint_msb(
    cpp_int value,
    int width,
    const SecretKey &sk)
{
    const auto bits = uint_to_bits_msb(value, width);

    std::vector<TLWE<lvl2param>> out(bits.size());
    for (size_t i = 0; i < bits.size(); ++i) {
        out[i] = bootsSymEncrypt<lvl2param>({bits[i]}, sk)[0];
    }
    return out;
}

static std::vector<uint8_t> decrypt_bits_msb(
    const std::vector<TLWE<lvl2param>> &bits,
    const SecretKey &sk)
{
    std::vector<uint8_t> plain(bits.size());
    for (size_t i = 0; i < bits.size(); ++i) {
        plain[i] = bootsSymDecrypt<lvl2param>({bits[i]}, sk)[0] & 1u;
    }
    return plain;
}

static double avg_ms(const std::vector<double> &times)
{
    if (times.empty()) return 0.0;
    return std::accumulate(times.begin(), times.end(), 0.0) /
        static_cast<double>(times.size());
}

static bool same_bits(
    const std::vector<uint8_t> &lhs,
    const std::vector<uint8_t> &rhs)
{
    return lhs == rhs;
}

static void print_bool(const char *label, bool ok)
{
    std::cout << label << "=" << (ok ? "PASS" : "FAIL");
}

static bool run_mul_epoch(
    int epoch,
    int epochs,
    int width,
    const cpp_int &a_plain,
    const cpp_int &b_plain,
    const SecretKey &sk,
    const EvalKey &ek,
    bool show_detail,
    EpochStats &stats)
{
    const int product_width = width * 2;
    const cpp_int expected = (a_plain * b_plain) & mask_bits(product_width);

    const auto a_ct = encrypt_uint_msb(a_plain, width, sk);
    const auto b_ct = encrypt_uint_msb(b_plain, width, sk);

    std::vector<TLWE<lvl2param>> serial_ct;
    std::vector<TLWE<lvl2param>> parallel_ct;
    double serial_ms = 0.0;
    double parallel_ms = 0.0;

    {
        AutoTimer timer(&serial_ms);
        serial_ct = HomMul(a_ct, b_ct, width, ek);
    }
    {
        AutoTimer timer(&parallel_ms);
        parallel_ct = HomMul_pal(a_ct, b_ct, width, ek, max_threads());
    }

    if (static_cast<int>(serial_ct.size()) != product_width ||
        static_cast<int>(parallel_ct.size()) != product_width) {
        throw std::runtime_error("HomMul must return a 2W-bit full product");
    }

    const auto serial_bits = decrypt_bits_msb(serial_ct, sk);
    const auto parallel_bits = decrypt_bits_msb(parallel_ct, sk);
    const cpp_int serial_plain = bits_to_uint_msb(serial_bits);
    const cpp_int parallel_plain = bits_to_uint_msb(parallel_bits);

    const bool serial_ok = (serial_plain == expected);
    const bool parallel_ok = (parallel_plain == expected);
    const bool sp_same = same_bits(serial_bits, parallel_bits);

    ++stats.total;
    stats.serial_ok += serial_ok ? 1 : 0;
    stats.parallel_ok += parallel_ok ? 1 : 0;
    stats.serial_parallel_same += sp_same ? 1 : 0;
    stats.serial_ms.push_back(serial_ms);
    stats.parallel_ms.push_back(parallel_ms);

    if (show_detail) {
        std::cout << "\n================= MUL EPOCH " << epoch << "/" << epochs
                  << " =================\n";
        std::cout << "input_width=" << width
                  << " output_width=" << product_width
                  << " semantics=full_2W_bits\n";
        std::cout << "A bits              : " << bits_to_string(uint_to_bits_msb(a_plain, width)) << "\n";
        std::cout << "B bits              : " << bits_to_string(uint_to_bits_msb(b_plain, width)) << "\n";
        std::cout << "Expected bits       : " << bits_to_string(uint_to_bits_msb(expected, product_width)) << "\n";
        std::cout << "HomMul bits         : " << bits_to_string(serial_bits) << "\n";
        std::cout << "HomMul_pal bits     : " << bits_to_string(parallel_bits) << "\n";
        std::cout << "A plaintext         : " << a_plain << "\n";
        std::cout << "B plaintext         : " << b_plain << "\n";
        std::cout << "A * B plaintext     : " << expected << "\n";
        std::cout << "HomMul decrypt      : " << serial_plain << "\n";
        std::cout << "HomMul_pal decrypt  : " << parallel_plain << "\n";
        std::cout << std::fixed << std::setprecision(3)
                  << "time_ms             : serial=" << serial_ms
                  << " parallel=" << parallel_ms << "\n";
        print_bool("serial_vs_plain", serial_ok);
        std::cout << "  ";
        print_bool("parallel_vs_plain", parallel_ok);
        std::cout << "  ";
        print_bool("serial_vs_parallel", sp_same);
        std::cout << "\n";
    }

    return serial_ok && parallel_ok && sp_same;
}

static bool run_div_epoch(
    int epoch,
    int epochs,
    int width,
    const cpp_int &a_plain,
    const cpp_int &b_plain,
    const SecretKey &sk,
    const EvalKey &ek,
    bool show_detail,
    EpochStats &stats)
{
    cpp_int expected_q = 0;
    cpp_int expected_r = a_plain & mask_bits(width);
    if (b_plain != 0) {
        expected_q = a_plain / b_plain;
        expected_r = a_plain % b_plain;
    }

    const auto a_ct = encrypt_uint_msb(a_plain, width, sk);
    const auto b_ct = encrypt_uint_msb(b_plain, width, sk);

    HomDivResult serial_ct;
    HomDivResult parallel_ct;
    double serial_ms = 0.0;
    double parallel_ms = 0.0;

    {
        AutoTimer timer(&serial_ms);
        serial_ct = HomDiv(a_ct, b_ct, width, ek);
    }
    {
        AutoTimer timer(&parallel_ms);
        parallel_ct = HomDiv_pal(a_ct, b_ct, width, ek, max_threads());
    }

    if (static_cast<int>(serial_ct.quotient.size()) != width ||
        static_cast<int>(serial_ct.remainder.size()) != width ||
        static_cast<int>(parallel_ct.quotient.size()) != width ||
        static_cast<int>(parallel_ct.remainder.size()) != width) {
        throw std::runtime_error("HomDiv must return W-bit quotient and W-bit remainder");
    }

    const auto serial_q_bits = decrypt_bits_msb(serial_ct.quotient, sk);
    const auto serial_r_bits = decrypt_bits_msb(serial_ct.remainder, sk);
    const auto parallel_q_bits = decrypt_bits_msb(parallel_ct.quotient, sk);
    const auto parallel_r_bits = decrypt_bits_msb(parallel_ct.remainder, sk);

    const cpp_int serial_q_plain = bits_to_uint_msb(serial_q_bits);
    const cpp_int serial_r_plain = bits_to_uint_msb(serial_r_bits);
    const cpp_int parallel_q_plain = bits_to_uint_msb(parallel_q_bits);
    const cpp_int parallel_r_plain = bits_to_uint_msb(parallel_r_bits);

    const bool serial_ok =
        (serial_q_plain == expected_q) && (serial_r_plain == expected_r);
    const bool parallel_ok =
        (parallel_q_plain == expected_q) && (parallel_r_plain == expected_r);
    const bool sp_same =
        same_bits(serial_q_bits, parallel_q_bits) &&
        same_bits(serial_r_bits, parallel_r_bits);

    ++stats.total;
    stats.serial_ok += serial_ok ? 1 : 0;
    stats.parallel_ok += parallel_ok ? 1 : 0;
    stats.serial_parallel_same += sp_same ? 1 : 0;
    stats.serial_ms.push_back(serial_ms);
    stats.parallel_ms.push_back(parallel_ms);

    if (show_detail) {
        std::cout << "\n================= DIV EPOCH " << epoch << "/" << epochs
                  << " =================\n";
        std::cout << "width=" << width
                  << " semantics=quotient_W_bits_remainder_W_bits";
        if (b_plain == 0) std::cout << " div0=q0_rA";
        std::cout << "\n";
        std::cout << "A bits              : " << bits_to_string(uint_to_bits_msb(a_plain, width)) << "\n";
        std::cout << "B bits              : " << bits_to_string(uint_to_bits_msb(b_plain, width)) << "\n";
        std::cout << "Expected Q bits     : " << bits_to_string(uint_to_bits_msb(expected_q, width)) << "\n";
        std::cout << "Expected R bits     : " << bits_to_string(uint_to_bits_msb(expected_r, width)) << "\n";
        std::cout << "HomDiv Q bits       : " << bits_to_string(serial_q_bits) << "\n";
        std::cout << "HomDiv R bits       : " << bits_to_string(serial_r_bits) << "\n";
        std::cout << "HomDiv_pal Q bits   : " << bits_to_string(parallel_q_bits) << "\n";
        std::cout << "HomDiv_pal R bits   : " << bits_to_string(parallel_r_bits) << "\n";
        std::cout << "A plaintext         : " << a_plain << "\n";
        std::cout << "B plaintext         : " << b_plain << "\n";
        std::cout << "Expected quotient   : " << expected_q << "\n";
        std::cout << "Expected remainder  : " << expected_r << "\n";
        std::cout << "HomDiv decrypt Q    : " << serial_q_plain << "\n";
        std::cout << "HomDiv decrypt R    : " << serial_r_plain << "\n";
        std::cout << "HomDiv_pal decrypt Q: " << parallel_q_plain << "\n";
        std::cout << "HomDiv_pal decrypt R: " << parallel_r_plain << "\n";
        std::cout << std::fixed << std::setprecision(3)
                  << "time_ms             : serial=" << serial_ms
                  << " parallel=" << parallel_ms << "\n";
        print_bool("serial_vs_plain", serial_ok);
        std::cout << "  ";
        print_bool("parallel_vs_plain", parallel_ok);
        std::cout << "  ";
        print_bool("serial_vs_parallel", sp_same);
        std::cout << "\n";
    }

    return serial_ok && parallel_ok && sp_same;
}

static void print_summary(
    const std::string &op,
    int width,
    const EpochStats &stats)
{
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\n== " << op << " SUMMARY ==\n";
    std::cout << "width=" << width
              << " epochs=" << stats.total
              << " threads=" << max_threads() << "\n";
    std::cout << "serial_vs_plain     : " << stats.serial_ok << "/" << stats.total << "\n";
    std::cout << "parallel_vs_plain   : " << stats.parallel_ok << "/" << stats.total << "\n";
    std::cout << "serial_vs_parallel  : " << stats.serial_parallel_same << "/" << stats.total << "\n";
    std::cout << "avg_time_ms         : serial=" << avg_ms(stats.serial_ms)
              << " parallel=" << avg_ms(stats.parallel_ms) << "\n";
    std::cout << "SUMMARY op=" << op
              << " width=" << width
              << " epochs=" << stats.total
              << " serial_correct=" << stats.serial_ok << "/" << stats.total
              << " parallel_correct=" << stats.parallel_ok << "/" << stats.total
              << " serial_parallel_same=" << stats.serial_parallel_same << "/" << stats.total
              << " avg_serial_ms=" << avg_ms(stats.serial_ms)
              << " avg_parallel_ms=" << avg_ms(stats.parallel_ms)
              << "\n";
}

} // namespace

int main(int argc, char **argv)
{
    try {
        if (argc < 3) {
            print_usage(argv[0]);
            return 1;
        }

        const std::string op_arg = argv[1];
        if (op_arg != "mul" && op_arg != "div" && op_arg != "both") {
            throw std::runtime_error("op must be mul, div, or both");
        }

        const int width = std::stoi(argv[2]);
        if (width <= 0 || width > 256) {
            throw std::runtime_error("bit-width must be in [1, 256]");
        }

        int epochs = 1;
        if (argc >= 4 && argv[3][0] != '-') {
            epochs = std::max(1, std::stoi(argv[3]));
        }

        bool show_detail = true;
        for (int i = 3; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--no-epoch-detail") {
                show_detail = false;
            } else if (arg == "--show-epochs") {
                show_detail = true;
            }
        }

#ifdef _OPENMP
        omp_set_num_threads(max_threads());
#endif

        std::cout << "--- Initializing Keys ---\n";
        SecretKey sk;
        EvalKey ek;
        ek.emplacebkfft<lvl02param>(sk);
        ek.emplaceiksk<lvl20param>(sk);
        init_constants(sk);

        std::default_random_engine rng(std::random_device{}());
        EpochStats mul_stats;
        EpochStats div_stats;
        bool all_ok = true;

        for (int ep = 1; ep <= epochs; ++ep) {
            const cpp_int a_plain = random_uint(width, rng);
            cpp_int b_plain = random_uint(width, rng);

            if ((op_arg == "div" || op_arg == "both") && ep == epochs) {
                b_plain = 0;
            }

            if (op_arg == "mul" || op_arg == "both") {
                all_ok = run_mul_epoch(
                    ep, epochs, width, a_plain, b_plain, sk, ek,
                    show_detail, mul_stats) && all_ok;
            }

            if (op_arg == "div" || op_arg == "both") {
                all_ok = run_div_epoch(
                    ep, epochs, width, a_plain, b_plain, sk, ek,
                    show_detail, div_stats) && all_ok;
            }
        }

        if (op_arg == "mul" || op_arg == "both") {
            print_summary("mul", width, mul_stats);
        }
        if (op_arg == "div" || op_arg == "both") {
            print_summary("div", width, div_stats);
        }

        return all_ok ? 0 : 2;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        print_usage(argv[0]);
        return 1;
    }
}
