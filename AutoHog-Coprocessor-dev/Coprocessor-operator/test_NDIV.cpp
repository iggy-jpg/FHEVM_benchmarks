

#include "Autohog_Cop.hpp"
#include "HomOP.h"
#include "timer.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace TFHEpp;
using boost::multiprecision::cpp_int;

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
    return (cpp_int(1) << width) - 1;
}

static std::vector<uint8_t> uint_to_bits_msb(cpp_int value, int width)
{
    value &= mask_bits(width);
    std::vector<uint8_t> bits(width);
    for (int i = 0; i < width; ++i)
        bits[i] = static_cast<uint8_t>(((value >> (width - 1 - i)) & 1) != 0);
    return bits;
}

static cpp_int bits_to_uint_msb(const std::vector<uint8_t> &bits)
{
    cpp_int v = 0;
    for (uint8_t b : bits) v = (v << 1) | (b & 1u);
    return v;
}

static cpp_int random_uint(int width, std::default_random_engine &rng)
{
    cpp_int v = 0;
    for (int produced = 0; produced < width; produced += 32) {
        v <<= 32;
        v |= static_cast<unsigned int>(rng());
    }
    return v & mask_bits(width);
}

static std::vector<TLWE<lvl2param>> encrypt_uint_msb(
    cpp_int value, int width, const SecretKey &sk)
{
    const auto bits = uint_to_bits_msb(value, width);
    std::vector<TLWE<lvl2param>> out(bits.size());
    for (size_t i = 0; i < bits.size(); ++i)
        out[i] = bootsSymEncrypt<lvl2param>({bits[i]}, sk)[0];
    return out;
}

static std::vector<uint8_t> decrypt_bits_msb(
    const std::vector<TLWE<lvl2param>> &bits, const SecretKey &sk)
{
    std::vector<uint8_t> plain(bits.size());
    for (size_t i = 0; i < bits.size(); ++i)
        plain[i] = bootsSymDecrypt<lvl2param>({bits[i]}, sk)[0] & 1u;
    return plain;
}

int main(int argc, char **argv)
{
    try {
        if (argc < 2) {
            std::cerr << "Usage: " << argv[0]
                      << " <bit-width> [epochs] [--show-epochs|--no-epoch-detail] [--no-ref]\n";
            return 1;
        }
        const int width = std::stoi(argv[1]);
        if (width <= 0 || width > 256) throw std::runtime_error("width in [1,256]");
        int epochs = 1;
        if (argc >= 3 && argv[2][0] != '-') epochs = std::max(1, std::stoi(argv[2]));
        bool show_detail = true;
        bool run_ref = true;
        for (int i = 2; i < argc; ++i) {
            const std::string a = argv[i];
            if (a == "--no-epoch-detail") show_detail = false;
            if (a == "--show-epochs") show_detail = true;
            if (a == "--no-ref") run_ref = false;
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
        int newton_ok = 0, ref_ok = 0, total = 0;
        std::vector<double> newton_ms, ref_ms;
        bool all_ok = true;

        for (int ep = 1; ep <= epochs; ++ep) {
            cpp_int a = random_uint(width, rng);
            cpp_int b = random_uint(width, rng);
            if (ep == epochs) b = 0;  
            if (ep == 1) b &= mask_bits(std::max(1, width / 4));  

            cpp_int eq = 0, er = a;
            if (b != 0) { eq = a / b; er = a % b; }

            const auto a_ct = encrypt_uint_msb(a, width, sk);
            const auto b_ct = encrypt_uint_msb(b, width, sk);

            HomDivResult nres, rres;
            double n_ms = 0.0, r_ms = 0.0;
            {
                AutoTimer t(&n_ms);
                nres = HomDivNewton(a_ct, b_ct, width, ek, max_threads());
            }
            if (run_ref) {
                AutoTimer t(&r_ms);
                rres = HomDiv_pal(a_ct, b_ct, width, ek, max_threads());
            }

            const cpp_int nq = bits_to_uint_msb(decrypt_bits_msb(nres.quotient, sk));
            const cpp_int nr = bits_to_uint_msb(decrypt_bits_msb(nres.remainder, sk));
            const bool n_ok = (nq == eq) && (nr == er);
            bool r_ok = true;
            if (run_ref) {
                const cpp_int rq = bits_to_uint_msb(decrypt_bits_msb(rres.quotient, sk));
                const cpp_int rr = bits_to_uint_msb(decrypt_bits_msb(rres.remainder, sk));
                r_ok = (rq == eq) && (rr == er);
            }

            ++total;
            newton_ok += n_ok ? 1 : 0;
            ref_ok += r_ok ? 1 : 0;
            newton_ms.push_back(n_ms);
            if (run_ref) ref_ms.push_back(r_ms);
            all_ok = all_ok && n_ok && r_ok;

            if (show_detail) {
                std::cout << "EPOCH " << ep << "/" << epochs
                          << " A=" << a << " B=" << b
                          << " expQ=" << eq << " gotQ=" << nq
                          << " expR=" << er << " gotR=" << nr
                          << " newton=" << (n_ok ? "PASS" : "FAIL");
                if (run_ref) std::cout << " restoring=" << (r_ok ? "PASS" : "FAIL");
                std::cout << std::fixed << std::setprecision(1)
                          << " newton_ms=" << n_ms;
                if (run_ref) std::cout << " restoring_ms=" << r_ms;
                std::cout << "\n";
            }
        }

        const double n_avg = newton_ms.empty() ? 0.0
            : std::accumulate(newton_ms.begin(), newton_ms.end(), 0.0) / newton_ms.size();
        const double r_avg = ref_ms.empty() ? 0.0
            : std::accumulate(ref_ms.begin(), ref_ms.end(), 0.0) / ref_ms.size();
        std::cout << std::fixed << std::setprecision(3)
                  << "SUMMARY op=ndiv width=" << width
                  << " epochs=" << total
                  << " newton_correct=" << newton_ok << "/" << total
                  << " restoring_correct=" << ref_ok << "/" << total
                  << " avg_newton_ms=" << n_avg
                  << " avg_restoring_ms=" << r_avg
                  << " threads=" << max_threads() << "\n";
        return all_ok ? 0 : 2;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
