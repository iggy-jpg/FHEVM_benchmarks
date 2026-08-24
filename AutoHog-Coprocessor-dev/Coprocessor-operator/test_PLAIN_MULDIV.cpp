

#include "Autohog_Cop.hpp"
#include "HomOP.h"
#include "timer.h"

#include <boost/multiprecision/cpp_int.hpp>

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

static cpp_int mask_bits(int w) { return (cpp_int(1) << w) - 1; }

static std::vector<uint8_t> to_bits_msb(cpp_int v, int w)
{
    v &= mask_bits(w);
    std::vector<uint8_t> bits(w);
    for (int i = 0; i < w; ++i)
        bits[i] = static_cast<uint8_t>(((v >> (w - 1 - i)) & 1) != 0);
    return bits;
}

static std::vector<TLWE<lvl2param>> enc(cpp_int v, int w, const SecretKey &sk)
{
    const auto bits = to_bits_msb(v, w);
    std::vector<TLWE<lvl2param>> out(w);
    for (int i = 0; i < w; ++i)
        out[i] = bootsSymEncrypt<lvl2param>({bits[i]}, sk)[0];
    return out;
}

static cpp_int dec(const std::vector<TLWE<lvl2param>> &ct, const SecretKey &sk)
{
    cpp_int v = 0;
    for (const auto &c : ct)
        v = (v << 1) | (bootsSymDecrypt<lvl2param>({c}, sk)[0] & 1);
    return v;
}

int main(int argc, char **argv)
{
    try {
        if (argc < 2) {
            std::cerr << "Usage: " << argv[0] << " <bit-width> [epochs] [--no-ref]\n";
            return 1;
        }
        const int W = std::stoi(argv[1]);
        int epochs = (argc >= 3 && argv[2][0] != '-') ? std::stoi(argv[2]) : 6;
        bool run_ref = true;
        for (int i = 2; i < argc; ++i)
            if (std::string(argv[i]) == "--no-ref") run_ref = false;

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
        auto rnd = [&](int w) {
            cpp_int v = 0;
            for (int p = 0; p < w; p += 32) { v <<= 32; v |= (unsigned)rng(); }
            return v & mask_bits(w);
        };

        
        std::vector<cpp_int> specials = {0, 1, 2, 3, 7,
                                         cpp_int(1) << (W - 1),
                                         (cpp_int(1) << (W - 1)) + 1,
                                         mask_bits(W)};

        int mul_ok = 0, div_ok = 0, total = 0;
        std::vector<double> mulp_ms, divp_ms, mulc_ms, divc_ms;
        bool all_ok = true;

        for (int ep = 0; ep < epochs; ++ep) {
            const cpp_int a = rnd(W);
            const cpp_int k = (ep < (int)specials.size()) ? specials[ep] : rnd(W);

            const auto A = enc(a, W, sk);
            const auto kb = to_bits_msb(k, W);

            double t_mp = 0, t_dp = 0, t_mc = 0, t_dc = 0;
            std::vector<TLWE<lvl2param>> prod;
            HomDivResult dres;
            {
                AutoTimer t(&t_mp);
                prod = HomMulPlain_pal(A, kb, W, ek, max_threads());
            }
            {
                AutoTimer t(&t_dp);
                dres = HomDivPlain_pal(A, kb, W, ek, max_threads());
            }
            if (run_ref && ep == 0) {
                const auto K = enc(k, W, sk);
                { AutoTimer t(&t_mc); (void)HomMul_pal(A, K, W, ek, max_threads()); }
                { AutoTimer t(&t_dc); (void)HomDiv_pal(A, K, W, ek, max_threads()); }
                mulc_ms.push_back(t_mc);
                divc_ms.push_back(t_dc);
            }

            const cpp_int exp_p = a * k;
            cpp_int exp_q = 0, exp_r = a;
            if (k != 0) { exp_q = a / k; exp_r = a % k; }

            const bool m_ok = dec(prod, sk) == exp_p;
            const bool d_ok = dec(dres.quotient, sk) == exp_q &&
                              dec(dres.remainder, sk) == exp_r;
            ++total;
            mul_ok += m_ok ? 1 : 0;
            div_ok += d_ok ? 1 : 0;
            mulp_ms.push_back(t_mp);
            divp_ms.push_back(t_dp);
            all_ok = all_ok && m_ok && d_ok;

            std::cout << "epoch " << ep + 1 << " k=" << k
                      << " mul=" << (m_ok ? "PASS" : "FAIL")
                      << " div=" << (d_ok ? "PASS" : "FAIL")
                      << std::fixed << std::setprecision(1)
                      << " mul_plain_ms=" << t_mp << " div_plain_ms=" << t_dp << "\n";
        }

        auto avg = [](const std::vector<double> &v) {
            return v.empty() ? 0.0
                             : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        };
        std::cout << std::fixed << std::setprecision(3)
                  << "SUMMARY op=plain_muldiv width=" << W << " epochs=" << total
                  << " mul_correct=" << mul_ok << "/" << total
                  << " div_correct=" << div_ok << "/" << total
                  << " avg_mul_plain_ms=" << avg(mulp_ms)
                  << " avg_div_plain_ms=" << avg(divp_ms)
                  << " ref_mul_ct_ms=" << avg(mulc_ms)
                  << " ref_div_ct_ms=" << avg(divc_ms)
                  << " threads=" << max_threads() << "\n";
        return all_ok ? 0 : 2;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
