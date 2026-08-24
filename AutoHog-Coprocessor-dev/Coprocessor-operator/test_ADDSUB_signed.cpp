

#include "Autohog_Cop.hpp"
#include "HomOP.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <iostream>
#include <random>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace TFHEpp;
using boost::multiprecision::cpp_int;

static cpp_int mask_bits(int w) { return (cpp_int(1) << w) - 1; }

static std::vector<TLWE<lvl2param>> enc(cpp_int v, int w, const SecretKey &sk)
{
    v &= mask_bits(w);
    std::vector<TLWE<lvl2param>> out(w);
    for (int i = 0; i < w; ++i) {
        uint8_t bit = static_cast<uint8_t>(((v >> (w - 1 - i)) & 1) != 0);
        out[i] = bootsSymEncrypt<lvl2param>({bit}, sk)[0];
    }
    return out;
}

static cpp_int dec(const std::vector<TLWE<lvl2param>> &ct, const SecretKey &sk)
{
    cpp_int v = 0;
    for (const auto &c : ct) v = (v << 1) | (bootsSymDecrypt<lvl2param>({c}, sk)[0] & 1);
    return v;
}

int main(int argc, char **argv)
{
    if (argc < 2) { std::cerr << "Usage: " << argv[0] << " <W> [epochs]\n"; return 1; }
    const int W = std::stoi(argv[1]);
    const int epochs = (argc >= 3) ? std::stoi(argv[2]) : 3;

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

    int ok = 0, total = 0;
    for (int ep = 0; ep < epochs; ++ep) {
        const cpp_int a = rnd(W), b = rnd(W);
        const auto A = enc(a, W, sk), B = enc(b, W, sk);
        const cpp_int esum = (a + b) & mask_bits(W);
        const cpp_int ediff = (a - b) & mask_bits(W);

        const cpp_int s_ser = dec(HomAdd(A, B, W, ek, false).SUM, sk);
        const cpp_int s_par = dec(HomAdd_pal(A, B, W, ek, false).SUM, sk);
        const cpp_int d_ser = dec(HomSub(A, B, W, ek, false).SUM, sk);
        const cpp_int d_par = dec(HomSub_pal(A, B, W, ek, false).SUM, sk);

        const bool pass = (s_ser == esum) && (s_par == esum) &&
                          (d_ser == ediff) && (d_par == ediff);
        ++total; ok += pass ? 1 : 0;
        std::cout << "epoch " << ep + 1 << ": add_ser=" << (s_ser == esum ? "PASS" : "FAIL")
                  << " add_par=" << (s_par == esum ? "PASS" : "FAIL")
                  << " sub_ser=" << (d_ser == ediff ? "PASS" : "FAIL")
                  << " sub_par=" << (d_par == ediff ? "PASS" : "FAIL") << "\n";
    }
    std::cout << "SUMMARY op=addsub_signed width=" << W << " correct=" << ok << "/" << total << "\n";
    return (ok == total) ? 0 : 2;
}
