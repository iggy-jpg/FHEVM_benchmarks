

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
    std::vector<TLWE<lvl2param>> out(w);
    for (int i = 0; i < w; ++i) {
        uint8_t bit = static_cast<uint8_t>(((v >> (w - 1 - i)) & 1) != 0);
        out[i] = bootsSymEncrypt<lvl2param>({bit}, sk)[0];
    }
    return out;
}

static int dec1(const TLWE<lvl2param> &c, const SecretKey &sk)
{
    return bootsSymDecrypt<lvl2param>({c}, sk)[0] & 1;
}

int main(int argc, char **argv)
{
    if (argc < 2) { std::cerr << "Usage: " << argv[0] << " <W> [epochs]\n"; return 1; }
    const int W = std::stoi(argv[1]);
    const int epochs = (argc >= 3) ? std::stoi(argv[2]) : 5;

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
        cpp_int a = rnd(W), b;
        if (ep == 0) b = a;                                    
        else if (ep == 1) b = a ^ (cpp_int(1) << (W - 1));     
        else if (ep == 2) b = a ^ 1;                           
        else if (ep == 3) { a = 0; b = mask_bits(W); }         
        else b = rnd(W);

        const auto A = enc(a, W, sk), B = enc(b, W, sk);
        const bool pass =
            dec1(HomEQ(A, B, W, ek), sk) == (a == b) &&
            dec1(HomEQ_pal(A, B, W, ek, 4), sk) == (a == b) &&
            dec1(HomNE(A, B, W, ek), sk) == (a != b);
        ++total; ok += pass ? 1 : 0;
        std::cout << "epoch " << ep + 1 << " -> " << (pass ? "PASS" : "FAIL") << "\n";
    }
    std::cout << "SUMMARY op=eq_ops width=" << W << " correct=" << ok << "/" << total << "\n";
    return (ok == total) ? 0 : 2;
}
