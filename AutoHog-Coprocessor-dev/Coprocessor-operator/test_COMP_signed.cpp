

#include "Autohog_Cop.hpp"
#include "HomOP.h"

#include <iostream>
#include <random>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace TFHEpp;

static std::vector<TLWE<lvl2param>> enc(long long v, int w, const SecretKey &sk)
{
    std::vector<TLWE<lvl2param>> out(w);
    for (int i = 0; i < w; ++i) {
        uint8_t bit = static_cast<uint8_t>((v >> (w - 1 - i)) & 1);
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
    const int epochs = (argc >= 3) ? std::stoi(argv[2]) : 4;

    std::cout << "--- Initializing Keys ---\n";
    SecretKey sk;
    EvalKey ek;
    ek.emplacebkfft<lvl02param>(sk);
    ek.emplaceiksk<lvl20param>(sk);
    init_constants(sk);

    std::default_random_engine rng(std::random_device{}());
    const long long lo = -(1LL << (W - 1)), hi = (1LL << (W - 1)) - 1;
    std::uniform_int_distribution<long long> dist(lo, hi);

    int ok = 0, total = 0;
    for (int ep = 0; ep < epochs; ++ep) {
        long long a = dist(rng), b = dist(rng);
        if (ep == 0) { a = lo; b = -1; }        // INT_MIN vs -1
        if (ep == 1) { a = hi; b = lo; }        // MAX vs MIN
        if (ep == 2) { b = a; }                 
        const auto A = enc(a, W, sk), B = enc(b, W, sk);

        const bool pass =
            dec1(HomGT(A, B, W, ek), sk) == (a > b) &&
            dec1(HomLT(A, B, W, ek), sk) == (a < b) &&
            dec1(HomGE(A, B, W, ek), sk) == (a >= b) &&
            dec1(HomLE(A, B, W, ek), sk) == (a <= b) &&
            dec1(HomGT_pal(A, B, W, ek, 4), sk) == (a > b) &&
            dec1(HomLT_pal(A, B, W, ek, 4), sk) == (a < b) &&
            dec1(HomGE_pal(A, B, W, ek, 4), sk) == (a >= b) &&
            dec1(HomLE_pal(A, B, W, ek, 4), sk) == (a <= b);

        ++total; ok += pass ? 1 : 0;
        std::cout << "epoch " << ep + 1 << " a=" << a << " b=" << b
                  << " -> " << (pass ? "PASS" : "FAIL") << "\n";
    }
    std::cout << "SUMMARY op=comp_signed width=" << W << " correct=" << ok << "/" << total << "\n";
    return (ok == total) ? 0 : 2;
}
