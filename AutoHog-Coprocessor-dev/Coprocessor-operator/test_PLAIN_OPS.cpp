

#include "Autohog_Cop.hpp"
#include "fhe_key.h"
#include "HomOP.h"
#include "timer.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <iomanip>
#include <iostream>
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

static int dec1(const TLWE<lvl2param> &c, const SecretKey &sk)
{
    return bootsSymDecrypt<lvl2param>({c}, sk)[0] & 1;
}

int main(int argc, char **argv)
{
    try {
        const int W = (argc >= 2) ? std::stoi(argv[1]) : 16;
        const int N = (argc >= 3) ? std::stoi(argv[2]) : 2 * max_threads();

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

        bool all_ok = true;
        int checks = 0, pass = 0;
        auto chk = [&](bool ok, const char *what) {
            ++checks;
            pass += ok ? 1 : 0;
            all_ok = all_ok && ok;
            if (!ok) std::cout << "FAIL: " << what << "\n";
        };

        
        {
            for (int trial = 0; trial < 50; ++trial) {
                const cpp_int v = rnd(W);
                Bytes be;
                const int nbytes = (W + 7) / 8;
                for (int i = nbytes - 1; i >= 0; --i)
                    be.push_back(static_cast<uint8_t>((v >> (i * 8)) & 0xff));
                chk(FheKey::plain_bytes_to_bits(be, W) == to_bits_msb(v, W),
                    "plain_bytes_to_bits");
            }
        }

        
        {
            const cpp_int a = rnd(W);
            std::vector<cpp_int> ks = {a, cpp_int(0), mask_bits(W),
                                       a == 0 ? cpp_int(1) : cpp_int(a - 1),
                                       cpp_int((a + 1) & mask_bits(W)), rnd(W)};
            const auto A = enc(a, W, sk);
            for (const auto &k : ks) {
                const auto kb = to_bits_msb(k, W);
                chk(dec1(HomEQPlain_pal(A, kb, W, ek, max_threads()), sk) == (a == k), "EQPlain");
                chk(dec1(HomGTPlain_pal(A, kb, W, ek, max_threads()), sk) == (a > k), "GTPlain");
                chk(dec1(HomGEPlain_pal(A, kb, W, ek, max_threads()), sk) == (a >= k), "GEPlain");
                chk(dec1(HomLTPlain_pal(A, kb, W, ek, max_threads()), sk) == (a < k), "LTPlain");
                chk(dec1(HomLEPlain_pal(A, kb, W, ek, max_threads()), sk) == (a <= k), "LEPlain");
            }
        }

        
        {
            const cpp_int a = rnd(W), k = rnd(W), x = rnd(W), y = rnd(W);
            const auto A = enc(a, W, sk);
            const auto kb = to_bits_msb(k, W);
            chk(dec(HomANDPlain(A, kb, W), sk) == (a & k), "ANDPlain");
            chk(dec(HomORPlain(A, kb, W), sk) == (a | k), "ORPlain");
            chk(dec(HomXORPlain(A, kb, W), sk) == (a ^ k), "XORPlain");
            const auto s1 = bootsSymEncrypt<lvl2param>({1}, sk)[0];
            const auto s0 = bootsSymEncrypt<lvl2param>({0}, sk)[0];
            chk(dec(HomSelectPlain(to_bits_msb(x, W), to_bits_msb(y, W), s1, W), sk) == y,
                "SelectPlain sel=1");
            chk(dec(HomSelectPlain(to_bits_msb(x, W), to_bits_msb(y, W), s0, W), sk) == x,
                "SelectPlain sel=0");
        }

        
        {
            const cpp_int k = rnd(W) | 1;
            const auto kb = to_bits_msb(k, W);
            std::vector<cpp_int> as(N);
            std::vector<std::vector<TLWE<lvl2param>>> As(N);
            for (int i = 0; i < N; ++i) { as[i] = rnd(W); As[i] = enc(as[i], W, sk); }

            double t_batch = 0, t_seq = 0;
            std::vector<std::vector<TLWE<lvl2param>>> prods;
            {
                AutoTimer t(&t_batch);
                prods = HomMulPlainBatch(As, kb, W, ek, max_threads());
            }
            {
                AutoTimer t(&t_seq);
                for (int i = 0; i < N; ++i)
                    (void)HomMulPlain_pal(As[i], kb, W, ek, max_threads());
            }
            for (int i = 0; i < N; ++i)
                chk(dec(prods[i], sk) == as[i] * k, "MulPlainBatch");

            double t_dbatch = 0;
            std::vector<HomDivResult> divs;
            {
                AutoTimer t(&t_dbatch);
                divs = HomDivPlainBatch(As, kb, W, ek, max_threads());
            }
            for (int i = 0; i < N; ++i) {
                chk(dec(divs[i].quotient, sk) == as[i] / k, "DivPlainBatch q");
                chk(dec(divs[i].remainder, sk) == as[i] % k, "DivPlainBatch r");
            }
            
            auto qs = HomDivPlainQuotientBatch(As, kb, W, ek, max_threads());
            for (int i = 0; i < N; ++i)
                chk(dec(qs[i], sk) == as[i] / k, "DivPlainQuotientBatch");
            for (const cpp_int &dv : {cpp_int(0), cpp_int(1), cpp_int(2),
                                      cpp_int(1) << (W - 1), mask_bits(W)}) {
                const auto dvb = to_bits_msb(dv, W);
                const cpp_int expq = (dv == 0) ? cpp_int(0) : cpp_int(as[0] / dv);
                chk(dec(HomDivPlainQuotient(As[0], dvb, W, ek), sk) == expq,
                    "DivPlainQuotient special");
            }
            
            double t_low = 0;
            std::vector<std::vector<TLWE<lvl2param>>> lows;
            {
                AutoTimer t(&t_low);
                lows = HomMulPlainLowBatch(As, kb, W, ek, max_threads());
            }
            for (int i = 0; i < N; ++i)
                chk(dec(lows[i], sk) == ((as[i] * k) & mask_bits(W)), "MulPlainLowBatch");

            
            double t_lowseq = 0;
            {
                AutoTimer t(&t_lowseq);
                for (int i = 0; i < N; ++i)
                    (void)HomMulPlainLow_pal(As[i], kb, W, ek, max_threads());
            }

            std::cout << std::fixed << std::setprecision(1)
                      << "batch N=" << N << " W=" << W
                      << " mul_low_seq_ms=" << t_lowseq
                      << " mul_batch_ms=" << t_batch
                      << " mul_seq_ms=" << t_seq
                      << " speedup=" << (t_batch > 0 ? t_seq / t_batch : 0)
                      << " mul_low_batch_ms=" << t_low
                      << " low_vs_full=" << (t_low > 0 ? t_batch / t_low : 0)
                      << " div_batch_ms=" << t_dbatch << "\n";
        }

        std::cout << "SUMMARY op=plain_ops width=" << W
                  << " checks=" << pass << "/" << checks
                  << " threads=" << max_threads() << "\n";
        return all_ok ? 0 : 2;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
