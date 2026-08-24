

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

struct Stat {
    std::vector<double> ser, par;
    int ok = 0, total = 0;
};

int main(int argc, char **argv)
{
    try {
        const int W = (argc >= 2) ? std::stoi(argv[1]) : 64;
        const int epochs = (argc >= 3) ? std::stoi(argv[2]) : 2;
        const std::string only = (argc >= 4) ? argv[3] : "all";
        auto want = [&](const char* n){ return only == "all" || only == n; };
        const int T = max_threads();

#ifdef _OPENMP
        omp_set_num_threads(T);
#endif
        std::cout << "--- Initializing Keys (W=" << W << ", threads=" << T << ") ---\n";
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

        Stat add, sub, sel, mulp, mull, divp, divq;

        for (int ep = 0; ep < epochs; ++ep) {
            const cpp_int a = rnd(W), b = rnd(W);
            const cpp_int k = rnd(W) | (cpp_int(1) << (W - 1));  
            const auto A = enc(a, W, sk), B = enc(b, W, sk);
            const auto kb = to_bits_msb(k, W);
            const auto S1 = bootsSymEncrypt<lvl2param>({1}, sk)[0];
            double t = 0;
            auto meas = [&](const char* nm, Stat& st, auto fn_ser, auto fn_par, auto ref) {
                if (!want(nm)) return;
                { double tt = 0; auto r = fn_ser(tt); st.ser.push_back(tt); st.ok += ref(r); st.total++; }
                { double tt = 0; auto r = fn_par(tt); st.par.push_back(tt); st.ok += ref(r); st.total++; }
            };
            const cpp_int sumv = (a + b) & mask_bits(W);
            const cpp_int subv = (a - b) & mask_bits(W);
            const cpp_int mulv = a * k;
            const cpp_int mullv = (a * k) & mask_bits(W);

            meas("add", add,
                 [&](double& tt){ HomAddResult r; { AutoTimer x(&tt); r = HomAdd_unsigned(A,B,W,ek,false); } return r; },
                 [&](double& tt){ HomAddResult r; { AutoTimer x(&tt); r = HomAdd_unsigned_pal(A,B,W,ek,false); } return r; },
                 [&](const HomAddResult& r){ return dec(r.SUM, sk) == sumv; });
            meas("sub", sub,
                 [&](double& tt){ HomAddResult r; { AutoTimer x(&tt); r = HomSub_unsigned(A,B,W,ek,false); } return r; },
                 [&](double& tt){ HomAddResult r; { AutoTimer x(&tt); r = HomSub_unsigned_pal(A,B,W,ek,false); } return r; },
                 [&](const HomAddResult& r){ return dec(r.SUM, sk) == subv; });
            meas("select", sel,
                 [&](double& tt){ std::vector<TLWE<lvl2param>> r; { AutoTimer x(&tt); r = HomSelect(A,B,S1,W,ek,1); } return r; },
                 [&](double& tt){ std::vector<TLWE<lvl2param>> r; { AutoTimer x(&tt); r = HomSelect(A,B,S1,W,ek,T); } return r; },
                 [&](const std::vector<TLWE<lvl2param>>& r){ return dec(r, sk) == b; });
            meas("mul_plain_full", mulp,
                 [&](double& tt){ std::vector<TLWE<lvl2param>> r; { AutoTimer x(&tt); r = HomMulPlain(A,kb,W,ek); } return r; },
                 [&](double& tt){ std::vector<TLWE<lvl2param>> r; { AutoTimer x(&tt); r = HomMulPlain_pal(A,kb,W,ek,T); } return r; },
                 [&](const std::vector<TLWE<lvl2param>>& r){ return dec(r, sk) == mulv; });
            meas("mul_plain_low", mull,
                 [&](double& tt){ std::vector<TLWE<lvl2param>> r; { AutoTimer x(&tt); r = HomMulPlainLow(A,kb,W,ek); } return r; },
                 [&](double& tt){ std::vector<TLWE<lvl2param>> r; { AutoTimer x(&tt); r = HomMulPlainLow_pal(A,kb,W,ek,T); } return r; },
                 [&](const std::vector<TLWE<lvl2param>>& r){ return dec(r, sk) == mullv; });
            meas("div_plain_q", divq,
                 [&](double& tt){ std::vector<TLWE<lvl2param>> r; { AutoTimer x(&tt); r = HomDivPlainQuotient(A,kb,W,ek); } return r; },
                 [&](double& tt){ std::vector<TLWE<lvl2param>> r; { AutoTimer x(&tt); r = HomDivPlainQuotient_pal(A,kb,W,ek,T); } return r; },
                 [&](const std::vector<TLWE<lvl2param>>& r){ return dec(r,sk)==a/k; });
            meas("div_plain", divp,
                 [&](double& tt){ HomDivResult r; { AutoTimer x(&tt); r = HomDivPlain(A,kb,W,ek); } return r; },
                 [&](double& tt){ HomDivResult r; { AutoTimer x(&tt); r = HomDivPlain_pal(A,kb,W,ek,T); } return r; },
                 [&](const HomDivResult& r){ return dec(r.quotient,sk)==a/k && dec(r.remainder,sk)==a%k; });
            std::cout << "epoch " << ep + 1 << "/" << epochs << " done\n";
        }

        auto avg = [](const std::vector<double> &v) {
            return v.empty() ? 0.0
                             : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        };
        auto row = [&](const char *name, const Stat &s) {
            std::cout << std::fixed << std::setprecision(1)
                      << "PERF op=" << name << " W=" << W
                      << " ser_ms=" << avg(s.ser)
                      << " par_ms=" << avg(s.par)
                      << " speedup=" << (avg(s.par) > 0 ? avg(s.ser) / avg(s.par) : 0)
                      << " correct=" << s.ok << "/" << s.total << "\n";
        };
        row("add", add);
        row("sub", sub);
        row("select", sel);
        row("mul_plain_full", mulp);
        row("mul_plain_low", mull);
        row("div_plain", divp);
        row("div_plain_q", divq);
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
