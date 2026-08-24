
//

//   test_PERF_ALL <bit-width> [epochs] [op|all]


//


#include "Autohog_Cop.hpp"
#include "HomOP.h"
#include "timer.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <functional>
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

static int dec1(const TLWE<lvl2param> &c, const SecretKey &sk)
{
    return bootsSymDecrypt<lvl2param>({c}, sk)[0] & 1;
}


static cpp_int to_signed(cpp_int v, int w)
{
    if ((v >> (w - 1)) & 1) v -= (cpp_int(1) << w);
    return v;
}

int main(int argc, char **argv)
{
    try {
        const int W = (argc >= 2) ? std::stoi(argv[1]) : 64;
        const int epochs = (argc >= 3) ? std::stoi(argv[2]) : 2;
        const std::string only = (argc >= 4) ? argv[3] : "all";
        const int T = max_threads();
#ifdef _OPENMP
        omp_set_num_threads(T);
#endif

        
        
        struct Op {
            std::string name;
            std::string cls;
            std::function<bool(double &)> run;   
        };
        std::vector<Op> ops;

        if (only == "list") {
            std::cout << "add_u sub_u add_s sub_s mul mul_low div div_newton "
                         "gt_u ge_u lt_u le_u gt_s ge_s lt_s le_s eq ne "
                         "and or xor not select cast "
                         "mul_plain mul_plain_low div_plain div_plain_q "
                         "eq_plain gt_plain and_plain or_plain xor_plain select_plain\n";
            return 0;
        }

        std::cout << "--- Initializing Keys (W=" << W << ", threads=" << T << ") ---\n";
        SecretKey sk;
        EvalKey ek;
        ek.emplacebkfft<lvl02param>(sk);
        ek.emplaceiksk<lvl20param>(sk);
        init_constants(sk);

        // Fixed inputs make planner-mode and backend A/B comparisons
        // reproducible; ciphertext encryption randomness remains independent.
        std::default_random_engine rng(0xC0FFEEu);
        auto rnd = [&](int w) {
            cpp_int v = 0;
            for (int p = 0; p < w; p += 32) { v <<= 32; v |= (unsigned)rng(); }
            return v & mask_bits(w);
        };

        struct Acc { std::vector<double> ms; int ok = 0, total = 0; std::string cls; };
        std::vector<std::pair<std::string, Acc>> results;
        auto slot = [&](const std::string &n, const std::string &cls) -> Acc & {
            for (auto &p : results) if (p.first == n) return p.second;
            results.push_back({n, Acc{}});
            results.back().second.cls = cls;
            return results.back().second;
        };

        // One unrecorded pass initializes the OpenMP worker pool, thread-local
        // FFTW plans and caches. Reported values therefore represent steady
        // operator execution rather than planner/startup cost.
        for (int ep = -1; ep < epochs; ++ep) {
            const bool warmup = ep < 0;
            const cpp_int a = rnd(W), b = rnd(W);
            const cpp_int k = rnd(W) | (cpp_int(1) << (W - 1));   
            const cpp_int x = rnd(W), y = rnd(W);
            const auto A = enc(a, W, sk), B = enc(b, W, sk);
            const auto kb = to_bits_msb(k, W);
            const auto S1 = bootsSymEncrypt<lvl2param>({1}, sk)[0];
            const cpp_int as = to_signed(a, W), bs = to_signed(b, W);

            auto run = [&](const char *name, const char *cls, auto body) {
                if (only != "all" && only != name) return;
                double t = 0;
                const bool ok = body(t);
                if (warmup) return;
                Acc &acc = slot(name, cls);
                acc.ms.push_back(t);
                acc.total++;
                acc.ok += ok ? 1 : 0;
                if (!ok) std::cout << "  FAIL: " << name << "\n";
            };

            
            run("add_u", "ct", [&](double &t) {
                HomAddResult r; { AutoTimer x2(&t); r = HomAdd_unsigned_pal(A, B, W, ek, false); }
                return dec(r.SUM, sk) == ((a + b) & mask_bits(W)); });
            run("sub_u", "ct", [&](double &t) {
                HomAddResult r; { AutoTimer x2(&t); r = HomSub_unsigned_pal(A, B, W, ek, false); }
                return dec(r.SUM, sk) == ((a - b) & mask_bits(W)); });
            run("add_s", "ct", [&](double &t) {
                HomAddResult r; { AutoTimer x2(&t); r = HomAdd_pal(A, B, W, ek, false, T); }
                return dec(r.SUM, sk) == ((a + b) & mask_bits(W)); });
            run("sub_s", "ct", [&](double &t) {
                HomAddResult r; { AutoTimer x2(&t); r = HomSub_pal(A, B, W, ek, false, T); }
                return dec(r.SUM, sk) == ((a - b) & mask_bits(W)); });
            run("mul", "ct", [&](double &t) {
                std::vector<TLWE<lvl2param>> r; { AutoTimer x2(&t); r = HomMul_pal(A, B, W, ek, T); }
                return dec(r, sk) == a * b; });
            run("mul_low", "ct", [&](double &t) {
                std::vector<TLWE<lvl2param>> r; { AutoTimer x2(&t); r = HomMulCompressorLow_pal(A, B, W, ek, T); }
                return dec(r, sk) == ((a * b) & mask_bits(W)); });
            run("div", "ct", [&](double &t) {
                HomDivResult r; { AutoTimer x2(&t); r = HomDiv_pal(A, B, W, ek, T); }
                if (b == 0) return dec(r.quotient, sk) == 0 && dec(r.remainder, sk) == a;
                return dec(r.quotient, sk) == a / b && dec(r.remainder, sk) == a % b; });
            run("div_newton", "ct", [&](double &t) {
                HomDivResult r; { AutoTimer x2(&t); r = HomDivNewton(A, B, W, ek, T); }
                if (b == 0) return dec(r.quotient, sk) == 0 && dec(r.remainder, sk) == a;
                return dec(r.quotient, sk) == a / b && dec(r.remainder, sk) == a % b; });

            run("gt_u", "ct", [&](double &t) { TLWE<lvl2param> r;
                { AutoTimer x2(&t); r = HomGT_unsigned_pal(A, B, W, ek, T); } return dec1(r, sk) == (a > b); });
            run("ge_u", "ct", [&](double &t) { TLWE<lvl2param> r;
                { AutoTimer x2(&t); r = HomGE_unsigned_pal(A, B, W, ek, T); } return dec1(r, sk) == (a >= b); });
            run("lt_u", "ct", [&](double &t) { TLWE<lvl2param> r;
                { AutoTimer x2(&t); r = HomLT_unsigned_pal(A, B, W, ek, T); } return dec1(r, sk) == (a < b); });
            run("le_u", "ct", [&](double &t) { TLWE<lvl2param> r;
                { AutoTimer x2(&t); r = HomLE_unsigned_pal(A, B, W, ek, T); } return dec1(r, sk) == (a <= b); });
            run("gt_s", "ct", [&](double &t) { TLWE<lvl2param> r;
                { AutoTimer x2(&t); r = HomGT_pal(A, B, W, ek, T); } return dec1(r, sk) == (as > bs); });
            run("ge_s", "ct", [&](double &t) { TLWE<lvl2param> r;
                { AutoTimer x2(&t); r = HomGE_pal(A, B, W, ek, T); } return dec1(r, sk) == (as >= bs); });
            run("lt_s", "ct", [&](double &t) { TLWE<lvl2param> r;
                { AutoTimer x2(&t); r = HomLT_pal(A, B, W, ek, T); } return dec1(r, sk) == (as < bs); });
            run("le_s", "ct", [&](double &t) { TLWE<lvl2param> r;
                { AutoTimer x2(&t); r = HomLE_pal(A, B, W, ek, T); } return dec1(r, sk) == (as <= bs); });
            run("eq", "ct", [&](double &t) { TLWE<lvl2param> r;
                { AutoTimer x2(&t); r = HomEQ_pal(A, B, W, ek, T); } return dec1(r, sk) == (a == b); });
            run("ne", "ct", [&](double &t) { TLWE<lvl2param> r;
                { AutoTimer x2(&t); r = HomNE_pal(A, B, W, ek, T); } return dec1(r, sk) == (a != b); });

            run("and", "ct", [&](double &t) { std::vector<TLWE<lvl2param>> r;
                { AutoTimer x2(&t); r = HomAND(A, B, W, ek, T); } return dec(r, sk) == (a & b); });
            run("or", "ct", [&](double &t) { std::vector<TLWE<lvl2param>> r;
                { AutoTimer x2(&t); r = HomOR(A, B, W, ek, T); } return dec(r, sk) == (a | b); });
            run("xor", "ct", [&](double &t) { std::vector<TLWE<lvl2param>> r;
                { AutoTimer x2(&t); r = HomXOR(A, B, W, ek, T); } return dec(r, sk) == (a ^ b); });
            run("not", "free", [&](double &t) { std::vector<TLWE<lvl2param>> r;
                { AutoTimer x2(&t); r = HomNOT(A, W, ek); } return dec(r, sk) == (~a & mask_bits(W)); });
            run("select", "ct", [&](double &t) { std::vector<TLWE<lvl2param>> r;
                { AutoTimer x2(&t); r = HomSelect(A, B, S1, W, ek, T); } return dec(r, sk) == b; });
            run("cast", "free", [&](double &t) { std::vector<TLWE<lvl2param>> r;
                { AutoTimer x2(&t); r = HomCast(A, W, 2 * W, false); } return dec(r, sk) == a; });

            
            run("mul_plain", "plain", [&](double &t) { std::vector<TLWE<lvl2param>> r;
                { AutoTimer x2(&t); r = HomMulPlain_pal(A, kb, W, ek, T); } return dec(r, sk) == a * k; });
            run("mul_plain_low", "plain", [&](double &t) { std::vector<TLWE<lvl2param>> r;
                { AutoTimer x2(&t); r = HomMulPlainLow_pal(A, kb, W, ek, T); }
                return dec(r, sk) == ((a * k) & mask_bits(W)); });
            run("div_plain", "plain", [&](double &t) { HomDivResult r;
                { AutoTimer x2(&t); r = HomDivPlain_pal(A, kb, W, ek, T); }
                return dec(r.quotient, sk) == a / k && dec(r.remainder, sk) == a % k; });
            run("div_plain_q", "plain", [&](double &t) { std::vector<TLWE<lvl2param>> r;
                { AutoTimer x2(&t); r = HomDivPlainQuotient_pal(A, kb, W, ek, T); }
                return dec(r, sk) == a / k; });
            run("eq_plain", "plain", [&](double &t) { TLWE<lvl2param> r;
                { AutoTimer x2(&t); r = HomEQPlain_pal(A, kb, W, ek, T); } return dec1(r, sk) == (a == k); });
            run("gt_plain", "plain", [&](double &t) { TLWE<lvl2param> r;
                { AutoTimer x2(&t); r = HomGTPlain_pal(A, kb, W, ek, T); } return dec1(r, sk) == (a > k); });
            run("and_plain", "free", [&](double &t) { std::vector<TLWE<lvl2param>> r;
                { AutoTimer x2(&t); r = HomANDPlain(A, kb, W); } return dec(r, sk) == (a & k); });
            run("or_plain", "free", [&](double &t) { std::vector<TLWE<lvl2param>> r;
                { AutoTimer x2(&t); r = HomORPlain(A, kb, W); } return dec(r, sk) == (a | k); });
            run("xor_plain", "free", [&](double &t) { std::vector<TLWE<lvl2param>> r;
                { AutoTimer x2(&t); r = HomXORPlain(A, kb, W); } return dec(r, sk) == (a ^ k); });
            run("select_plain", "free", [&](double &t) { std::vector<TLWE<lvl2param>> r;
                { AutoTimer x2(&t); r = HomSelectPlain(to_bits_msb(x, W), to_bits_msb(y, W), S1, W); }
                return dec(r, sk) == y; });

            if (warmup)
                std::cout << "warmup done\n";
            else
                std::cout << "epoch " << ep + 1 << "/" << epochs << " done\n";
        }

        auto avg = [](const std::vector<double> &v) {
            return v.empty() ? 0.0 : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        };
        for (const auto &p : results) {
            std::cout << std::fixed << std::setprecision(2)
                      << "PERF op=" << p.first << " W=" << W
                      << " ms=" << avg(p.second.ms)
                      << " pbs_class=" << p.second.cls
                      << " correct=" << p.second.ok << "/" << p.second.total << "\n";
        }
        int bad = 0;
        for (const auto &p : results) bad += (p.second.ok != p.second.total);
        return bad == 0 ? 0 : 2;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
