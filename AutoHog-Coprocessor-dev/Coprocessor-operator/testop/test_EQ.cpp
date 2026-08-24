#include "HomOP.h"
#include "HomUtil.h"
#include "timer.h"

#include <tfhe++.hpp>
#include <random>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <boost/multiprecision/cpp_int.hpp>


#include <thread>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace TFHEpp;
using namespace boost::multiprecision;


template <class T>
static double avg(const std::vector<T> &v)
{
    if (v.empty())
        return 0.0;
    double s = std::accumulate(v.begin(), v.end(), 0.0);
    return s / static_cast<double>(v.size());
}

static void print_usage(const char *prog)
{
    std::cerr << "Usage: " << prog << " <bit-width | comma-list | all> [epochs] "
              << "[--no-epoch-detail|--show-epochs]\n"
              << "  Examples:\n"
              << "    " << prog << " 128 100 --show-epochs\n"
              << "    " << prog << " 8,16,32,64,128 20 --no-epoch-detail\n"
              << "    " << prog << " all 50\n";
}

static std::vector<int> parse_widths(const std::string &s)
{
    if (s == "all" || s == "ALL")
        return {8, 16, 32, 64, 128};
    std::vector<int> out;
    if (s.find(',') == std::string::npos)
    {
        int w = std::atoi(s.c_str());
        if (w > 0)
            out.push_back(w);
        return out;
    }
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ','))
    {
        if (tok.empty())
            continue;
        int w = std::atoi(tok.c_str());
        if (w > 0)
            out.push_back(w);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

static TLWE<lvl2param> enc_bit(uint8_t b, const SecretKey &sk)
{
    return tlweSymEncrypt<lvl2param>(b ? lvl2param::μ : -lvl2param::μ,
                                     lvl2param::α,
                                     sk.key.get<lvl2param>());
}

static std::vector<TLWE<lvl2param>> encrypt_bits(const std::vector<uint8_t> &v,
                                                 const SecretKey &sk)
{
    std::vector<TLWE<lvl2param>> ct(v.size());
    for (size_t i = 0; i < v.size(); ++i)
        ct[i] = enc_bit(v[i], sk);
    return ct;
}

static void print_bits(const char *tag, const std::vector<uint8_t> &v)
{
    std::cout << tag;
    for (uint8_t b : v)
        std::cout << int(b) << ' ';
    std::cout << '\n';
}


int main(int argc, char *argv[])
{
    try
    {
        if (argc < 2)
        {
            print_usage(argv[0]);
            return 1;
        }

        
#ifdef _OPENMP
        {
            omp_set_dynamic(0);
            int hw = std::max(1u, std::thread::hardware_concurrency());
            omp_set_num_threads(hw);
        }
#endif

        
        std::vector<int> widths = parse_widths(argv[1]);
        if (widths.empty())
        {
            std::cerr << "Error: invalid <bit-width | comma-list | all>.\n";
            return 1;
        }

        
        int epochs = 1;
        if (argc >= 3 && argv[2][0] != '-')
        {
            epochs = std::max(1, std::atoi(argv[2]));
        }

        
        bool show_epoch_detail = true;
        for (int i = 2; i < argc; ++i)
        {
            std::string s = argv[i];
            if (s == "--no-epoch-detail")
                show_epoch_detail = false;
            else if (s == "--show-epochs")
                show_epoch_detail = true;
        }

        
        for (size_t wi = 0; wi < widths.size(); ++wi)
        {
            int W_req = widths[wi];
            std::cout << "\n================= BIT-WIDTH " << W_req
                      << " (epochs=" << epochs << ") =================\n";

            
            std::vector<double> t_ms_rand_ser, t_ms_equal_ser, t_ms_noteq_ser;
            std::vector<double> t_ms_rand_par, t_ms_equal_par, t_ms_noteq_par;

            
            
            size_t ser_total = 0, ser_correct = 0;
            size_t par_total = 0, par_correct = 0;

            for (int ep = 0; ep < epochs; ++ep)
            {
                if (show_epoch_detail)
                {
                    std::cout << "\n----- EPOCH " << (ep + 1) << "/" << epochs << " -----\n";
                }

                
                SecretKey sk;
                EvalKey ek;
                ek.emplacebkfft<lvl02param>(sk);
                ek.emplaceiksk<lvl20param>(sk);
                init_constants(sk);

                
                CT_ZERO = tlweSymEncrypt<lvl2param>(-lvl2param::μ, lvl2param::α, sk.key.get<lvl2param>());
                CT_ONE = tlweSymEncrypt<lvl2param>(+lvl2param::μ, lvl2param::α, sk.key.get<lvl2param>());

                
                int W = W_req;
                if (W > static_cast<int>(lvl2param::n))
                {
                    if (show_epoch_detail)
                    {
                        std::cerr << "[Warn] given W=" << W
                                  << " > lvl2param::n=" << lvl2param::n
                                  << ", truncated.\n";
                    }
                    W = static_cast<int>(lvl2param::n);
                }

                // RNG
                std::default_random_engine eng(std::random_device{}());
                std::uniform_int_distribution<int> bit01(0, 1);

                
                auto report_cmp = [&](const std::string &title,
                                      const std::vector<uint8_t> &Apl,
                                      const std::vector<uint8_t> &Bpl,
                                      double t_ser, double t_par,
                                      uint8_t eq_ser, uint8_t eq_par)
                {
                    
                    cpp_int A_s = bits_to_int_signed(Apl, W);
                    cpp_int B_s = bits_to_int_signed(Bpl, W);
                    bool truth = (A_s == B_s);

                    bool ser_ok = (eq_ser == static_cast<uint8_t>(truth));
                    bool par_ok = (eq_par == static_cast<uint8_t>(truth));
                    bool ser_par_same = (eq_ser == eq_par);

                    
                    ser_total += 1;
                    par_total += 1;
                    if (ser_ok)
                        ser_correct += 1;
                    if (par_ok)
                        par_correct += 1;

                    if (show_epoch_detail)
                    {
                        std::cout << title << "\n";
                        print_bits("A   (msb->lsb): ", Apl);
                        print_bits("B   (msb->lsb): ", Bpl);
                        std::cout << "-----------------------------\n";
                        std::cout << "A (signed)      = " << A_s << "\n";
                        std::cout << "B (signed)      = " << B_s << "\n";

                        
                        std::cout << "(A == B) [SER]  = " << int(eq_ser)
                                  << "  " << (ser_ok ? "✅" : "❌") << "  (truth=" << (truth ? 1 : 0) << ")\n";
                        std::cout << "(A == B) [PAR]  = " << int(eq_par)
                                  << "  " << (par_ok ? "✅" : "❌") << "  (truth=" << (truth ? 1 : 0) << ")\n";

                        std::cout << "SER vs PAR equal? " << (ser_par_same ? "YES ✅" : "NO ❌") << "\n";
                        std::cout << "Circuit execution [SER]: " << t_ser << " ms\n";
                        std::cout << "Circuit execution [PAR]: " << t_par << " ms\n\n";
                    }
                };

                
                std::vector<uint8_t> Apl(W), Bpl(W);
                for (int i = 0; i < W; ++i)
                {
                    Apl[i] = static_cast<uint8_t>(bit01(eng));
                    Bpl[i] = static_cast<uint8_t>(bit01(eng));
                }
                auto A = encrypt_bits(Apl, sk);
                auto B = encrypt_bits(Bpl, sk);

                TLWE<lvl2param> eq_ser1, eq_par1;
                double t_ser1 = 0.0, t_par1 = 0.0;
                {
                    AutoTimer timer(&t_ser1);
                    eq_ser1 = HomEQ(A, B, W, ek);
                }
                {
                    AutoTimer timer(&t_par1);
                    eq_par1 = HomEQ_pal(A, B, W, ek);
                }
                uint8_t eq1_ser = tlweSymDecrypt<lvl2param>(eq_ser1, sk.key.get<lvl2param>());
                uint8_t eq1_par = tlweSymDecrypt<lvl2param>(eq_par1, sk.key.get<lvl2param>());
                t_ms_rand_ser.push_back(t_ser1);
                t_ms_rand_par.push_back(t_par1);
                report_cmp("SignedInt CMP(A==B) Bit width = " + std::to_string(W) + "  (Random vs Random)",
                           Apl, Bpl, t_ser1, t_par1, eq1_ser, eq1_par);

                
                std::vector<uint8_t> B_eq_pl = Apl;
                auto B_eq = encrypt_bits(B_eq_pl, sk);

                TLWE<lvl2param> eq_ser2, eq_par2;
                double t_ser2 = 0.0, t_par2 = 0.0;
                {
                    AutoTimer timer(&t_ser2);
                    eq_ser2 = HomEQ(A, B_eq, W, ek);
                }
                {
                    AutoTimer timer(&t_par2);
                    eq_par2 = HomEQ_pal(A, B_eq, W, ek);
                }
                uint8_t eq2_ser = tlweSymDecrypt<lvl2param>(eq_ser2, sk.key.get<lvl2param>());
                uint8_t eq2_par = tlweSymDecrypt<lvl2param>(eq_par2, sk.key.get<lvl2param>());
                t_ms_equal_ser.push_back(t_ser2);
                t_ms_equal_par.push_back(t_par2);
                report_cmp("SignedInt CMP(A==B) Bit width = " + std::to_string(W) + "  (Guaranteed Equal)",
                           Apl, B_eq_pl, t_ser2, t_par2, eq2_ser, eq2_par);

                
                std::vector<uint8_t> B_ne_pl = Apl;
                if (!B_ne_pl.empty())
                    B_ne_pl.back() ^= 1; 
                auto B_ne = encrypt_bits(B_ne_pl, sk);

                TLWE<lvl2param> eq_ser3, eq_par3;
                double t_ser3 = 0.0, t_par3 = 0.0;
                {
                    AutoTimer timer(&t_ser3);
                    eq_ser3 = HomEQ(A, B_ne, W, ek);
                }
                {
                    AutoTimer timer(&t_par3);
                    eq_par3 = HomEQ_pal(A, B_ne, W, ek);
                }
                uint8_t eq3_ser = tlweSymDecrypt<lvl2param>(eq_ser3, sk.key.get<lvl2param>());
                uint8_t eq3_par = tlweSymDecrypt<lvl2param>(eq_par3, sk.key.get<lvl2param>());
                t_ms_noteq_ser.push_back(t_ser3);
                t_ms_noteq_par.push_back(t_par3);
                report_cmp("SignedInt CMP(A==B) Bit width = " + std::to_string(W) + "  (Guaranteed Not-Equal)",
                           Apl, B_ne_pl, t_ser3, t_par3, eq3_ser, eq3_par);

                
                double avg_ser_epoch = (t_ser1 + t_ser2 + t_ser3) / 3.0;
                double avg_par_epoch = (t_par1 + t_par2 + t_par3) / 3.0;
                if (show_epoch_detail)
                {
                    std::cout << "== Epoch " << (ep + 1) << " AVERAGE (W=" << W << ") ==\n";
                    std::cout << "HomEQ [SER] AVG: " << avg_ser_epoch << " ms\n";
                    std::cout << "HomEQ [PAR] AVG: " << avg_par_epoch << " ms\n";
                }
            } // end epochs

            
            if (epochs > 1)
            {
                std::cout << "\n---- SUMMARY (Bit-width " << W_req
                          << ", average over " << epochs << " epochs) ----\n";
                std::cout << "HomEQ [SER] (Random vs Random)      AVG: " << avg(t_ms_rand_ser) << " ms\n";
                std::cout << "HomEQ [PAR] (Random vs Random)      AVG: " << avg(t_ms_rand_par) << " ms\n";
                std::cout << "HomEQ [SER] (Guaranteed Equal)      AVG: " << avg(t_ms_equal_ser) << " ms\n";
                std::cout << "HomEQ [PAR] (Guaranteed Equal)      AVG: " << avg(t_ms_equal_par) << " ms\n";
                std::cout << "HomEQ [SER] (Guaranteed NotEqual)   AVG: " << avg(t_ms_noteq_ser) << " ms\n";
                std::cout << "HomEQ [PAR] (Guaranteed NotEqual)   AVG: " << avg(t_ms_noteq_par) << " ms\n";
            }

            
            auto pct = [](size_t ok, size_t tot) -> double
            {
                if (tot == 0)
                    return 0.0;
                return 100.0 * static_cast<double>(ok) / static_cast<double>(tot);
            };
            size_t overall_total = ser_total + par_total;  
            size_t overall_ok = ser_correct + par_correct; 

            std::cout << "\n==== ACCURACY (Bit-width " << W_req << ") ====\n";
            std::cout << "SER: " << ser_correct << "/" << ser_total
                      << " (" << pct(ser_correct, ser_total) << "%) "
                      << (ser_correct == ser_total ? "✅" : "❌") << "\n";
            std::cout << "PAR: " << par_correct << "/" << par_total
                      << " (" << pct(par_correct, par_total) << "%) "
                      << (par_correct == par_total ? "✅" : "❌") << "\n";
            std::cout << "ALL: " << overall_ok << "/" << overall_total
                      << " (" << pct(overall_ok, overall_total) << "%) "
                      << (overall_ok == overall_total ? "✅" : "❌") << "\n";

        } // end widths

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return -1;
    }
}