#include "HomOP.h"
#include "HomUtil.h"
#include "timer.h"

#include <tfhe++.hpp>
#include <random>
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <boost/multiprecision/cpp_int.hpp>

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

        
        long long global_total_ser = 0, global_correct_ser = 0;
        long long global_total_par = 0, global_correct_par = 0;

        
        for (size_t wi = 0; wi < widths.size(); ++wi)
        {
            int W_req = widths[wi];
            std::cout << "\n================= BIT-WIDTH " << W_req
                      << " (epochs=" << epochs << ") =================\n";

            
            std::vector<double> t_ser_ms_vec, t_par_ms_vec;
            int mismatches_ser = 0, mismatches_par = 0;

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

                
                std::default_random_engine eng(std::random_device{}());
                std::uniform_int_distribution<int> bit01(0, 1);

                std::vector<uint8_t> Apl(W), Bpl(W);
                for (int i = 0; i < W; ++i)
                {
                    Apl[i] = bit01(eng);
                    Bpl[i] = bit01(eng);
                }

                
                std::vector<TLWE<lvl2param>> A(W), B(W);
                for (int i = 0; i < W; ++i)
                {
                    A[i] = tlweSymEncrypt<lvl2param>(
                        Apl[i] ? lvl2param::μ : -lvl2param::μ,
                        lvl2param::α,
                        sk.key.get<lvl2param>());

                    B[i] = tlweSymEncrypt<lvl2param>(
                        Bpl[i] ? lvl2param::μ : -lvl2param::μ,
                        lvl2param::α,
                        sk.key.get<lvl2param>());
                }

                
                TLWE<lvl2param> ge_ser;
                double t_ser_ms = 0.0;
                {
                    AutoTimer timer(&t_ser_ms);
                    ge_ser = HomGE(A, B, W, ek);
                }
                t_ser_ms_vec.push_back(t_ser_ms);
                uint8_t ge_ser_pl = tlweSymDecrypt<lvl2param>(ge_ser, sk.key.get<lvl2param>()) ? 1 : 0;

                
                TLWE<lvl2param> ge_par;
                double t_par_ms = 0.0;
                {
                    AutoTimer timer(&t_par_ms);
                    ge_par = HomGE_pal(A, B, W, ek, -1); 
                }
                t_par_ms_vec.push_back(t_par_ms);
                uint8_t ge_par_pl = tlweSymDecrypt<lvl2param>(ge_par, sk.key.get<lvl2param>()) ? 1 : 0;

                
                cpp_int A_s = bits_to_int_signed(Apl, W);
                cpp_int B_s = bits_to_int_signed(Bpl, W);
                const uint8_t ref_ge = (A_s >= B_s) ? 1 : 0;

                
                bool ok_ser = (ge_ser_pl == ref_ge);
                bool ok_par = (ge_par_pl == ref_ge);

                global_total_ser += 1;
                global_total_par += 1;
                if (ok_ser)
                    global_correct_ser += 1;
                else
                    mismatches_ser += 1;
                if (ok_par)
                    global_correct_par += 1;
                else
                    mismatches_par += 1;

                
                if (show_epoch_detail)
                {
                    std::cout << "SignedInt GE (A >= B)  Bit width = " << W << "\n";
                    print_bits("A   (msb->lsb): ", Apl);
                    print_bits("B   (msb->lsb): ", Bpl);
                    std::cout << "-----------------------------\n";
                    std::cout << "A (signed) = " << A_s << "\n";
                    std::cout << "B (signed) = " << B_s << "\n";
                    std::cout << "[SER] GE = " << int(ge_ser_pl)
                              << "  | t = " << t_ser_ms << " ms"
                              << (ok_ser ? "" : "  (mismatch!)") << "\n";
                    std::cout << "[PAR] GE = " << int(ge_par_pl)
                              << "  | t = " << t_par_ms << " ms"
                              << (ok_par ? "" : "  (mismatch!)") << "\n";
                }
            } // end epochs

            
            const long long width_total = epochs;
            const long long width_correct_ser = width_total - mismatches_ser;
            const long long width_correct_par = width_total - mismatches_par;
            const double acc_ser = width_total ? (100.0 * width_correct_ser / width_total) : 0.0;
            const double acc_par = width_total ? (100.0 * width_correct_par / width_total) : 0.0;
            const double avg_ser = avg(t_ser_ms_vec);
            const double avg_par = avg(t_par_ms_vec);
            const double speedup = (avg_par > 0.0) ? (avg_ser / avg_par) : 0.0;

            std::cout << std::fixed << std::setprecision(3);
            std::cout << "\n---- SUMMARY (Bit-width " << W_req
                      << ", average over " << epochs << " epochs) ----\n";
            std::cout << "HomGE [SER] execution AVG: " << avg_ser << " ms\n";
            std::cout << "HomGE [PAR] execution AVG: " << avg_par << " ms\n";
            std::cout << "Speedup (SER/ PAR): " << std::setprecision(2) << speedup << "x\n";
            std::cout << std::setprecision(2);
            std::cout << "Accuracy [SER]: " << width_correct_ser << " / " << width_total
                      << "  (" << acc_ser << "%)" << (width_correct_ser == width_total ? " ✅" : "") << "\n";
            std::cout << "Accuracy [PAR]: " << width_correct_par << " / " << width_total
                      << "  (" << acc_par << "%)" << (width_correct_par == width_total ? " ✅" : "") << "\n";
        } // end widths

        
        const double global_acc_ser = global_total_ser ? (100.0 * global_correct_ser / global_total_ser) : 0.0;
        const double global_acc_par = global_total_par ? (100.0 * global_correct_par / global_total_par) : 0.0;
        std::cout << "\n================= GLOBAL ACCURACY =================\n";
        std::cout << "SER: " << global_correct_ser << " / " << global_total_ser
                  << "  (" << std::fixed << std::setprecision(2) << global_acc_ser << "%)"
                  << (global_correct_ser == global_total_ser ? " ✅" : "") << "\n";
        std::cout << "PAR: " << global_correct_par << " / " << global_total_par
                  << "  (" << std::fixed << std::setprecision(2) << global_acc_par << "%)"
                  << (global_correct_par == global_total_par ? " ✅" : "") << "\n";

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return -1;
    }
}