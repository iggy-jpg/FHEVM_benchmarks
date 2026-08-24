#include "HomOP.h"
#include "HomUtil.h"
#include "timer.h"

#include <tfhe++.hpp>
#include <random>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <sstream>
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

        
        uint64_t global_correct_bits = 0;
        uint64_t global_total_bits = 0;

        
        for (size_t wi = 0; wi < widths.size(); ++wi)
        {
            const int W_req = widths[wi];
            std::cout << "\n================= BIT-WIDTH " << W_req
                      << " (epochs=" << epochs << ") =================\n";

            
            std::vector<double> t_or_ms_vec;

            for (int ep = 0; ep < epochs; ++ep)
            {
                if (show_epoch_detail)
                {
                    std::cout << "\n----- EPOCH " << (ep + 1) << "/" << epochs << " -----\n";
                }

                
                SecretKey sk;
                EvalKey ek;
                {
                    ek.emplacebkfft<lvl02param>(sk);
                    ek.emplaceiksk<lvl20param>(sk);
                    init_constants(sk);
                }

                
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
                for (int i = 0; i < W; i++)
                {
                    Apl[i] = bit01(eng);
                    Bpl[i] = bit01(eng);
                }

                
                std::vector<TLWE<lvl2param>> A(W), B(W);
                for (int i = 0; i < W; i++)
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

                
                std::vector<TLWE<lvl2param>> or_cts;
                double t_ms = 0.0;
                {
                    AutoTimer timer(&t_ms);
                    or_cts = HomOR(A, B, W, ek); 
                }
                t_or_ms_vec.push_back(t_ms);

                
                std::vector<uint8_t> ORpl(W);
                for (int i = 0; i < W; i++)
                    ORpl[i] = tlweSymDecrypt<lvl2param>(or_cts[i], sk.key.get<lvl2param>());

                std::vector<uint8_t> refOR(W);
                for (int i = 0; i < W; i++)
                    refOR[i] = (Apl[i] | Bpl[i]);

                bool ok = true;
                uint64_t correct_bits = 0;
                for (int i = 0; i < W; i++)
                {
                    if (ORpl[i] == refOR[i])
                        ++correct_bits;
                    else
                        ok = false;
                }

                
                global_correct_bits += correct_bits;
                global_total_bits += static_cast<uint64_t>(W);

                if (show_epoch_detail)
                {
                    auto print_bits = [](const char *tag, const std::vector<uint8_t> &v)
                    {
                        std::cout << tag;
                        for (uint8_t b : v)
                            std::cout << int(b) << ' ';
                        std::cout << '\n';
                    };
                    std::cout << "Bitwise OR (A | B) Bit width = " << W << "\n";
                    print_bits("A   (msb->lsb): ", Apl);
                    print_bits("B   (msb->lsb): ", Bpl);
                    std::cout << "-----------------------------\n";
                    print_bits("A | B (cipher): ", ORpl);
                    print_bits("A | B (plain) : ", refOR);
                    std::cout << "Correct bits this epoch: " << correct_bits << " / " << W << "\n";
                    if (!ok)
                        std::cout << "Warning: OR result mismatch!\n";
                    std::cout << "Circuit execution: " << t_ms << " ms\n";
                }
                else
                {
                    if (!ok)
                        std::cout << "[Epoch " << (ep + 1) << "] OR mismatch!\n";
                }
            } // end epochs

            
            if (epochs > 1)
            {
                std::cout << "\n---- SUMMARY (Bit-width " << W_req
                          << ", average over " << epochs << " epochs) ----\n";
                std::cout << "HomOR execution AVG: " << avg(t_or_ms_vec) << " ms\n";
            }
            else
            {
                
                std::cout << "\nHomOR execution AVG: " << avg(t_or_ms_vec) << " ms\n";
            }
        } // end widths

        
        if (global_total_bits > 0)
        {
            double acc = 100.0 * double(global_correct_bits) / double(global_total_bits);
            std::cout << "\n================= GLOBAL ACCURACY =================\n";
            std::cout << "Correct bits: " << global_correct_bits << " / " << global_total_bits
                      << "  (" << acc << "%)\n";
        }

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return -1;
    }
}