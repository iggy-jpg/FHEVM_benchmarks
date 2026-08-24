#include "HomOP.h"
#include "HomUtil.h"
#include "timer.h"

#include <random>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>

#include <tfhe++.hpp>
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
template <typename T>
static inline void deep_free(std::vector<T> &v) { std::vector<T>().swap(v); }

static void print_usage(const char *prog)
{
    std::cerr << "Usage: " << prog << " <bit-width> [epochs] [--no-epoch-detail|--show-epochs]\n";
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
        int W = std::atoi(argv[1]);
        if (W <= 0)
        {
            std::cerr << "Error: bit-width must be positive.\n";
            return 1;
        }

        int epochs = 1;
        if (argc >= 3 && argv[2][0] != '-')
            epochs = std::max(1, std::atoi(argv[2]));

        bool show_epoch_detail = true;
        for (int i = 2; i < argc; ++i)
        {
            std::string s = argv[i];
            if (s == "--no-epoch-detail")
                show_epoch_detail = false;
            else if (s == "--show-epochs")
                show_epoch_detail = true;
        }

        
        if (W > static_cast<int>(lvl2param::n))
        {
            std::cerr << "[Warn] given W=" << W << " > lvl2param::n=" << lvl2param::n
                      << ", truncated.\n";
            W = static_cast<int>(lvl2param::n);
        }

        
        SecretKey sk;
        EvalKey ek;
        ek.emplacebkfft<lvl02param>(sk);
        ek.emplaceiksk<lvl20param>(sk);
        init_constants(sk);

        
        std::default_random_engine eng(std::random_device{}());
        std::uniform_int_distribution<int> bit01(0, 1);

        long long global_total_bits = 0;
        long long global_correct_bits = 0;
        std::vector<double> t_each_ms;

        if (show_epoch_detail)
        {
            std::cout << "Bitwise XOR (A ^ B) Bit width = " << W
                      << "  (epochs=" << epochs << ")\n";
        }

        for (int ep = 0; ep < epochs; ++ep)
        {
            if (show_epoch_detail)
            {
                std::cout << "\n----- EPOCH " << (ep + 1) << "/" << epochs << " -----\n";
            }

            
            std::vector<uint8_t> Apl(W), Bpl(W);
            for (int i = 0; i < W; ++i)
            {
                Apl[i] = static_cast<uint8_t>(bit01(eng));
                Bpl[i] = static_cast<uint8_t>(bit01(eng));
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

            
            std::vector<TLWE<lvl2param>> xor_cts;
            double t_ms = 0.0;
            {
                AutoTimer timer(&t_ms);
                xor_cts = HomXOR(A, B, W, ek);
            }
            t_each_ms.push_back(t_ms);

            
            std::vector<uint8_t> XORpl(W), refXOR(W);
            long long correct_bits = 0;
            for (int i = 0; i < W; ++i)
            {
                XORpl[i] = tlweSymDecrypt<lvl2param>(xor_cts[i], sk.key.get<lvl2param>());
                refXOR[i] = static_cast<uint8_t>(Apl[i] ^ Bpl[i]);
                if (XORpl[i] == refXOR[i])
                    ++correct_bits;
            }
            global_total_bits += W;
            global_correct_bits += correct_bits;
            const bool epoch_ok = (correct_bits == W);

            if (show_epoch_detail)
            {
                auto print_bits = [](const char *tag, const std::vector<uint8_t> &v)
                {
                    std::cout << tag;
                    for (uint8_t b : v)
                        std::cout << int(b) << ' ';
                    std::cout << '\n';
                };
                print_bits("A   (msb->lsb): ", Apl);
                print_bits("B   (msb->lsb): ", Bpl);
                std::cout << "-----------------------------\n";
                print_bits("A ^ B (cipher): ", XORpl);
                print_bits("A ^ B (plain) : ", refXOR);
                std::cout << "Correct this epoch: " << correct_bits << "/" << W
                          << (epoch_ok ? " ✅" : " ❌") << "\n";
                std::cout << "Circuit execution: " << t_ms << " ms\n";
            }
            else
            {
                if (!epoch_ok)
                    std::cout << "[Epoch " << (ep + 1) << "] XOR mismatch!\n";
            }

            
            deep_free(Apl);
            deep_free(Bpl);
            deep_free(A);
            deep_free(B);
            deep_free(xor_cts);
            deep_free(XORpl);
            deep_free(refXOR);
        }

        
        double acc = (global_total_bits > 0)
                         ? (100.0 * double(global_correct_bits) / double(global_total_bits))
                         : 0.0;
        std::cout << "\n==== SUMMARY (XOR) ====\n";
        std::cout << "Average time over " << epochs << " epoch(s): " << avg(t_each_ms) << " ms\n";
        std::cout << "Global accuracy: " << global_correct_bits << "/" << global_total_bits
                  << " (" << acc << "%) " << (global_correct_bits == global_total_bits ? "✅" : "❌")
                  << "\n";

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return -1;
    }
}