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


static cpp_int trunc_signed(cpp_int x, int W)
{
    cpp_int mod = cpp_int(1);
    mod <<= W;              // 2^W
    cpp_int mask = mod - 1; 
    cpp_int y = x & mask;   
    
    if (((y >> (W - 1)) & 1) != 0)
        y -= mod;
    return y;
}

static void print_usage(const char *prog)
{
    std::cerr << "Usage: " << prog
              << " <bit-width | comma-list | all> [epochs] "
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

        
        SecretKey sk;
        EvalKey ek;
        ek.emplacebkfft<lvl02param>(sk);
        ek.emplaceiksk<lvl20param>(sk);
        init_constants(sk); 

        
        for (size_t wi = 0; wi < widths.size(); ++wi)
        {
            int W_req = widths[wi];
            std::cout << "\n================= BIT-WIDTH " << W_req
                      << " (epochs=" << epochs << ") =================\n";

            
            std::vector<double> t_ms_ser_vec, t_ms_par_vec;
            long long total_bits = 0;
            long long correct_bits_ser = 0;
            long long correct_bits_par = 0;

            int result_ok_ser = 0; 
            int result_ok_par = 0; 
            int ovf_ok_ser = 0;    
            int ovf_ok_par = 0;    
            int ser_par_equal = 0; 

            for (int ep = 0; ep < epochs; ++ep)
            {
                if (show_epoch_detail)
                {
                    std::cout << "\n----- EPOCH " << (ep + 1) << "/" << epochs << " -----\n";
                }

                
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

                // -------------------- Plaintext Generation --------------------
                std::vector<uint8_t> Apl(W), Bpl(W);
                std::default_random_engine eng(std::random_device{}());
                std::uniform_int_distribution<int> bit01(0, 1);
                for (int i = 0; i < W; ++i)
                {
                    Apl[i] = bit01(eng);
                    Bpl[i] = bit01(eng);
                }

                // -------------------- Encryption（0→-μ，1→+μ） --------------------
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

                
                cpp_int A_s = bits_to_int_signed(Apl, W);
                cpp_int B_s = bits_to_int_signed(Bpl, W);
                cpp_int plain_sub = A_s - B_s;
                cpp_int plain_sub_mod = trunc_signed(plain_sub, W);
                const uint8_t expected_ovf = (plain_sub != plain_sub_mod) ? 1 : 0;

                
                HomAddResult res_ser;
                double t_ser_ms = 0.0;
                {
                    AutoTimer timer(&t_ser_ms);
                    res_ser = HomSub(A, B, W, ek, /*compute_overflow=*/true);
                }
                t_ms_ser_vec.push_back(t_ser_ms);

                
                HomAddResult res_par;
                double t_par_ms = 0.0;
                {
                    AutoTimer timer(&t_par_ms);
                    res_par = HomSub_pal(A, B, W, ek, /*compute_overflow=*/true);
                }
                t_ms_par_vec.push_back(t_par_ms);

                
                std::vector<uint8_t> Diff_ser(W), Diff_par(W);
                long long correct_epoch_bits_ser = 0;
                long long correct_epoch_bits_par = 0;

                for (int i = 0; i < W; ++i)
                {
                    Diff_ser[i] = tlweSymDecrypt<lvl2param>(res_ser.SUM[i], sk.key.get<lvl2param>());
                    Diff_par[i] = tlweSymDecrypt<lvl2param>(res_par.SUM[i], sk.key.get<lvl2param>());
                }

                cpp_int S_ser = bits_to_int_signed(Diff_ser, W);
                cpp_int S_par = bits_to_int_signed(Diff_par, W);

                
                for (int i = 0; i < W; ++i)
                {
                    uint8_t refbit = static_cast<uint8_t>((plain_sub_mod >> (W - 1 - i)) & 1);
                    if (Diff_ser[i] == refbit)
                        ++correct_epoch_bits_ser;
                    if (Diff_par[i] == refbit)
                        ++correct_epoch_bits_par;
                }

                total_bits += W;
                correct_bits_ser += correct_epoch_bits_ser;
                correct_bits_par += correct_epoch_bits_par;

                
                const bool ok_ser = (S_ser == plain_sub_mod);
                const bool ok_par = (S_par == plain_sub_mod);
                if (ok_ser)
                    ++result_ok_ser;
                if (ok_par)
                    ++result_ok_par;

                
                uint8_t ovf_ser = 0, ovf_par = 0;
                if (res_ser.ovf_ct)
                    ovf_ser = tlweSymDecrypt<lvl2param>(*res_ser.ovf_ct, sk.key.get<lvl2param>());
                if (res_par.ovf_ct)
                    ovf_par = tlweSymDecrypt<lvl2param>(*res_par.ovf_ct, sk.key.get<lvl2param>());
                if (ovf_ser == expected_ovf)
                    ++ovf_ok_ser;
                if (ovf_par == expected_ovf)
                    ++ovf_ok_par;

                
                bool ser_eq_par = true;
                for (int i = 0; i < W; ++i)
                    if (Diff_ser[i] != Diff_par[i])
                    {
                        ser_eq_par = false;
                        break;
                    }
                if (ser_eq_par && ovf_ser == ovf_par)
                    ++ser_par_equal;

                
                if (show_epoch_detail)
                {
                    print_bits("A   (msb->lsb): ", Apl);
                    print_bits("B   (msb->lsb): ", Bpl);
                    std::cout << "-----------------------------\n";
                    print_bits("Diff(serial):    ", Diff_ser);
                    print_bits("Diff(parallel):  ", Diff_par);
                    std::cout << "A (signed)       = " << A_s << "\n";
                    std::cout << "B (signed)       = " << B_s << "\n";
                    std::cout << "A - B (signed)   = " << plain_sub << "  (trunc@" << W
                              << " -> " << plain_sub_mod << ")\n";
                    std::cout << "OVF expected     = " << int(expected_ovf) << "\n";
                    std::cout << "OVF(serial/par)  = " << int(ovf_ser) << " / " << int(ovf_par) << "\n";
                    std::cout << "SER ok? " << (ok_ser ? "✅" : "❌")
                              << "   PAR ok? " << (ok_par ? "✅" : "❌")
                              << "   SER==PAR? " << (ser_eq_par ? "✅" : "❌") << "\n";
                    std::cout << "Time (ms): SER=" << t_ser_ms << "  PAR=" << t_par_ms << "\n";
                }

                
                deep_free(Apl);
                deep_free(Bpl);
                deep_free(A);
                deep_free(B);
                deep_free(res_ser.SUM);
                if (res_ser.ovf_ct)
                    res_ser.ovf_ct.reset();
                deep_free(res_par.SUM);
                if (res_par.ovf_ct)
                    res_par.ovf_ct.reset();
                deep_free(Diff_ser);
                deep_free(Diff_par);
            } // end epochs

            
            const double acc_ser = (total_bits > 0)
                                       ? (100.0 * double(correct_bits_ser) / double(total_bits))
                                       : 0.0;
            const double acc_par = (total_bits > 0)
                                       ? (100.0 * double(correct_bits_par) / double(total_bits))
                                       : 0.0;

            std::cout << "\n---- SUMMARY (Bit-width " << W_req
                      << ", average over " << epochs << " epochs) ----\n";
            std::cout << "SER avg time: " << avg(t_ms_ser_vec) << " ms\n";
            std::cout << "PAR avg time: " << avg(t_ms_par_vec) << " ms\n";

            std::cout << "Result correctness (SER): " << result_ok_ser << "/" << epochs
                      << (result_ok_ser == epochs ? " ✅" : " ❌") << "\n";
            std::cout << "Result correctness (PAR): " << result_ok_par << "/" << epochs
                      << (result_ok_par == epochs ? " ✅" : " ❌") << "\n";

            std::cout << "Overflow-bit correctness (SER): " << ovf_ok_ser << "/" << epochs
                      << (ovf_ok_ser == epochs ? " ✅" : " ❌") << "\n";
            std::cout << "Overflow-bit correctness (PAR): " << ovf_ok_par << "/" << epochs
                      << (ovf_ok_par == epochs ? " ✅" : " ❌") << "\n";

            std::cout << "SER vs PAR identical epochs: " << ser_par_equal << "/" << epochs
                      << (ser_par_equal == epochs ? " ✅" : " ❌") << "\n";

            std::cout << "Global bit-accuracy (SER): " << correct_bits_ser << "/" << total_bits
                      << " (" << acc_ser << "%) " << (correct_bits_ser == total_bits ? "✅" : "❌")
                      << "\n";
            std::cout << "Global bit-accuracy (PAR): " << correct_bits_par << "/" << total_bits
                      << " (" << acc_par << "%) " << (correct_bits_par == total_bits ? "✅" : "❌")
                      << "\n";
        } // end widths

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return -1;
    }
}