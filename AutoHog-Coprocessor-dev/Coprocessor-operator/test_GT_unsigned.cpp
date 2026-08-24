#include "Autohog_Cop.hpp"
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


static cpp_int bits_to_int_unsigned_local(const std::vector<uint8_t> &bits, int W) {
    cpp_int res = 0;
    for (int i = 0; i < W; ++i) {
        if (bits[i]) {
            res |= (cpp_int(1) << (W - 1 - i));
        }
    }
    return res;
}

template <class T>
static double avg(const std::vector<T> &v) {
    if (v.empty()) return 0.0;
    double s = std::accumulate(v.begin(), v.end(), 0.0);
    return s / static_cast<double>(v.size());
}

static void print_usage(const char *prog) {
    std::cerr << "Usage: " << prog << " <bit-width | comma-list | all> [epochs] "
              << "[--no-epoch-detail|--show-epochs]\n"
              << "  Examples:\n"
              << "    " << prog << " 32 100 --show-epochs\n"
              << "    " << prog << " 8,16,32,64 20 --no-epoch-detail\n";
}

static std::vector<int> parse_widths(const std::string &s) {
    if (s == "all" || s == "ALL") return {8, 16, 32, 64};
    std::vector<int> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        int w = std::atoi(tok.c_str());
        if (w > 0) out.push_back(w);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

int main(int argc, char *argv[]) {
    try {
        if (argc < 2) {
            print_usage(argv[0]);
            return 1;
        }

        std::vector<int> widths = parse_widths(argv[1]);
        int epochs = (argc >= 3 && argv[2][0] != '-') ? std::max(1, std::atoi(argv[2])) : 1;
        bool show_epoch_detail = true;
        for (int i = 2; i < argc; ++i) {
            if (std::string(argv[i]) == "--no-epoch-detail") show_epoch_detail = false;
        }

        long long global_total = 0;
        long long global_correct_ser = 0, global_correct_par = 0;

        for (int W_req : widths) {
            std::cout << "\n================= UNSIGNED BIT-WIDTH " << W_req
                      << " (epochs=" << epochs << ") =================\n";

            std::vector<double> t_ser_vec, t_par_vec;
            int current_mismatches_ser = 0, current_mismatches_par = 0;

            for (int ep = 0; ep < epochs; ++ep) {
                
                SecretKey sk;
                EvalKey ek;
                ek.emplacebkfft<lvl02param>(sk);
                ek.emplaceiksk<lvl20param>(sk);
                init_constants(sk);

                int W = W_req;
                std::default_random_engine eng(std::random_device{}());
                std::uniform_int_distribution<int> bit01(0, 1);

                std::vector<uint8_t> Apl(W), Bpl(W);
                for (int i = 0; i < W; ++i) {
                    Apl[i] = bit01(eng);
                    Bpl[i] = bit01(eng);
                }

                
                std::vector<TLWE<lvl2param>> A(W), B(W);
                for (int i = 0; i < W; ++i) {
                    A[i] = tlweSymEncrypt<lvl2param>(Apl[i] ? lvl2param::μ : -lvl2param::μ, lvl2param::α, sk.key.get<lvl2param>());
                    B[i] = tlweSymEncrypt<lvl2param>(Bpl[i] ? lvl2param::μ : -lvl2param::μ, lvl2param::α, sk.key.get<lvl2param>());
                }

                
                TLWE<lvl2param> gt_ser;
                double t_ser = 0;
                {
                    AutoTimer timer(&t_ser);
                    gt_ser = HomGT_unsigned(A, B, W, ek);
                }
                t_ser_vec.push_back(t_ser);
                uint8_t res_ser = tlweSymDecrypt<lvl2param>(gt_ser, sk.key.get<lvl2param>()) ? 1 : 0;

                
                TLWE<lvl2param> gt_par;
                double t_par = 0;
                {
                    AutoTimer timer(&t_par);
                    gt_par = HomGT_unsigned_pal(A, B, W, ek, -1);
                }
                t_par_vec.push_back(t_par);
                uint8_t res_par = tlweSymDecrypt<lvl2param>(gt_par, sk.key.get<lvl2param>()) ? 1 : 0;

                
                cpp_int A_u = bits_to_int_unsigned_local(Apl, W);
                cpp_int B_u = bits_to_int_unsigned_local(Bpl, W);
                uint8_t ref_gt = (A_u > B_u) ? 1 : 0;

                if (res_ser == ref_gt) global_correct_ser++; else current_mismatches_ser++;
                if (res_par == ref_gt) global_correct_par++; else current_mismatches_par++;
                global_total++;

                if (show_epoch_detail) {
                    std::cout << "Epoch " << ep+1 << " | A: " << A_u << " B: " << B_u 
                              << " | GT_u: " << (int)res_ser << " (Ref: " << (int)ref_gt << ")"
                              << (res_ser == ref_gt ? " ✅" : " ❌") << " | " << t_ser << "ms" << std::endl;
                }
            }

            
            double avg_ser = avg(t_ser_vec);
            double avg_par = avg(t_par_vec);
            std::cout << std::fixed << std::setprecision(2);
            std::cout << "\n[Summary W=" << W_req << "]\n";
            std::cout << "  Avg Time Ser: " << avg_ser << " ms\n";
            std::cout << "  Avg Time Par: " << avg_par << " ms\n";
            std::cout << "  Speedup:      " << (avg_par > 0 ? avg_ser / avg_par : 0) << "x\n";
            std::cout << "  Acc Ser:      " << (100.0 * (epochs - current_mismatches_ser) / epochs) << "%\n";
            std::cout << "  Acc Par:      " << (100.0 * (epochs - current_mismatches_par) / epochs) << "%\n";
        }

        
        std::cout << "\n================= GLOBAL SUMMARY =================\n";
        std::cout << "Total Tests: " << global_total << "\n";
        std::cout << "Global Acc Ser: " << (100.0 * global_correct_ser / global_total) << "%\n";
        std::cout << "Global Acc Par: " << (100.0 * global_correct_par / global_total) << "%\n";

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }
    return 0;
}