#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <numeric>

#include <tfhe++.hpp>
#include "HomPack.h"

using namespace TFHEpp;


template <class T>
static double avg(const std::vector<T>& v) {
    if (v.empty()) return 0.0;
    double s = std::accumulate(v.begin(), v.end(), 0.0);
    return s / static_cast<double>(v.size());
}

static void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " <bit-width | comma-list | all> [num-test] [epochs] "
              << "[--no-epoch-detail|--show-epochs]\n"
              << "  Examples:\n"
              << "    " << prog << " 16 200 10 --show-epochs\n"
              << "    " << prog << " 8,16,32 100 --no-epoch-detail\n"
              << "    " << prog << " all 50\n";
}

static std::vector<int> parse_widths(const std::string& s) {
    if (s == "all" || s == "ALL") return {8,16,32,64,128};
    std::vector<int> out;
    if (s.find(',') == std::string::npos) {
        int w = std::atoi(s.c_str());
        if (w > 0) out.push_back(w);
        return out;
    }
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (tok.empty()) continue;
        int w = std::atoi(tok.c_str());
        if (w > 0) out.push_back(w);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) { print_usage(argv[0]); return 1; }

        
        std::vector<int> widths = parse_widths(argv[1]);
        if (widths.empty()) {
            std::cerr << "Error: invalid <bit-width | comma-list | all>.\n";
            return 1;
        }
        int num_test = 100;
        if (argc >= 3 && argv[2][0] != '-') {
            num_test = std::max(1, std::atoi(argv[2]));
        }
        int epochs = 1;
        if (argc >= 4 && argv[3][0] != '-') {
            epochs = std::max(1, std::atoi(argv[3]));
        }

        bool show_epoch_detail = true;
        for (int i = 2; i < argc; ++i) {
            std::string s = argv[i];
            if (s == "--no-epoch-detail") show_epoch_detail = false;
            else if (s == "--show-epochs") show_epoch_detail = true;
        }

        std::random_device rd;
        std::mt19937       eng(rd());
        std::uniform_int_distribution<uint32_t> bernoulli(0, 1);

        
        for (int W : widths) {
            std::cout << "\n================= BIT-WIDTH " << W
                      << " (num_test=" << num_test
                      << ", epochs=" << epochs << ") =================\n";

            std::vector<double> t_ms_all;
            int mismatches = 0;

            for (int ep = 0; ep < epochs; ++ep) {
                if (show_epoch_detail) {
                    std::cout << "\n----- EPOCH " << (ep + 1)
                              << "/" << epochs << " -----\n";
                }

                
                SecretKey sk;
                AnnihilateKey<lvl1param> ahk;
                annihilatekeygen<lvl1param>(ahk, sk);

                
                const uint numtlwe = W;
                std::vector<std::vector<uint8_t>> pin(num_test, std::vector<uint8_t>(numtlwe));
                for (auto& vec : pin) {
                    for (auto& b : vec) b = static_cast<uint8_t>(bernoulli(eng));
                }

                std::vector<std::vector<TLWE<lvl1param>>> ct(num_test);
                for (int i = 0; i < num_test; ++i) {
                    ct[i] = bootsSymEncrypt<lvl1param>(pin[i], sk);
                }

                
                std::vector<TRLWE<lvl1param>> packed(num_test);
                auto t0 = std::chrono::high_resolution_clock::now();
                for (int i = 0; i < num_test; ++i) {
                    HomPack::TLWE2TRLWEPacking<lvl1param>(packed[i], ct[i], ahk);
                }
                auto t1 = std::chrono::high_resolution_clock::now();

                const double ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
                                  / static_cast<double>(num_test);
                t_ms_all.push_back(ms);

                
                for (int i = 0; i < num_test; ++i) {
                    const auto dec = trlweSymDecrypt<lvl1param>(packed[i], sk.key.get<lvl1param>());
                    for (int j = 0; j < numtlwe; ++j) {
                        if (static_cast<uint8_t>(dec[j]) != (pin[i][j] > 0)) {
                            mismatches++;
                        }
                    }
                }

                if (show_epoch_detail) {
                    std::cout << "PACK Passed (Epoch " << (ep+1) << ")\n";
                    std::cout << "Avg: " << ms << " ms\n";
                }
            } // epochs

            if (epochs > 1) {
                std::cout << "\n---- SUMMARY (Bit-width " << W << ") ----\n";
                std::cout << "Average packing time: " << avg(t_ms_all) << " ms\n";
                std::cout << "Mismatches: " << mismatches << "\n";
            }
        }

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return -1;
    }
}