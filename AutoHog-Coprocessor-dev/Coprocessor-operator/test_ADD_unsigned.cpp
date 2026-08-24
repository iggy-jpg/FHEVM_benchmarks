#include "Autohog_Cop.hpp"
#include "timer.h"
#include <tfhe++.hpp>

#include <random>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>

using namespace TFHEpp;
using namespace boost::multiprecision;


static cpp_int bits_to_int_unsigned(const std::vector<uint8_t> &bits, int W) {
    cpp_int res = 0;
    for (int i = 0; i < W; ++i) {
        if (bits[i]) {
            res |= (cpp_int(1) << (W - 1 - i));
        }
    }
    return res;
}

int main(int argc, char *argv[]) {
    
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <bit-width> [epochs]\n";
        return 1;
    }

    int W_arg = std::atoi(argv[1]);
    int epochs = (argc >= 3) ? std::max(1, std::atoi(argv[2])) : 1;

    
    std::cout << "--- Initializing Keys ---" << std::endl;
    std::cerr << "[TRACE] SecretKey begin" << std::endl;
    SecretKey sk;
    std::cerr << "[TRACE] SecretKey done; EvalKey begin" << std::endl;
    EvalKey ek;
    std::cerr << "[TRACE] EvalKey done; bkfft begin" << std::endl;
    ek.emplacebkfft<lvl02param>(sk);
    std::cerr << "[TRACE] bkfft done; iksk begin" << std::endl;
    ek.emplaceiksk<lvl20param>(sk);
    std::cerr << "[TRACE] iksk done; constants begin" << std::endl;
    init_constants(sk);
    std::cerr << "[TRACE] constants done" << std::endl;

    std::default_random_engine eng(std::random_device{}());
    std::uniform_int_distribution<int> bit01(0, 1);

    
    std::vector<double> t_add_ser, t_add_par, t_sub_ser, t_sub_par;
    
    int ok_add_ser = 0, ok_add_par = 0;
    int ok_sub_ser = 0, ok_sub_par = 0;

    cpp_int modW = cpp_int(1) << W_arg;

    for (int ep = 0; ep < epochs; ++ep) {
        std::cout << "\n--- Epoch " << ep + 1 << "/" << epochs << " ---" << std::endl;

        
        std::vector<uint8_t> Apl(W_arg), Bpl(W_arg);
        for (int i = 0; i < W_arg; ++i) {
            Apl[i] = bit01(eng);
            Bpl[i] = bit01(eng);
        }

        
        std::vector<TLWE<lvl2param>> A(W_arg), B(W_arg);
        for (int i = 0; i < W_arg; ++i) {
            A[i] = bootsSymEncrypt<lvl2param>(std::vector<uint8_t>{Apl[i]}, sk)[0];
            B[i] = bootsSymEncrypt<lvl2param>(std::vector<uint8_t>{Bpl[i]}, sk)[0];
        }

        
        cpp_int A_val = bits_to_int_unsigned(Apl, W_arg);
        cpp_int B_val = bits_to_int_unsigned(Bpl, W_arg);
        
        cpp_int expected_sum = (A_val + B_val) % modW;
        uint8_t expected_carry = (A_val + B_val >= modW) ? 1 : 0;
        
        cpp_int expected_diff;
        if (A_val >= B_val) {
            expected_diff = A_val - B_val;
        } else {
            expected_diff = modW + A_val - B_val;
        }
        uint8_t expected_borrow = (A_val < B_val) ? 1 : 0;

        // ==========================================
        
        // ==========================================
        HomAddResult res_add_s, res_add_p;
        double d_add_s = 0, d_add_p = 0;

        { 
            AutoTimer timer(&d_add_s);
            res_add_s = HomAdd_unsigned(A, B, W_arg, ek, true);
        }
        t_add_ser.push_back(d_add_s);

        { 
            AutoTimer timer(&d_add_p);
            res_add_p = HomAdd_unsigned_pal(A, B, W_arg, ek, true);
        }
        t_add_par.push_back(d_add_p);

        
        auto verify_add = [&](const HomAddResult& res) {
            std::vector<uint8_t> dec_bits(W_arg);
            for(int i=0; i<W_arg; ++i) dec_bits[i] = bootsSymDecrypt<lvl2param>({res.SUM[i]}, sk)[0];
            uint8_t carry = bootsSymDecrypt<lvl2param>({*res.ovf_ct}, sk)[0];
            return (bits_to_int_unsigned(dec_bits, W_arg) == expected_sum && carry == expected_carry);
        };

        if(verify_add(res_add_s)) ok_add_ser++;
        if(verify_add(res_add_p)) ok_add_par++;

        // ==========================================
        
        // ==========================================
        HomAddResult res_sub_s, res_sub_p;
        double d_sub_s = 0, d_sub_p = 0;

        { 
            AutoTimer timer(&d_sub_s);
            res_sub_s = HomSub_unsigned(A, B, W_arg, ek, true);
        }
        t_sub_ser.push_back(d_sub_s);

        { 
            AutoTimer timer(&d_sub_p);
            res_sub_p = HomSub_unsigned_pal(A, B, W_arg, ek, true);
        }
        t_sub_par.push_back(d_sub_p);

        
        auto verify_sub = [&](const HomAddResult& res) {
            std::vector<uint8_t> dec_bits(W_arg);
            for(int i=0; i<W_arg; ++i) dec_bits[i] = bootsSymDecrypt<lvl2param>({res.SUM[i]}, sk)[0];
            uint8_t borrow = bootsSymDecrypt<lvl2param>({*res.ovf_ct}, sk)[0];
            return (bits_to_int_unsigned(dec_bits, W_arg) == expected_diff && borrow == expected_borrow);
        };

        if(verify_sub(res_sub_s)) ok_sub_ser++;
        if(verify_sub(res_sub_p)) ok_sub_par++;

        std::cout << "[ADD] Serial: " << (verify_add(res_add_s)?"✅":"❌") << " | Parallel: " << (verify_add(res_add_p)?"✅":"❌") << std::endl;
        std::cout << "[SUB] Serial: " << (verify_sub(res_sub_s)?"✅":"❌") << " | Parallel: " << (verify_sub(res_sub_p)?"✅":"❌") << std::endl;
    }

    
    auto get_avg = [](const std::vector<double>& v) {
        return v.empty() ? 0.0 : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    };

    std::cout << "\n" << std::string(55, '=') << "\n";
    std::cout << "              FINAL PERFORMANCE SUMMARY\n";
    std::cout << std::string(55, '=') << "\n";
    std::cout << "Operation | Accuracy (Ser/Par) | Avg Time (Ser/Par)\n";
    std::cout << "----------|--------------------|--------------------\n";
    printf("ADD       | %d/%d (%.1f%%)      | %.2f / %.2f ms\n", 
           ok_add_ser, ok_add_par, (double)ok_add_par/epochs*100, get_avg(t_add_ser), get_avg(t_add_par));
    printf("SUB       | %d/%d (%.1f%%)      | %.2f / %.2f ms\n", 
           ok_sub_ser, ok_sub_par, (double)ok_sub_par/epochs*100, get_avg(t_sub_ser), get_avg(t_sub_par));
    std::cout << std::string(55, '=') << std::endl;

    return 0;
}
