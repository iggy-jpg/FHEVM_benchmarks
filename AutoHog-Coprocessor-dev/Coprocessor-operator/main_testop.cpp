#include "Autohog_Cop.hpp" 
#include "fhe_key.h"
#include "timer.h"
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace TFHEpp;

void print_usage(const char* prog) {
    cout << "Usage: " << prog << " <op> <bit_width> <mode> [epochs] [detail]\n"
         << "\nOptions:\n"
         << "  op         : FheAdd, FheSub, FheBitAnd, FheBitOr, FheBitXor, FheBitNot, FheEq, FheNe, FheGe, FheGt, FheLe, FheLt, FheIfThenElse, TrivialEncrypt\n"
         << "  bit_width  : number of bits for the operands (must be a multiple of 8, e.g., 8, 16, 32, 64)\n"
         << "  mode       : 'std' (Standard Encryption) or 'triv' (Trivial Encryption)\n"
         << "  epochs     : (Optional) iterations for benchmarking, default 1\n"
         << "  detail     : (Optional) 1 to show per-epoch details, 0 to hide, default 1\n"
         << "\nExample:\n"
         << "  " << prog << " FheAdd 32 triv 5 0\n";
}

// ====================================================================

// ====================================================================
static boost::multiprecision::cpp_int bytes_to_cpp_int(const Bytes& data) {
    boost::multiprecision::cpp_int res = 0;
    for (uint8_t b : data) {
        res = (res << 8) | b;
    }
    return res;
}

static bool are_ciphertexts_identical(const tlwe_ptr& a, const tlwe_ptr& b) {
    if (a->data.size() != b->data.size()) return false;
    return std::memcmp(a->data.data(), b->data.data(), sizeof(a->data)) == 0;
}

static bool are_ciphertext_vecs_identical(const tlwe_vec_ptr& a, const tlwe_vec_ptr& b) {
    if (a->W != b->W) return false;
    if (a->data.size() != b->data.size()) return false;
    for(size_t i = 0; i < a->data.size(); ++i) {
        if (std::memcmp(a->data[i].data(), b->data[i].data(), sizeof(a->data[i])) != 0) {
            return false;
        }
    }
    return true;
}
// ====================================================================

// ====================================================================
static boost::multiprecision::cpp_int bytes_to_signed_cpp_int(const Bytes& data) {
    int W = data.size() * 8;
    boost::multiprecision::cpp_int res = 0;
    
    
    for (uint8_t b : data) {
        res = (res << 8) | b;
    }
    
    
    if (data.size() > 0 && (data[0] & 0x80)) {
        
        boost::multiprecision::cpp_int mod = boost::multiprecision::cpp_int(1) << W;
        res -= mod;
    }
    
    return res;
}


void test_FheAdd(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail);
void test_FheSub(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail);


void test_FheBitAnd(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail);
void test_FheBitOr(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail);
void test_FheBitXor(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail);
void test_FheBitNot(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail); 



void test_FheEq(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail); // A == B
void test_FheNe(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail); // A != B
void test_FheGt(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail); // A > B
void test_FheGe(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail); // A >= B
void test_FheLt(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail); // A < B
void test_FheLe(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail); // A <= B



void test_FheIfThenElse(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const tlwe_ptr& ctCond, const Bytes& pA, const Bytes& pB, uint8_t pCond, int epochs, int detail);

void test_TrivialEncrypt(FheKey* key, const Bytes& pA, int epochs, int detail);


// ====================================================================

// ====================================================================



int main(int argc, char* argv[]) {
    if (argc < 4) {
        print_usage(argv[0]);
        return 1;
    }

    string op_str    = argv[1];
    int bit_width    = stoi(argv[2]);
    string mode      = argv[3];
    int epochs       = (argc >= 5) ? stoi(argv[4]) : 1;
    int detail       = (argc >= 6) ? stoi(argv[5]) : 1; 

    if (bit_width <= 0 || bit_width % 8 != 0) {
        cerr << "Error: bit_width must be > 0 and a multiple of 8 (e.g., 8, 16, 32, 64)" << endl; 
        return 1;
    }
    
    
    int byte_count = bit_width / 8;

    if (mode != "std" && mode != "triv") {
        cerr << "Error: mode must be 'std' or 'triv'" << endl; return 1;
    }

#ifdef _OPENMP
    int hw_threads = max(1u, thread::hardware_concurrency());
    omp_set_num_threads(hw_threads);
#endif

    
    cout << "--- Initializing FheKey ---" << endl;
    auto key = create_fhe_key();

    
    default_random_engine eng(random_device{}());
    uniform_int_distribution<int> dist(0, 255);
    uniform_int_distribution<int> bit_dist(0, 1);
    
    Bytes pA(byte_count), pB(byte_count);
    for(int i = 0; i < byte_count; ++i) {
        pA[i] = dist(eng);
        pB[i] = dist(eng);
    }
    uint8_t pCond = bit_dist(eng); 

    
    tlwe_vec_ptr ctA, ctB;
    tlwe_ptr ctCond;

    if (op_str != "TrivialEncrypt") {
        cout << "--- Encrypting Inputs (Mode: " << mode << ", Bits: " << bit_width << ") ---" << endl;
        if (mode == "triv") {
            ctA = decompress_ciphertexts(key->trivialEncryptBatch(pA));
            ctB = decompress_ciphertexts(key->trivialEncryptBatch(pB));
            ctCond = to_tlwe_ptr(decompress_ciphertexts(key->trivialEncrypt(pCond)));
        } else {
            ctA = decompress_ciphertexts(key->encryptBatch(pA));
            ctB = decompress_ciphertexts(key->encryptBatch(pB));
            ctCond = to_tlwe_ptr(decompress_ciphertexts(key->encrypt(pCond)));
        }
    }

    
    cout << "\n=======================================================\n";
    cout << " TARGET OP : " << op_str << "\n";
    cout << " MODE      : " << mode << "\n";
    cout << " EPOCHS    : " << epochs << "\n";
    cout << " DETAIL    : " << (detail ? "ON" : "OFF") << "\n";
    cout << "=======================================================\n";

    if      (op_str == "FheAdd")         test_FheAdd(key.get(), ctA, ctB, pA, pB, epochs, detail);
    else if (op_str == "FheSub")         test_FheSub(key.get(), ctA, ctB, pA, pB, epochs, detail);
    else if (op_str == "FheBitAnd")      test_FheBitAnd(key.get(), ctA, ctB, pA, pB, epochs, detail);
    else if (op_str == "FheBitOr")       test_FheBitOr(key.get(), ctA, ctB, pA, pB, epochs, detail);
    else if (op_str == "FheBitXor")      test_FheBitXor(key.get(), ctA, ctB, pA, pB, epochs, detail);
    else if (op_str == "FheBitNot")      test_FheBitNot(key.get(), ctA, ctB, pA, pB, epochs, detail);
    
    else if (op_str == "FheEq")          test_FheEq(key.get(), ctA, ctB, pA, pB, epochs, detail);
    else if (op_str == "FheNe")          test_FheNe(key.get(), ctA, ctB, pA, pB, epochs, detail);
    else if (op_str == "FheGt")          test_FheGt(key.get(), ctA, ctB, pA, pB, epochs, detail);
    else if (op_str == "FheGe")          test_FheGe(key.get(), ctA, ctB, pA, pB, epochs, detail);
    else if (op_str == "FheLt")          test_FheLt(key.get(), ctA, ctB, pA, pB, epochs, detail);
    else if (op_str == "FheLe")          test_FheLe(key.get(), ctA, ctB, pA, pB, epochs, detail);
    
    else if (op_str == "FheIfThenElse")  test_FheIfThenElse(key.get(), ctA, ctB, ctCond, pA, pB, pCond, epochs, detail);
    else if (op_str == "TrivialEncrypt") test_TrivialEncrypt(key.get(), pA, epochs, detail);
    else {
        cerr << "Unknown operation: " << op_str << endl;
        cerr << "Please use the full name (e.g., FheAdd, FheBitAnd)." << endl;
        return 1;
    }

    return 0;
}

// ====================================================================

// ====================================================================
void test_FheAdd(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail) {
    std::cout << "\n[*] Running FheAdd test (Serial vs Parallel)..." << std::endl;

    int W = pA.size() * 8; 
    int hw_threads = std::max(1u, std::thread::hardware_concurrency());

    std::vector<double> t_ser, t_par;
    int ok_serial_vs_plain = 0;
    int ok_parallel_vs_serial = 0;

    
    boost::multiprecision::cpp_int a_val = bytes_to_cpp_int(pA);
    boost::multiprecision::cpp_int b_val = bytes_to_cpp_int(pB);
    boost::multiprecision::cpp_int plain_add = a_val + b_val;
    
    
    boost::multiprecision::cpp_int mask = (boost::multiprecision::cpp_int(1) << W) - 1;
    boost::multiprecision::cpp_int expected_sum = plain_add & mask;

    for (int ep = 0; ep < epochs; ++ep) {
        double ms_ser = 0.0, ms_par = 0.0;

        
        tlwe_vec_ptr res_ser;
        {
            AutoTimer t(&ms_ser);
            res_ser = key->add(ctA, ctB, false); 
        }

        
        tlwe_vec_ptr res_par;
        {
            AutoTimer t(&ms_par);
            res_par = key->add_pal(ctA, ctB, false, hw_threads);
        }

        t_ser.push_back(ms_ser);
        t_par.push_back(ms_par);

        
        Bytes dec_ser = key->decryptBatch(res_ser);
        Bytes dec_par = key->decryptBatch(res_par);

        boost::multiprecision::cpp_int sum_ser = bytes_to_cpp_int(dec_ser);

        
        bool pass_serial_vs_plain = (sum_ser == expected_sum);
        bool pass_parallel_vs_serial = (dec_ser == dec_par); 

        if (pass_serial_vs_plain) ok_serial_vs_plain++;
        if (pass_parallel_vs_serial) ok_parallel_vs_serial++;

        
        if (detail) {
            std::cout << "  Epoch [" << ep + 1 << "/" << epochs << "] "
                      << "| Serial: " << ms_ser << " ms "
                      << "| Parallel (" << hw_threads << "T): " << ms_par << " ms\n"
                      << "    -> Correctness: " << (pass_serial_vs_plain ? "PASS ✅" : "FAIL ❌") 
                      << " | Par==Ser: " << (pass_parallel_vs_serial ? "PASS ✅" : "FAIL ❌") << std::endl;
        }
    }

    
    double avg_ser = std::accumulate(t_ser.begin(), t_ser.end(), 0.0) / epochs;
    double avg_par = std::accumulate(t_par.begin(), t_par.end(), 0.0) / epochs;

    auto percent = [](int num, int den) -> double {
        return (den > 0) ? (100.0 * num / den) : 0.0;
    };

    std::cout << "\n================= SUMMARY (FheAdd) =================\n";
    std::cout << "Data Size     : " << pA.size() << " Bytes (" << W << " bits)\n";
    std::cout << "HomAdd (ser)  : " << avg_ser << " ms (avg)\n";
    std::cout << "HomAdd (pal)  : " << avg_par << " ms (avg)\n";
    if (avg_par > 0) {
        std::cout << "Speedup       : x" << (avg_ser / avg_par) << "\n";
    }
    std::cout << "----------------------------------------------------\n";
    std::cout << "Serial vs Plain 2^W : " << ok_serial_vs_plain << "/" << epochs 
              << " (" << percent(ok_serial_vs_plain, epochs) << "%)\n";
    std::cout << "Parallel vs Serial  : " << ok_parallel_vs_serial << "/" << epochs 
              << " (" << percent(ok_parallel_vs_serial, epochs) << "%)\n";
    std::cout << "====================================================\n";
}

// ====================================================================

// ====================================================================
void test_FheSub(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail) {
    std::cout << "\n[*] Running FheSub test (Serial vs Parallel)..." << std::endl;

    int W = pA.size() * 8; 
    int hw_threads = std::max(1u, std::thread::hardware_concurrency());

    std::vector<double> t_ser, t_par;
    int ok_serial_vs_plain = 0;
    int ok_parallel_vs_serial = 0;

    
    boost::multiprecision::cpp_int a_val = bytes_to_cpp_int(pA);
    boost::multiprecision::cpp_int b_val = bytes_to_cpp_int(pB);
    boost::multiprecision::cpp_int plain_sub = a_val - b_val;
    
    
    boost::multiprecision::cpp_int mask = (boost::multiprecision::cpp_int(1) << W) - 1;
    boost::multiprecision::cpp_int expected_sub = plain_sub & mask;

    for (int ep = 0; ep < epochs; ++ep) {
        double ms_ser = 0.0, ms_par = 0.0;

        
        tlwe_vec_ptr res_ser;
        {
            AutoTimer t(&ms_ser);
            
            res_ser = key->sub(ctA, ctB, false); 
        }

        
        tlwe_vec_ptr res_par;
        {
            AutoTimer t(&ms_par);
            res_par = key->sub_pal(ctA, ctB, false, hw_threads);
        }

        t_ser.push_back(ms_ser);
        t_par.push_back(ms_par);

        
        Bytes dec_ser = key->decryptBatch(res_ser);
        Bytes dec_par = key->decryptBatch(res_par);

        boost::multiprecision::cpp_int sub_ser = bytes_to_cpp_int(dec_ser);

        
        bool pass_serial_vs_plain = (sub_ser == expected_sub);
        bool pass_parallel_vs_serial = (dec_ser == dec_par); 

        if (pass_serial_vs_plain) ok_serial_vs_plain++;
        if (pass_parallel_vs_serial) ok_parallel_vs_serial++;

        
        if (detail) {
            std::cout << "  Epoch [" << ep + 1 << "/" << epochs << "] "
                      << "| Serial: " << ms_ser << " ms "
                      << "| Parallel (" << hw_threads << "T): " << ms_par << " ms\n"
                      << "    -> Correctness: " << (pass_serial_vs_plain ? "PASS ✅" : "FAIL ❌") 
                      << " | Par==Ser: " << (pass_parallel_vs_serial ? "PASS ✅" : "FAIL ❌") << std::endl;
        }
    }

    
    double avg_ser = std::accumulate(t_ser.begin(), t_ser.end(), 0.0) / epochs;
    double avg_par = std::accumulate(t_par.begin(), t_par.end(), 0.0) / epochs;

    auto percent = [](int num, int den) -> double {
        return (den > 0) ? (100.0 * num / den) : 0.0;
    };

    std::cout << "\n================= SUMMARY (FheSub) =================\n";
    std::cout << "Data Size     : " << pA.size() << " Bytes (" << W << " bits)\n";
    std::cout << "HomSub (ser)  : " << avg_ser << " ms (avg)\n";
    std::cout << "HomSub (pal)  : " << avg_par << " ms (avg)\n";
    if (avg_par > 0) {
        std::cout << "Speedup       : x" << (avg_ser / avg_par) << "\n";
    }
    std::cout << "----------------------------------------------------\n";
    std::cout << "Serial vs Plain 2^W : " << ok_serial_vs_plain << "/" << epochs 
              << " (" << percent(ok_serial_vs_plain, epochs) << "%)\n";
    std::cout << "Parallel vs Serial  : " << ok_parallel_vs_serial << "/" << epochs 
              << " (" << percent(ok_parallel_vs_serial, epochs) << "%)\n";
    std::cout << "====================================================\n";
}

// ====================================================================

// ====================================================================
void test_FheBitAnd(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail) {
    std::cout << "\n[*] Running FheBitAnd test (Serial vs Parallel)..." << std::endl;

    int W = pA.size() * 8; 
    int hw_threads = std::max(1u, std::thread::hardware_concurrency());

    std::vector<double> t_ser, t_par;
    int ok_serial_vs_plain = 0;
    int ok_parallel_vs_serial = 0;

    
    Bytes expected_and(pA.size());
    for (size_t i = 0; i < pA.size(); ++i) {
        expected_and[i] = pA[i] & pB[i];
    }

    for (int ep = 0; ep < epochs; ++ep) {
        double ms_ser = 0.0, ms_par = 0.0;

        
        tlwe_vec_ptr res_ser;
        {
            AutoTimer t(&ms_ser);
            res_ser = key->bit_and(ctA, ctB, 1); 
        }

        
        tlwe_vec_ptr res_par;
        {
            AutoTimer t(&ms_par);
            res_par = key->bit_and(ctA, ctB, hw_threads);
        }

        t_ser.push_back(ms_ser);
        t_par.push_back(ms_par);

        
        Bytes dec_ser = key->decryptBatch(res_ser);
        Bytes dec_par = key->decryptBatch(res_par);

        
        bool pass_serial_vs_plain = (dec_ser == expected_and);
        bool pass_parallel_vs_serial = (dec_ser == dec_par); 

        if (pass_serial_vs_plain) ok_serial_vs_plain++;
        if (pass_parallel_vs_serial) ok_parallel_vs_serial++;

        
        if (detail) {
            std::cout << "  Epoch [" << ep + 1 << "/" << epochs << "] "
                      << "| Serial: " << ms_ser << " ms "
                      << "| Parallel (" << hw_threads << "T): " << ms_par << " ms\n"
                      << "    -> Correctness: " << (pass_serial_vs_plain ? "PASS ✅" : "FAIL ❌") 
                      << " | Par==Ser: " << (pass_parallel_vs_serial ? "PASS ✅" : "FAIL ❌") << std::endl;
        }
    }

    
    double avg_ser = std::accumulate(t_ser.begin(), t_ser.end(), 0.0) / epochs;
    double avg_par = std::accumulate(t_par.begin(), t_par.end(), 0.0) / epochs;

    auto percent = [](int num, int den) -> double {
        return (den > 0) ? (100.0 * num / den) : 0.0;
    };

    std::cout << "\n================= SUMMARY (FheBitAnd) =================\n";
    std::cout << "Data Size     : " << pA.size() << " Bytes (" << W << " bits)\n";
    std::cout << "HomAND (ser)  : " << avg_ser << " ms (avg)\n";
    std::cout << "HomAND (pal)  : " << avg_par << " ms (avg)\n";
    if (avg_par > 0) {
        std::cout << "Speedup       : x" << (avg_ser / avg_par) << "\n";
    }
    std::cout << "----------------------------------------------------\n";
    std::cout << "Serial vs Plain AND : " << ok_serial_vs_plain << "/" << epochs 
              << " (" << percent(ok_serial_vs_plain, epochs) << "%)\n";
    std::cout << "Parallel vs Serial  : " << ok_parallel_vs_serial << "/" << epochs 
              << " (" << percent(ok_parallel_vs_serial, epochs) << "%)\n";
    std::cout << "====================================================\n";
}

// ====================================================================

// ====================================================================
void test_FheBitOr(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail) {
    std::cout << "\n[*] Running FheBitOr test (Serial vs Parallel)..." << std::endl;

    int W = pA.size() * 8; 
    int hw_threads = std::max(1u, std::thread::hardware_concurrency());

    std::vector<double> t_ser, t_par;
    int ok_serial_vs_plain = 0;
    int ok_parallel_vs_serial = 0;

    
    Bytes expected_or(pA.size());
    for (size_t i = 0; i < pA.size(); ++i) {
        expected_or[i] = pA[i] | pB[i];
    }

    for (int ep = 0; ep < epochs; ++ep) {
        double ms_ser = 0.0, ms_par = 0.0;

        
        tlwe_vec_ptr res_ser;
        {
            AutoTimer t(&ms_ser);
            res_ser = key->bit_or(ctA, ctB, 1); 
        }

        
        tlwe_vec_ptr res_par;
        {
            AutoTimer t(&ms_par);
            res_par = key->bit_or(ctA, ctB, hw_threads);
        }

        t_ser.push_back(ms_ser);
        t_par.push_back(ms_par);

        
        Bytes dec_ser = key->decryptBatch(res_ser);
        Bytes dec_par = key->decryptBatch(res_par);

        
        bool pass_serial_vs_plain = (dec_ser == expected_or);
        bool pass_parallel_vs_serial = (dec_ser == dec_par); 

        if (pass_serial_vs_plain) ok_serial_vs_plain++;
        if (pass_parallel_vs_serial) ok_parallel_vs_serial++;

        
        if (detail) {
            std::cout << "  Epoch [" << ep + 1 << "/" << epochs << "] "
                      << "| Serial: " << ms_ser << " ms "
                      << "| Parallel (" << hw_threads << "T): " << ms_par << " ms\n"
                      << "    -> Correctness: " << (pass_serial_vs_plain ? "PASS ✅" : "FAIL ❌") 
                      << " | Par==Ser: " << (pass_parallel_vs_serial ? "PASS ✅" : "FAIL ❌") << std::endl;
        }
    }

    
    double avg_ser = std::accumulate(t_ser.begin(), t_ser.end(), 0.0) / epochs;
    double avg_par = std::accumulate(t_par.begin(), t_par.end(), 0.0) / epochs;

    auto percent = [](int num, int den) -> double {
        return (den > 0) ? (100.0 * num / den) : 0.0;
    };

    std::cout << "\n================= SUMMARY (FheBitOr) =================\n";
    std::cout << "Data Size     : " << pA.size() << " Bytes (" << W << " bits)\n";
    std::cout << "HomOR (ser)   : " << avg_ser << " ms (avg)\n";
    std::cout << "HomOR (pal)   : " << avg_par << " ms (avg)\n";
    if (avg_par > 0) {
        std::cout << "Speedup       : x" << (avg_ser / avg_par) << "\n";
    }
    std::cout << "----------------------------------------------------\n";
    std::cout << "Serial vs Plain OR  : " << ok_serial_vs_plain << "/" << epochs 
              << " (" << percent(ok_serial_vs_plain, epochs) << "%)\n";
    std::cout << "Parallel vs Serial  : " << ok_parallel_vs_serial << "/" << epochs 
              << " (" << percent(ok_parallel_vs_serial, epochs) << "%)\n";
    std::cout << "====================================================\n";
}

// ====================================================================

// ====================================================================
void test_FheBitXor(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail) {
    std::cout << "\n[*] Running FheBitXor test (Serial vs Parallel)..." << std::endl;

    int W = pA.size() * 8; 
    int hw_threads = std::max(1u, std::thread::hardware_concurrency());

    std::vector<double> t_ser, t_par;
    int ok_serial_vs_plain = 0;
    int ok_parallel_vs_serial = 0;

    
    Bytes expected_xor(pA.size());
    for (size_t i = 0; i < pA.size(); ++i) {
        expected_xor[i] = pA[i] ^ pB[i];
    }

    for (int ep = 0; ep < epochs; ++ep) {
        double ms_ser = 0.0, ms_par = 0.0;

        
        tlwe_vec_ptr res_ser;
        {
            AutoTimer t(&ms_ser);
            res_ser = key->bit_xor(ctA, ctB, 1); 
        }

        
        tlwe_vec_ptr res_par;
        {
            AutoTimer t(&ms_par);
            res_par = key->bit_xor(ctA, ctB, hw_threads);
        }

        t_ser.push_back(ms_ser);
        t_par.push_back(ms_par);

        
        Bytes dec_ser = key->decryptBatch(res_ser);
        Bytes dec_par = key->decryptBatch(res_par);

        
        bool pass_serial_vs_plain = (dec_ser == expected_xor);
        bool pass_parallel_vs_serial = (dec_ser == dec_par); 

        if (pass_serial_vs_plain) ok_serial_vs_plain++;
        if (pass_parallel_vs_serial) ok_parallel_vs_serial++;

        
        if (detail) {
            std::cout << "  Epoch [" << ep + 1 << "/" << epochs << "] "
                      << "| Serial: " << ms_ser << " ms "
                      << "| Parallel (" << hw_threads << "T): " << ms_par << " ms\n"
                      << "    -> Correctness: " << (pass_serial_vs_plain ? "PASS ✅" : "FAIL ❌") 
                      << " | Par==Ser: " << (pass_parallel_vs_serial ? "PASS ✅" : "FAIL ❌") << std::endl;
        }
    }

    
    double avg_ser = std::accumulate(t_ser.begin(), t_ser.end(), 0.0) / epochs;
    double avg_par = std::accumulate(t_par.begin(), t_par.end(), 0.0) / epochs;

    auto percent = [](int num, int den) -> double {
        return (den > 0) ? (100.0 * num / den) : 0.0;
    };

    std::cout << "\n================= SUMMARY (FheBitXor) =================\n";
    std::cout << "Data Size     : " << pA.size() << " Bytes (" << W << " bits)\n";
    std::cout << "HomXOR (ser)  : " << avg_ser << " ms (avg)\n";
    std::cout << "HomXOR (pal)  : " << avg_par << " ms (avg)\n";
    if (avg_par > 0) {
        std::cout << "Speedup       : x" << (avg_ser / avg_par) << "\n";
    }
    std::cout << "----------------------------------------------------\n";
    std::cout << "Serial vs Plain XOR : " << ok_serial_vs_plain << "/" << epochs 
              << " (" << percent(ok_serial_vs_plain, epochs) << "%)\n";
    std::cout << "Parallel vs Serial  : " << ok_parallel_vs_serial << "/" << epochs 
              << " (" << percent(ok_parallel_vs_serial, epochs) << "%)\n";
    std::cout << "====================================================\n";
}

// ====================================================================

// ====================================================================
void test_FheBitNot(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail) {
    std::cout << "\n[*] Running FheBitNot test (Serial Only)..." << std::endl;

    int W = pA.size() * 8; 
    std::vector<double> t_ser;
    int ok_serial_vs_plain = 0;

    
    Bytes expected_not(pA.size());
    for (size_t i = 0; i < pA.size(); ++i) {
        expected_not[i] = ~pA[i];
    }

    for (int ep = 0; ep < epochs; ++ep) {
        double ms_ser = 0.0;

        
        tlwe_vec_ptr res_ser;
        {
            AutoTimer t(&ms_ser);
            res_ser = key->bit_not(ctA); 
        }

        t_ser.push_back(ms_ser);

        
        Bytes dec_ser = key->decryptBatch(res_ser);

        
        bool pass_serial_vs_plain = (dec_ser == expected_not);
        if (pass_serial_vs_plain) ok_serial_vs_plain++;

        
        if (detail) {
            std::cout << "  Epoch [" << ep + 1 << "/" << epochs << "] "
                      << "| Serial: " << ms_ser << " ms\n"
                      << "    -> Correctness: " << (pass_serial_vs_plain ? "PASS ✅" : "FAIL ❌") << std::endl;
        }
    }

    
    double avg_ser = std::accumulate(t_ser.begin(), t_ser.end(), 0.0) / epochs;

    auto percent = [](int num, int den) -> double {
        return (den > 0) ? (100.0 * num / den) : 0.0;
    };

    std::cout << "\n================= SUMMARY (FheBitNot) =================\n";
    std::cout << "Data Size     : " << pA.size() << " Bytes (" << W << " bits)\n";
    std::cout << "HomNOT (ser)  : " << avg_ser << " ms (avg)\n";
    std::cout << "HomNOT (pal)  : N/A (Unary operation, extremely fast)\n";
    std::cout << "Speedup       : N/A\n";
    std::cout << "----------------------------------------------------\n";
    std::cout << "Serial vs Plain NOT : " << ok_serial_vs_plain << "/" << epochs 
              << " (" << percent(ok_serial_vs_plain, epochs) << "%)\n";
    std::cout << "Parallel vs Serial  : N/A\n";
    std::cout << "====================================================\n";
}


// ====================================================================

// ====================================================================
void test_FheEq(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail) {
    std::cout << "\n[*] Running FheEq test (Serial vs Parallel)..." << std::endl;

    int W = pA.size() * 8; 
    int hw_threads = std::max(1u, std::thread::hardware_concurrency());

    
    std::vector<double> t_rand_ser, t_eq_ser, t_neq_ser;
    std::vector<double> t_rand_par, t_eq_par, t_neq_par;

    
    int ok_serial_vs_plain = 0;
    int ok_parallel_vs_serial = 0;
    int total_cases = epochs * 3; 

    
    Bytes pB_eq = pA; 
    Bytes pB_neq = pA;
    if (!pB_neq.empty()) pB_neq.back() ^= 1; 

    for (int ep = 0; ep < epochs; ++ep) {
        if (detail) {
            std::cout << "\n  ----- Epoch [" << ep + 1 << "/" << epochs << "] -----";
        }

        // ==========================================
        
        // ==========================================
        auto run_case = [&](const std::string& case_name, const tlwe_vec_ptr& ct_B_target, const Bytes& pB_target, 
                            std::vector<double>& t_ser_vec, std::vector<double>& t_par_vec) {
            
            
            boost::multiprecision::cpp_int a_val = bytes_to_cpp_int(pA);
            boost::multiprecision::cpp_int b_val = bytes_to_cpp_int(pB_target);
            uint8_t truth_bit = (a_val == b_val) ? 1 : 0;

            
            tlwe_ptr res_ser, res_par;
            double ms_ser = 0.0, ms_par = 0.0;
            {
                AutoTimer t(&ms_ser);
                res_ser = to_tlwe_ptr(key->eq(ctA, ct_B_target));
            }
            {
                AutoTimer t(&ms_par);
                res_par = to_tlwe_ptr(key->eq_pal(ctA, ct_B_target, hw_threads));
            }
            t_ser_vec.push_back(ms_ser);
            t_par_vec.push_back(ms_par);

            
            uint8_t dec_ser = key->decrypt(res_ser);
            uint8_t dec_par = key->decrypt(res_par);

            
            bool pass_ser_ok = (dec_ser == truth_bit);
            bool pass_par_vs_ser = (dec_ser == dec_par);

            if (pass_ser_ok) ok_serial_vs_plain++;
            if (pass_par_vs_ser) ok_parallel_vs_serial++;

            if (detail) {
                std::cout << "\n    [" << case_name << "] "
                          << "| Ser: " << ms_ser << " ms "
                          << "| Par (" << hw_threads << "T): " << ms_par << " ms\n"
                          << "      -> SerOutput: " << (int)dec_ser << " (Truth: " << (int)truth_bit << ") " 
                          << (pass_ser_ok ? "✅" : "❌")
                          << " | Par==Ser: " << (pass_par_vs_ser ? "✅" : "❌");
            }
        };

        // ==========================================
        
        // ==========================================
        
        run_case("Rand vs Rand", ctB, pB, t_rand_ser, t_rand_par);
        
        
        tlwe_vec_ptr ctB_eq = decompress_ciphertexts(key->encryptBatch(pB_eq));
        run_case("Guarantee ==", ctB_eq, pB_eq, t_eq_ser, t_eq_par);

        
        tlwe_vec_ptr ctB_neq = decompress_ciphertexts(key->encryptBatch(pB_neq));
        run_case("Guarantee !=", ctB_neq, pB_neq, t_neq_ser, t_neq_par);
        
        if (detail) std::cout << "\n";
    }

    
    
    auto calc_avg = [](const std::vector<double>& v1, const std::vector<double>& v2, const std::vector<double>& v3) {
        double sum = std::accumulate(v1.begin(), v1.end(), 0.0) + 
                     std::accumulate(v2.begin(), v2.end(), 0.0) + 
                     std::accumulate(v3.begin(), v3.end(), 0.0);
        return sum / (v1.size() + v2.size() + v3.size());
    };

    double avg_ser = calc_avg(t_rand_ser, t_eq_ser, t_neq_ser);
    double avg_par = calc_avg(t_rand_par, t_eq_par, t_neq_par);

    auto percent = [](int num, int den) -> double {
        return (den > 0) ? (100.0 * num / den) : 0.0;
    };

    std::cout << "\n================= SUMMARY (FheEq) =================\n";
    std::cout << "Data Size     : " << pA.size() << " Bytes (" << W << " bits)\n";
    std::cout << "HomEQ (ser)   : " << avg_ser << " ms (avg over all cases)\n";
    std::cout << "HomEQ (pal)   : " << avg_par << " ms (avg over all cases)\n";
    if (avg_par > 0) {
        std::cout << "Speedup       : x" << (avg_ser / avg_par) << "\n";
    }
    std::cout << "----------------------------------------------------\n";
    std::cout << "Serial vs Truth Bit : " << ok_serial_vs_plain << "/" << total_cases 
              << " (" << percent(ok_serial_vs_plain, total_cases) << "%)\n";
    std::cout << "Parallel vs Serial  : " << ok_parallel_vs_serial << "/" << total_cases 
              << " (" << percent(ok_parallel_vs_serial, total_cases) << "%)\n";
    std::cout << "====================================================\n";
}

// ====================================================================

// ====================================================================
void test_FheNe(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail) {
    std::cout << "\n[*] Running FheNe test (Serial vs Parallel)..." << std::endl;

    int W = pA.size() * 8; 
    int hw_threads = std::max(1u, std::thread::hardware_concurrency());

    
    std::vector<double> t_rand_ser, t_eq_ser, t_neq_ser;
    std::vector<double> t_rand_par, t_eq_par, t_neq_par;

    
    int ok_serial_vs_plain = 0;
    int ok_parallel_vs_serial = 0;
    int total_cases = epochs * 3; 

    
    Bytes pB_eq = pA; 
    Bytes pB_neq = pA;
    if (!pB_neq.empty()) pB_neq.back() ^= 1; 

    for (int ep = 0; ep < epochs; ++ep) {
        if (detail) {
            std::cout << "\n  ----- Epoch [" << ep + 1 << "/" << epochs << "] -----";
        }

        // ==========================================
        
        // ==========================================
        auto run_case = [&](const std::string& case_name, const tlwe_vec_ptr& ct_B_target, const Bytes& pB_target, 
                            std::vector<double>& t_ser_vec, std::vector<double>& t_par_vec) {
            
            
            boost::multiprecision::cpp_int a_val = bytes_to_cpp_int(pA);
            boost::multiprecision::cpp_int b_val = bytes_to_cpp_int(pB_target);
            uint8_t truth_bit = (a_val != b_val) ? 1 : 0;

            
            tlwe_ptr res_ser, res_par;
            double ms_ser = 0.0, ms_par = 0.0;
            {
                AutoTimer t(&ms_ser);
                res_ser = to_tlwe_ptr(key->ne(ctA, ct_B_target));
            }
            {
                AutoTimer t(&ms_par);
                res_par = to_tlwe_ptr(key->ne_pal(ctA, ct_B_target, hw_threads));
            }
            t_ser_vec.push_back(ms_ser);
            t_par_vec.push_back(ms_par);

            
            uint8_t dec_ser = key->decrypt(res_ser);
            uint8_t dec_par = key->decrypt(res_par);

            
            bool pass_ser_ok = (dec_ser == truth_bit);
            bool pass_par_vs_ser = (dec_ser == dec_par);

            if (pass_ser_ok) ok_serial_vs_plain++;
            if (pass_par_vs_ser) ok_parallel_vs_serial++;

            if (detail) {
                std::cout << "\n    [" << case_name << "] "
                          << "| Ser: " << ms_ser << " ms "
                          << "| Par (" << hw_threads << "T): " << ms_par << " ms\n"
                          << "      -> SerOutput: " << (int)dec_ser << " (Truth: " << (int)truth_bit << ") " 
                          << (pass_ser_ok ? "✅" : "❌")
                          << " | Par==Ser: " << (pass_par_vs_ser ? "✅" : "❌");
            }
        };

        // ==========================================
        
        // ==========================================
        
        run_case("Rand vs Rand", ctB, pB, t_rand_ser, t_rand_par);
        
        
        tlwe_vec_ptr ctB_eq = decompress_ciphertexts(key->encryptBatch(pB_eq));
        run_case("Guarantee == (Expect 0)", ctB_eq, pB_eq, t_eq_ser, t_eq_par);

        
        tlwe_vec_ptr ctB_neq = decompress_ciphertexts(key->encryptBatch(pB_neq));
        run_case("Guarantee != (Expect 1)", ctB_neq, pB_neq, t_neq_ser, t_neq_par);
        
        if (detail) std::cout << "\n";
    }

    
    auto calc_avg = [](const std::vector<double>& v1, const std::vector<double>& v2, const std::vector<double>& v3) {
        double sum = std::accumulate(v1.begin(), v1.end(), 0.0) + 
                     std::accumulate(v2.begin(), v2.end(), 0.0) + 
                     std::accumulate(v3.begin(), v3.end(), 0.0);
        return sum / (v1.size() + v2.size() + v3.size());
    };

    double avg_ser = calc_avg(t_rand_ser, t_eq_ser, t_neq_ser);
    double avg_par = calc_avg(t_rand_par, t_eq_par, t_neq_par);

    auto percent = [](int num, int den) -> double {
        return (den > 0) ? (100.0 * num / den) : 0.0;
    };

    std::cout << "\n================= SUMMARY (FheNe) =================\n";
    std::cout << "Data Size     : " << pA.size() << " Bytes (" << W << " bits)\n";
    std::cout << "HomNE (ser)   : " << avg_ser << " ms (avg over all cases)\n";
    std::cout << "HomNE (pal)   : " << avg_par << " ms (avg over all cases)\n";
    if (avg_par > 0) {
        std::cout << "Speedup       : x" << (avg_ser / avg_par) << "\n";
    }
    std::cout << "----------------------------------------------------\n";
    std::cout << "Serial vs Truth Bit : " << ok_serial_vs_plain << "/" << total_cases 
              << " (" << percent(ok_serial_vs_plain, total_cases) << "%)\n";
    std::cout << "Parallel vs Serial  : " << ok_parallel_vs_serial << "/" << total_cases 
              << " (" << percent(ok_parallel_vs_serial, total_cases) << "%)\n";
    std::cout << "====================================================\n";
}

// ====================================================================

// ====================================================================
void test_FheGe(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail) {
    std::cout << "\n[*] Running FheGe (Greater Than or Equal, Signed) test (Serial vs Parallel)..." << std::endl;

    int W = pA.size() * 8; 
    int hw_threads = std::max(1u, std::thread::hardware_concurrency());

    
    std::vector<double> t_ser, t_par;

    
    int ok_serial_vs_plain = 0;
    int ok_parallel_vs_serial = 0;
    int total_cases = epochs * 3; 

    for (int ep = 0; ep < epochs; ++ep) {
        if (detail) {
            std::cout << "\n  ----- Epoch [" << ep + 1 << "/" << epochs << "] -----";
        }

        // ==========================================
        
        // ==========================================
        auto run_case = [&](const std::string& case_name, const tlwe_vec_ptr& ct_A_target, const tlwe_vec_ptr& ct_B_target, 
                            const Bytes& pA_target, const Bytes& pB_target) {
            
            
            boost::multiprecision::cpp_int a_val = bytes_to_signed_cpp_int(pA_target);
            boost::multiprecision::cpp_int b_val = bytes_to_signed_cpp_int(pB_target);
            uint8_t truth_bit = (a_val >= b_val) ? 1 : 0;

            
            tlwe_ptr res_ser, res_par;
            double ms_ser = 0.0, ms_par = 0.0;
            {
                AutoTimer t(&ms_ser);
                res_ser = to_tlwe_ptr(key->ge(ct_A_target, ct_B_target));
            }
            {
                AutoTimer t(&ms_par);
                res_par = to_tlwe_ptr(key->ge_pal(ct_A_target, ct_B_target, hw_threads));
            }
            t_ser.push_back(ms_ser);
            t_par.push_back(ms_par);

            
            uint8_t dec_ser = key->decrypt(res_ser);
            uint8_t dec_par = key->decrypt(res_par);

            
            bool pass_ser_ok = (dec_ser == truth_bit);
            bool pass_par_vs_ser = (dec_ser == dec_par);

            if (pass_ser_ok) ok_serial_vs_plain++;
            if (pass_par_vs_ser) ok_parallel_vs_serial++;

            if (detail) {
                std::cout << "\n    [" << case_name << "]"
                          << "\n      A(signed): " << a_val << ", B(signed): " << b_val
                          << "\n      | Ser: " << ms_ser << " ms | Par (" << hw_threads << "T): " << ms_par << " ms"
                          << "\n      -> SerOutput: " << (int)dec_ser << " (Truth: " << (int)truth_bit << ") " 
                          << (pass_ser_ok ? "✅" : "❌")
                          << " | Par==Ser: " << (pass_par_vs_ser ? "✅" : "❌") << "\n";
            }
        };

        // ==========================================
        
        // ==========================================
        
        run_case("Random vs Random", ctA, ctB, pA, pB);
        
        
        
        Bytes pA_gt = pA; pA_gt[0] &= 0x7F; 
        Bytes pB_lt = pB; pB_lt[0] |= 0x80; 
        tlwe_vec_ptr ctA_gt = decompress_ciphertexts(key->encryptBatch(pA_gt));
        tlwe_vec_ptr ctB_lt = decompress_ciphertexts(key->encryptBatch(pB_lt));
        run_case("Guaranteed A > B (Positive >= Negative, Output 1)", ctA_gt, ctB_lt, pA_gt, pB_lt);

        
        
        run_case("Guaranteed A == B (Should output 1)", ctA, ctA, pA, pA);
    }

    
    double avg_ser = std::accumulate(t_ser.begin(), t_ser.end(), 0.0) / t_ser.size();
    double avg_par = std::accumulate(t_par.begin(), t_par.end(), 0.0) / t_par.size();

    auto percent = [](int num, int den) -> double {
        return (den > 0) ? (100.0 * num / den) : 0.0;
    };

    std::cout << "\n================= SUMMARY (FheGe) =================\n";
    std::cout << "Data Size     : " << pA.size() << " Bytes (" << W << " bits)\n";
    std::cout << "HomGE (ser)   : " << avg_ser << " ms (avg over all cases)\n";
    std::cout << "HomGE (pal)   : " << avg_par << " ms (avg over all cases)\n";
    if (avg_par > 0) {
        std::cout << "Speedup       : x" << (avg_ser / avg_par) << "\n";
    }
    std::cout << "----------------------------------------------------\n";
    std::cout << "Serial vs Truth Bit : " << ok_serial_vs_plain << "/" << total_cases 
              << " (" << percent(ok_serial_vs_plain, total_cases) << "%)\n";
    std::cout << "Parallel vs Serial  : " << ok_parallel_vs_serial << "/" << total_cases 
              << " (" << percent(ok_parallel_vs_serial, total_cases) << "%)\n";
    std::cout << "====================================================\n";
}

// ====================================================================

// ====================================================================
void test_FheGt(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail) {
    std::cout << "\n[*] Running FheGt (Greater Than, Signed) test (Serial vs Parallel)..." << std::endl;

    int W = pA.size() * 8; 
    int hw_threads = std::max(1u, std::thread::hardware_concurrency());

    
    std::vector<double> t_ser, t_par;

    
    int ok_serial_vs_plain = 0;
    int ok_parallel_vs_serial = 0;
    int total_cases = epochs * 3; 

    for (int ep = 0; ep < epochs; ++ep) {
        if (detail) {
            std::cout << "\n  ----- Epoch [" << ep + 1 << "/" << epochs << "] -----";
        }

        // ==========================================
        
        // ==========================================
        auto run_case = [&](const std::string& case_name, const tlwe_vec_ptr& ct_A_target, const tlwe_vec_ptr& ct_B_target, 
                            const Bytes& pA_target, const Bytes& pB_target) {
            
            
            boost::multiprecision::cpp_int a_val = bytes_to_signed_cpp_int(pA_target);
            boost::multiprecision::cpp_int b_val = bytes_to_signed_cpp_int(pB_target);
            uint8_t truth_bit = (a_val > b_val) ? 1 : 0;

            
            tlwe_ptr res_ser, res_par;
            double ms_ser = 0.0, ms_par = 0.0;
            {
                AutoTimer t(&ms_ser);
                res_ser = to_tlwe_ptr(key->gt(ct_A_target, ct_B_target));
            }
            {
                AutoTimer t(&ms_par);
                res_par = to_tlwe_ptr(key->gt_pal(ct_A_target, ct_B_target, hw_threads));
            }
            t_ser.push_back(ms_ser);
            t_par.push_back(ms_par);

            
            uint8_t dec_ser = key->decrypt(res_ser);
            uint8_t dec_par = key->decrypt(res_par);

            
            bool pass_ser_ok = (dec_ser == truth_bit);
            bool pass_par_vs_ser = (dec_ser == dec_par);

            if (pass_ser_ok) ok_serial_vs_plain++;
            if (pass_par_vs_ser) ok_parallel_vs_serial++;

            if (detail) {
                std::cout << "\n    [" << case_name << "]"
                          << "\n      A(signed): " << a_val << ", B(signed): " << b_val
                          << "\n      | Ser: " << ms_ser << " ms | Par (" << hw_threads << "T): " << ms_par << " ms"
                          << "\n      -> SerOutput: " << (int)dec_ser << " (Truth: " << (int)truth_bit << ") " 
                          << (pass_ser_ok ? "✅" : "❌")
                          << " | Par==Ser: " << (pass_par_vs_ser ? "✅" : "❌") << "\n";
            }
        };

        // ==========================================
        
        // ==========================================
        
        run_case("Random vs Random", ctA, ctB, pA, pB);
        
        
        Bytes pA_gt = pA; pA_gt[0] &= 0x7F; // MSB = 0
        Bytes pB_lt = pB; pB_lt[0] |= 0x80; // MSB = 1
        tlwe_vec_ptr ctA_gt = decompress_ciphertexts(key->encryptBatch(pA_gt));
        tlwe_vec_ptr ctB_lt = decompress_ciphertexts(key->encryptBatch(pB_lt));
        run_case("Guaranteed A > B (Positive > Negative)", ctA_gt, ctB_lt, pA_gt, pB_lt);

        
        run_case("Guaranteed A == B (Should output 0)", ctA, ctA, pA, pA);
    }

    
    double avg_ser = std::accumulate(t_ser.begin(), t_ser.end(), 0.0) / t_ser.size();
    double avg_par = std::accumulate(t_par.begin(), t_par.end(), 0.0) / t_par.size();

    auto percent = [](int num, int den) -> double {
        return (den > 0) ? (100.0 * num / den) : 0.0;
    };

    std::cout << "\n================= SUMMARY (FheGt) =================\n";
    std::cout << "Data Size     : " << pA.size() << " Bytes (" << W << " bits)\n";
    std::cout << "HomGT (ser)   : " << avg_ser << " ms (avg over all cases)\n";
    std::cout << "HomGT (pal)   : " << avg_par << " ms (avg over all cases)\n";
    if (avg_par > 0) {
        std::cout << "Speedup       : x" << (avg_ser / avg_par) << "\n";
    }
    std::cout << "----------------------------------------------------\n";
    std::cout << "Serial vs Truth Bit : " << ok_serial_vs_plain << "/" << total_cases 
              << " (" << percent(ok_serial_vs_plain, total_cases) << "%)\n";
    std::cout << "Parallel vs Serial  : " << ok_parallel_vs_serial << "/" << total_cases 
              << " (" << percent(ok_parallel_vs_serial, total_cases) << "%)\n";
    std::cout << "====================================================\n";
}

// ====================================================================

// ====================================================================
void test_FheLe(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail) {
    std::cout << "\n[*] Running FheLe (Less Than or Equal, Signed) test (Serial vs Parallel)..." << std::endl;

    int W = pA.size() * 8; 
    int hw_threads = std::max(1u, std::thread::hardware_concurrency());

    
    std::vector<double> t_ser, t_par;

    
    int ok_serial_vs_plain = 0;
    int ok_parallel_vs_serial = 0;
    int total_cases = epochs * 3; 

    for (int ep = 0; ep < epochs; ++ep) {
        if (detail) {
            std::cout << "\n  ----- Epoch [" << ep + 1 << "/" << epochs << "] -----";
        }

        // ==========================================
        
        // ==========================================
        auto run_case = [&](const std::string& case_name, const tlwe_vec_ptr& ct_A_target, const tlwe_vec_ptr& ct_B_target, 
                            const Bytes& pA_target, const Bytes& pB_target) {
            
            
            boost::multiprecision::cpp_int a_val = bytes_to_signed_cpp_int(pA_target);
            boost::multiprecision::cpp_int b_val = bytes_to_signed_cpp_int(pB_target);
            uint8_t truth_bit = (a_val <= b_val) ? 1 : 0;

            
            tlwe_ptr res_ser, res_par;
            double ms_ser = 0.0, ms_par = 0.0;
            {
                AutoTimer t(&ms_ser);
                res_ser = to_tlwe_ptr(key->le(ct_A_target, ct_B_target));
            }
            {
                AutoTimer t(&ms_par);
                res_par = to_tlwe_ptr(key->le_pal(ct_A_target, ct_B_target, hw_threads));
            }
            t_ser.push_back(ms_ser);
            t_par.push_back(ms_par);

            
            uint8_t dec_ser = key->decrypt(res_ser);
            uint8_t dec_par = key->decrypt(res_par);

            
            bool pass_ser_ok = (dec_ser == truth_bit);
            bool pass_par_vs_ser = (dec_ser == dec_par);

            if (pass_ser_ok) ok_serial_vs_plain++;
            if (pass_par_vs_ser) ok_parallel_vs_serial++;

            if (detail) {
                std::cout << "\n    [" << case_name << "]"
                          << "\n      A(signed): " << a_val << ", B(signed): " << b_val
                          << "\n      | Ser: " << ms_ser << " ms | Par (" << hw_threads << "T): " << ms_par << " ms"
                          << "\n      -> SerOutput: " << (int)dec_ser << " (Truth: " << (int)truth_bit << ") " 
                          << (pass_ser_ok ? "✅" : "❌")
                          << " | Par==Ser: " << (pass_par_vs_ser ? "✅" : "❌") << "\n";
            }
        };

        // ==========================================
        
        // ==========================================
        
        run_case("Random vs Random", ctA, ctB, pA, pB);
        
        
        
        Bytes pA_lt = pA; pA_lt[0] |= 0x80; 
        Bytes pB_gt = pB; pB_gt[0] &= 0x7F; 
        tlwe_vec_ptr ctA_lt = decompress_ciphertexts(key->encryptBatch(pA_lt));
        tlwe_vec_ptr ctB_gt = decompress_ciphertexts(key->encryptBatch(pB_gt));
        run_case("Guaranteed A < B (Negative <= Positive, Output 1)", ctA_lt, ctB_gt, pA_lt, pB_gt);

        
        
        run_case("Guaranteed A == B (Should output 1)", ctA, ctA, pA, pA);
    }

    
    double avg_ser = std::accumulate(t_ser.begin(), t_ser.end(), 0.0) / t_ser.size();
    double avg_par = std::accumulate(t_par.begin(), t_par.end(), 0.0) / t_par.size();

    auto percent = [](int num, int den) -> double {
        return (den > 0) ? (100.0 * num / den) : 0.0;
    };

    std::cout << "\n================= SUMMARY (FheLe) =================\n";
    std::cout << "Data Size     : " << pA.size() << " Bytes (" << W << " bits)\n";
    std::cout << "HomLE (ser)   : " << avg_ser << " ms (avg over all cases)\n";
    std::cout << "HomLE (pal)   : " << avg_par << " ms (avg over all cases)\n";
    if (avg_par > 0) {
        std::cout << "Speedup       : x" << (avg_ser / avg_par) << "\n";
    }
    std::cout << "----------------------------------------------------\n";
    std::cout << "Serial vs Truth Bit : " << ok_serial_vs_plain << "/" << total_cases 
              << " (" << percent(ok_serial_vs_plain, total_cases) << "%)\n";
    std::cout << "Parallel vs Serial  : " << ok_parallel_vs_serial << "/" << total_cases 
              << " (" << percent(ok_parallel_vs_serial, total_cases) << "%)\n";
    std::cout << "====================================================\n";
}

// ====================================================================

// ====================================================================
void test_FheLt(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const Bytes& pA, const Bytes& pB, int epochs, int detail) {
    std::cout << "\n[*] Running FheLt (Less Than, Signed) test (Serial vs Parallel)..." << std::endl;

    int W = pA.size() * 8; 
    int hw_threads = std::max(1u, std::thread::hardware_concurrency());

    
    std::vector<double> t_ser, t_par;

    
    int ok_serial_vs_plain = 0;
    int ok_parallel_vs_serial = 0;
    int total_cases = epochs * 3; 

    for (int ep = 0; ep < epochs; ++ep) {
        if (detail) {
            std::cout << "\n  ----- Epoch [" << ep + 1 << "/" << epochs << "] -----";
        }

        // ==========================================
        
        // ==========================================
        auto run_case = [&](const std::string& case_name, const tlwe_vec_ptr& ct_A_target, const tlwe_vec_ptr& ct_B_target, 
                            const Bytes& pA_target, const Bytes& pB_target) {
            
            
            boost::multiprecision::cpp_int a_val = bytes_to_signed_cpp_int(pA_target);
            boost::multiprecision::cpp_int b_val = bytes_to_signed_cpp_int(pB_target);
            uint8_t truth_bit = (a_val < b_val) ? 1 : 0;

            
            tlwe_ptr res_ser, res_par;
            double ms_ser = 0.0, ms_par = 0.0;
            {
                AutoTimer t(&ms_ser);
                res_ser = to_tlwe_ptr(key->lt(ct_A_target, ct_B_target));
            }
            {
                AutoTimer t(&ms_par);
                res_par = to_tlwe_ptr(key->lt_pal(ct_A_target, ct_B_target, hw_threads));
            }
            t_ser.push_back(ms_ser);
            t_par.push_back(ms_par);

            
            uint8_t dec_ser = key->decrypt(res_ser);
            uint8_t dec_par = key->decrypt(res_par);

            
            bool pass_ser_ok = (dec_ser == truth_bit);
            bool pass_par_vs_ser = (dec_ser == dec_par);

            if (pass_ser_ok) ok_serial_vs_plain++;
            if (pass_par_vs_ser) ok_parallel_vs_serial++;

            if (detail) {
                std::cout << "\n    [" << case_name << "]"
                          << "\n      A(signed): " << a_val << ", B(signed): " << b_val
                          << "\n      | Ser: " << ms_ser << " ms | Par (" << hw_threads << "T): " << ms_par << " ms"
                          << "\n      -> SerOutput: " << (int)dec_ser << " (Truth: " << (int)truth_bit << ") " 
                          << (pass_ser_ok ? "✅" : "❌")
                          << " | Par==Ser: " << (pass_par_vs_ser ? "✅" : "❌") << "\n";
            }
        };

        // ==========================================
        
        // ==========================================
        
        run_case("Random vs Random", ctA, ctB, pA, pB);
        
        
        Bytes pA_lt = pA; pA_lt[0] |= 0x80; 
        Bytes pB_gt = pB; pB_gt[0] &= 0x7F; 
        tlwe_vec_ptr ctA_lt = decompress_ciphertexts(key->encryptBatch(pA_lt));
        tlwe_vec_ptr ctB_gt = decompress_ciphertexts(key->encryptBatch(pB_gt));
        run_case("Guaranteed A < B (Negative < Positive)", ctA_lt, ctB_gt, pA_lt, pB_gt);

        
        run_case("Guaranteed A == B (Should output 0)", ctA, ctA, pA, pA);
    }

    
    double avg_ser = std::accumulate(t_ser.begin(), t_ser.end(), 0.0) / t_ser.size();
    double avg_par = std::accumulate(t_par.begin(), t_par.end(), 0.0) / t_par.size();

    auto percent = [](int num, int den) -> double {
        return (den > 0) ? (100.0 * num / den) : 0.0;
    };

    std::cout << "\n================= SUMMARY (FheLt) =================\n";
    std::cout << "Data Size     : " << pA.size() << " Bytes (" << W << " bits)\n";
    std::cout << "HomLT (ser)   : " << avg_ser << " ms (avg over all cases)\n";
    std::cout << "HomLT (pal)   : " << avg_par << " ms (avg over all cases)\n";
    if (avg_par > 0) {
        std::cout << "Speedup       : x" << (avg_ser / avg_par) << "\n";
    }
    std::cout << "----------------------------------------------------\n";
    std::cout << "Serial vs Truth Bit : " << ok_serial_vs_plain << "/" << total_cases 
              << " (" << percent(ok_serial_vs_plain, total_cases) << "%)\n";
    std::cout << "Parallel vs Serial  : " << ok_parallel_vs_serial << "/" << total_cases 
              << " (" << percent(ok_parallel_vs_serial, total_cases) << "%)\n";
    std::cout << "====================================================\n";
}


// ====================================================================


// ====================================================================
void test_FheIfThenElse(FheKey* key, const tlwe_vec_ptr& ctA, const tlwe_vec_ptr& ctB, const tlwe_ptr& ctCond, 
                        const Bytes& pA, const Bytes& pB, uint8_t pCond, int epochs, int detail) {
    std::cout << "\n[*] Running FheIfThenElse (Select) test (Serial vs Parallel)..." << std::endl;

    int W = pA.size() * 8; 
    int hw_threads = std::max(1u, std::thread::hardware_concurrency());

    std::vector<double> t_ser, t_par;
    int ok_serial_vs_plain = 0;
    int ok_parallel_vs_serial = 0;

    
    Bytes expected_out = pCond ? pB : pA;

    if (detail) {
        std::cout << "  [Info] Plaintext Selection Bit (pCond) = " << (int)pCond 
                  << " (Expecting to output " << (pCond ? "B" : "A") << ")\n";
    }

    for (int ep = 0; ep < epochs; ++ep) {
        double ms_ser = 0.0, ms_par = 0.0;

        
        tlwe_vec_ptr res_ser;
        {
            AutoTimer t(&ms_ser);
            res_ser = key->select(ctA, ctB, to_tlwe_vec_ptr(ctCond), 1); 
        }

        
        tlwe_vec_ptr res_par;
        {
            AutoTimer t(&ms_par);
            res_par = key->select(ctA, ctB, to_tlwe_vec_ptr(ctCond), hw_threads);
        }

        t_ser.push_back(ms_ser);
        t_par.push_back(ms_par);

        
        Bytes dec_ser = key->decryptBatch(res_ser);
        Bytes dec_par = key->decryptBatch(res_par);

        
        bool pass_serial_vs_plain = (dec_ser == expected_out);
        bool pass_parallel_vs_serial = (dec_ser == dec_par); 

        if (pass_serial_vs_plain) ok_serial_vs_plain++;
        if (pass_parallel_vs_serial) ok_parallel_vs_serial++;

        
        if (detail) {
            std::cout << "  Epoch [" << ep + 1 << "/" << epochs << "] "
                      << "| Serial: " << ms_ser << " ms "
                      << "| Parallel (" << hw_threads << "T): " << ms_par << " ms\n"
                      << "    -> Correctness: " << (pass_serial_vs_plain ? "PASS ✅" : "FAIL ❌") 
                      << " | Par==Ser: " << (pass_parallel_vs_serial ? "PASS ✅" : "FAIL ❌") << std::endl;
        }
    }

    
    double avg_ser = std::accumulate(t_ser.begin(), t_ser.end(), 0.0) / epochs;
    double avg_par = std::accumulate(t_par.begin(), t_par.end(), 0.0) / epochs;

    auto percent = [](int num, int den) -> double {
        return (den > 0) ? (100.0 * num / den) : 0.0;
    };

    std::cout << "\n================= SUMMARY (FheSelect) =================\n";
    std::cout << "Data Size     : " << pA.size() << " Bytes (" << W << " bits)\n";
    std::cout << "HomSelect(ser): " << avg_ser << " ms (avg)\n";
    std::cout << "HomSelect(pal): " << avg_par << " ms (avg)\n";
    if (avg_par > 0) {
        std::cout << "Speedup       : x" << (avg_ser / avg_par) << "\n";
    }
    std::cout << "----------------------------------------------------\n";
    std::cout << "Serial vs Plain MUX : " << ok_serial_vs_plain << "/" << epochs 
              << " (" << percent(ok_serial_vs_plain, epochs) << "%)\n";
    std::cout << "Parallel vs Serial  : " << ok_parallel_vs_serial << "/" << epochs 
              << " (" << percent(ok_parallel_vs_serial, epochs) << "%)\n";
    std::cout << "====================================================\n";
}

// ====================================================================

// ====================================================================
void test_TrivialEncrypt(FheKey* key, const Bytes& pA, int epochs, int detail) {
    std::cout << "\n[*] Running TrivialEncrypt test (Trivial vs Standard)..." << std::endl;

    int ByteCount = pA.size();
    int W = ByteCount * 8; 

    
    std::vector<double> t_triv_single, t_std_single;
    std::vector<double> t_triv_batch, t_std_batch;

    int ok_correctness = 0;
    int ok_determinism = 0;

    for (int ep = 0; ep < epochs; ++ep) {
        uint8_t plain_bit = ep % 2; 

        double ms_triv_s = 0.0, ms_std_s = 0.0;
        double ms_triv_b = 0.0, ms_std_b = 0.0;

        
        tlwe_ptr ct_triv1;
        { AutoTimer t(&ms_triv_s); ct_triv1 = to_tlwe_ptr(decompress_ciphertexts(key->trivialEncrypt(plain_bit))); }
        tlwe_ptr ct_triv2 = to_tlwe_ptr(decompress_ciphertexts(key->trivialEncrypt(plain_bit))); 
        
        tlwe_ptr ct_std;
        { AutoTimer t(&ms_std_s);  ct_std = to_tlwe_ptr(decompress_ciphertexts(key->encrypt(plain_bit))); }

        
        tlwe_vec_ptr vec_triv1;
        { AutoTimer t(&ms_triv_b); vec_triv1 = decompress_ciphertexts(key->trivialEncryptBatch(pA)); }
        tlwe_vec_ptr vec_triv2 = decompress_ciphertexts(key->trivialEncryptBatch(pA)); 

        tlwe_vec_ptr vec_std;
        { AutoTimer t(&ms_std_b);  vec_std = decompress_ciphertexts(key->encryptBatch(pA)); }

        
        Bytes decrypted_bytes = key->decryptBatch(vec_triv1);

        
        bool single_correct = (key->decrypt(ct_triv1) == plain_bit);
        bool single_det = are_ciphertexts_identical(ct_triv1, ct_triv2);
        
        bool batch_correct = (pA == decrypted_bytes);
        bool batch_det = are_ciphertext_vecs_identical(vec_triv1, vec_triv2);

        bool epoch_correct = single_correct && batch_correct;
        bool epoch_det = single_det && batch_det;

        if (epoch_correct) ok_correctness++;
        if (epoch_det) ok_determinism++;

        t_triv_single.push_back(ms_triv_s);
        t_std_single.push_back(ms_std_s);
        t_triv_batch.push_back(ms_triv_b);
        t_std_batch.push_back(ms_std_b);

        if (detail) {
            std::cout << "  Epoch [" << ep + 1 << "/" << epochs << "] "
                      << "| TrivBatch: " << ms_triv_b << " ms "
                      << "| StdBatch: " << ms_std_b << " ms\n"
                      << "    -> Correctness: " << (epoch_correct ? "PASS ✅" : "FAIL ❌") 
                      << " | Determinism: " << (epoch_det ? "PASS ✅" : "FAIL ❌") << std::endl;
        }
    }

    
    double avg_triv_s = std::accumulate(t_triv_single.begin(), t_triv_single.end(), 0.0) / epochs;
    double avg_std_s  = std::accumulate(t_std_single.begin(), t_std_single.end(), 0.0) / epochs;
    double avg_triv_b = std::accumulate(t_triv_batch.begin(), t_triv_batch.end(), 0.0) / epochs;
    double avg_std_b  = std::accumulate(t_std_batch.begin(), t_std_batch.end(), 0.0) / epochs;

    auto percent = [](int num, int den) -> double {
        return (den > 0) ? (100.0 * num / den) : 0.0;
    };

    std::cout << "\n================= SUMMARY (Trivial Encrypt) =================\n";
    std::cout << "Data Size     : " << ByteCount << " Bytes (" << W << " bits)\n";
    
    std::cout << "Std Enc  (1b) : " << avg_std_s << " ms (avg)\n";
    std::cout << "Triv Enc (1b) : " << avg_triv_s << " ms (avg)\n";
    if (avg_triv_s > 0) {
        std::cout << "Speedup  (1b) : x" << (avg_std_s / avg_triv_s) << "\n";
    }
    std::cout << "----------------------------------------------------\n";
    
    std::cout << "Std Enc (Bat) : " << avg_std_b << " ms (avg)\n";
    std::cout << "Triv Enc(Bat) : " << avg_triv_b << " ms (avg)\n";
    if (avg_triv_b > 0) {
        std::cout << "Speedup (Bat) : x" << (avg_std_b / avg_triv_b) << "\n";
    }
    std::cout << "----------------------------------------------------\n";
    
    std::cout << "Correctness   : " << ok_correctness << "/" << epochs 
              << " (" << percent(ok_correctness, epochs) << "%)\n";
    std::cout << "Determinism   : " << ok_determinism << "/" << epochs 
              << " (" << percent(ok_determinism, epochs) << "%)\n";
    std::cout << "====================================================\n";
}
