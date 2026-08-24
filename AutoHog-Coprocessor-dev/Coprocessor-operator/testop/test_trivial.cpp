#include "Autohog_Cop.hpp" 
#include "timer.h"
#include "fhe_key.h"

#include <tfhe++.hpp>
#include <random>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cstring> // for memcmp

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace TFHEpp;


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
    std::cerr << "Usage: " << prog << " <byte-count> [epochs] "
              << "[--no-epoch-detail|--show-epochs]\n";
}


bool are_ciphertexts_identical(const tlwe_ptr& a, const tlwe_ptr& b) {
    
    if (a->data.size() != b->data.size()) return false;
    
    
    return std::memcmp(a->data.data(), b->data.data(), sizeof(a->data)) == 0;
}


bool are_ciphertexts_identical(const trlwe_ptr& a, const trlwe_ptr& b) {
    auto serialized_a = serialize_compressed_ciphertexts(a);
    auto serialized_b = serialize_compressed_ciphertexts(b);
    return serialized_a == serialized_b;
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
        int ByteCount = std::atoi(argv[1]); 
        if (ByteCount <= 0)
        {
            std::cerr << "Error: byte-count must be positive.\n";
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
            if (s == "--no-epoch-detail") show_epoch_detail = false;
            else if (s == "--show-epochs") show_epoch_detail = true;
        }

#ifdef _OPENMP
        {
            int hw = std::max(1u, std::thread::hardware_concurrency());
            omp_set_num_threads(hw);
        }
#endif

        
        std::vector<double> t_triv_enc_single, t_triv_enc_batch;
        std::vector<double> t_std_enc_single; 
        std::vector<double> t_dec_batch;

        int ok_correctness = 0;
        int ok_determinism = 0;
        int total_epochs = epochs;

        
        std::cout << "Initializing FheKey (KeyGen)..." << std::endl;
        auto key = create_fhe_key(); 
        std::cout << "KeyGen finished." << std::endl;

        std::default_random_engine eng(std::random_device{}());
        std::uniform_int_distribution<int> rand_byte(0, 255);

        for (int ep = 0; ep < epochs; ++ep)
        {
            if (show_epoch_detail)
            {
                std::cout << "\n================= EPOCH " << (ep + 1) << "/" << epochs << " =================\n";
            }

            
            double t_single_ms = 0.0;
            double t_batch_ms = 0.0;
            double t_std_ms = 0.0;
            double t_dec_ms = 0.0;

            // ==========================================
            
            // ==========================================
            uint8_t plain_bit = ep % 2; 
            
            
            trlwe_ptr ct1;
            {
                AutoTimer t(&t_single_ms);
                ct1 = key->trivialEncrypt(plain_bit);
            }

            
            trlwe_ptr ct2 = key->trivialEncrypt(plain_bit);

            
            trlwe_ptr ct_std;
            {
                AutoTimer t(&t_std_ms);
                ct_std = key->encrypt(plain_bit);
            }

            
            tlwe_ptr ct1_0(new tlwe(decompress_ciphertexts(ct1)->data[0]));
            bool single_correct = (key->decrypt(ct1_0) == plain_bit);
            bool single_det = are_ciphertexts_identical(ct1, ct2);
            bool std_diff = !are_ciphertexts_identical(ct1, ct_std); 

            if (show_epoch_detail) {
                std::cout << "[Single Bit " << (int)plain_bit << "] Correctness: " << (single_correct ? "OK" : "FAIL") << "\n";
                std::cout << "[Single Bit " << (int)plain_bit << "] Determinism: " << (single_det ? "OK (Identical)" : "FAIL (Different)") << "\n";
                if (!single_det) std::cout << "  ❌ Error: Trivial encrypt must produce identical ciphertexts!\n";
            }

            // ==========================================
            
            // ==========================================
            std::vector<uint8_t> plain_bytes(ByteCount);
            for(int i=0; i<ByteCount; ++i) plain_bytes[i] = rand_byte(eng);

            
            trlwe_ptr vec1;
            {
                AutoTimer t(&t_batch_ms);
                vec1 = key->trivialEncryptBatch(plain_bytes);
            }

            
            trlwe_ptr vec2 = key->trivialEncryptBatch(plain_bytes);

            
            std::vector<uint8_t> decrypted_bytes;
            {
                AutoTimer t(&t_dec_ms);
                decrypted_bytes = key->decryptBatch(decompress_ciphertexts(vec1));
            }

            
            bool batch_correct = (plain_bytes == decrypted_bytes);
            bool batch_det = are_ciphertexts_identical(vec1, vec2);

            if (show_epoch_detail) {
                std::cout << "[Batch " << ByteCount << " Bytes] Correctness: " << (batch_correct ? "OK" : "FAIL") << "\n";
                std::cout << "[Batch " << ByteCount << " Bytes] Determinism: " << (batch_det ? "OK" : "FAIL") << "\n";
            }

            // ==========================================
            
            // ==========================================
            if (single_correct && single_det && batch_correct && batch_det) {
                ok_correctness++;
                ok_determinism++;
            } else {
                if (single_correct && batch_correct) ok_correctness++;
                // Determinism fail logic handled implicitly
            }
            
            t_triv_enc_single.push_back(t_single_ms);
            t_triv_enc_batch.push_back(t_batch_ms);
            t_std_enc_single.push_back(t_std_ms);
            t_dec_batch.push_back(t_dec_ms);
        }

        
        if (epochs >= 1)
        {
            std::cout << "\n================= SUMMARY (average over " << epochs << " epochs) =================\n";
            std::cout << "\n== Timing AVG (ms) ==\n";
            std::cout << "Trivial Encrypt (1 bit):   " << avg(t_triv_enc_single) << " ms\n";
            std::cout << "Standard Encrypt (1 bit):  " << avg(t_std_enc_single) << " ms\n";
            std::cout << "  -> Speedup (Std/Triv):   x" << (avg(t_std_enc_single)/avg(t_triv_enc_single)) << "\n";
            std::cout << "--------------------------------\n";
            std::cout << "Trivial Batch (" << ByteCount << " B):  " << avg(t_triv_enc_batch) << " ms\n";
            std::cout << "Decrypt Batch (" << ByteCount << " B):  " << avg(t_dec_batch) << " ms\n";

            auto percent = [](int num, int den) -> double {
                return (den > 0) ? (100.0 * num / den) : 0.0;
            };

            std::cout << "\n== Verification Results ==\n";
            std::cout << "All Correctness Checks: " << ok_correctness << "/" << total_epochs 
                      << " (" << percent(ok_correctness, total_epochs) << "%)\n";
            std::cout << "All Determinism Checks: " << ok_determinism << "/" << total_epochs 
                      << " (" << percent(ok_determinism, total_epochs) << "%)\n";
            
            if (ok_determinism == total_epochs) {
                std::cout << "\n✅ SUCCESS: Trivial Encryption is deterministic and correct.\n";
            } else {
                std::cout << "\n❌ FAILURE: Trivial Encryption is NOT deterministic.\n";
            }
        }

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return -1;
    }
}
