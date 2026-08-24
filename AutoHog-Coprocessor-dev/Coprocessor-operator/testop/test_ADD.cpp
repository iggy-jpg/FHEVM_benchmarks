// test_add_pack.cpp
#include "Autohog_Cop.hpp"
#include "timer.h"

#include <tfhe++.hpp>
#include "HomPack.h"
#include "SecureRandom.h"

#include <random>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <array>
#include <stdexcept>
#include <cassert>
#include <thread>    // hardware_concurrency
#include <algorithm> // std::max
#include <numeric>   // std::accumulate
#include <boost/multiprecision/cpp_int.hpp>
#include "mem_probe.hpp"


// #if defined(__GLIBC__)
// #include <malloc.h>
// #endif

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace TFHEpp;
using namespace boost::multiprecision;


template <typename T>
inline void deep_free_vector(std::vector<T> &v)
{
    std::vector<T>().swap(v); 
}

template <typename T>
inline void free_object(T &obj)
{
    T tmp;
    using std::swap;
    swap(obj, tmp); 
}

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
    std::cerr << "Usage: " << prog << " <bit-width> [epochs] "
              << "[--no-epoch-detail|--show-epochs]\n";
}


static cpp_int normalize_twos_signed(const cpp_int &x, int W)
{
    cpp_int mod = cpp_int(1) << W;        // 2^W
    cpp_int half = cpp_int(1) << (W - 1); // 2^(W-1)
    cpp_int y = x % mod;
    if (y < 0)
        y += mod;
    if (y >= half)
        return y - mod; 
    return y;           
}

int main(int argc, char *argv[])
{
    using P = TFHEpp::lvl2param;              // TLWE/TRLWE level-2
    using BitsArray = std::array<bool, P::n>; 

    try
    {
        
        if (argc < 2)
        {
            print_usage(argv[0]);
            return 1;
        }
        int W_arg = std::atoi(argv[1]);
        if (W_arg <= 0)
        {
            std::cerr << "Error: bit-width must be positive.\n";
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
            {
                show_epoch_detail = false;
            }
            else if (s == "--show-epochs")
            {
                show_epoch_detail = true;
            }
        }

        
        // #if defined(__GLIBC__)
        // mallopt(M_ARENA_MAX, 2);
        // mallopt(M_TRIM_THRESHOLD, 64*1024);
        // #endif

#ifdef _OPENMP
        {
            int hw = std::max(1u, std::thread::hardware_concurrency());
            omp_set_num_threads(hw);
        }
#endif

        
        std::vector<double> t_encA, t_encB, t_deA, t_deB;
        std::vector<double> t_add_ser, t_add_par;
        std::vector<double> t_pack, t_save, t_load, t_verify;

        
        int ok_serial_vs_plain = 0;    
        int ok_parallel_vs_serial = 0; 
        int total_epochs = epochs;

        
        // MEM_SNAPSHOT("program start");

        
        double t_keygen_bk = 0.0;
        double t_keygen_iks = 0.0;
        double t_keygen_const = 0.0;
        double t_keygen_pk = 0.0;

        autohog::SecureRng eng;

        // [MOD]
        // MEM_SNAPSHOT("before keygen (once)");
        // MemScope key_scope("KeyGen (once)");

        SecretKey sk;
        EvalKey ek;
        {
            AutoTimer t_bk(&t_keygen_bk);
            ek.emplacebkfft<lvl02param>(sk);
        }
        {
            AutoTimer t_iks(&t_keygen_iks);
            ek.emplaceiksk<lvl20param>(sk);
        }
        {
            AutoTimer t_const(&t_keygen_const);
            init_constants(sk);
        }

        
        using PkPair = PKE::SeededRlwePkPair<P>;
        PkPair spk;
        {
            AutoTimer t_pk(&t_keygen_pk);
            spk = PKE::generate_rlwe_pk_pair_seeded<P>(sk, eng);
        }
        const std::string pkfile = "seeded_rlwe_public_key.bin";
        TFHEpp::saveSeededRlwePublicKey<P>(spk, pkfile);
        if (show_epoch_detail)
            TFHEpp::printFileSize(pkfile);

        PkPair spk2;
        TFHEpp::loadSeededRlwePublicKey<P>(spk2, pkfile);

        // [MOD]
        // MEM_SNAPSHOT("after keygen (once)");

        for (int ep = 0; ep < epochs; ++ep)
        {
            if (show_epoch_detail)
            {
                std::cout << "\n================= EPOCH " << (ep + 1) << "/" << epochs << " =================\n";
            }
            std::ostringstream ep_tag;
            ep_tag << "epoch " << (ep + 1);
            // [MOD]
            

            
            double t_enc_trlwe_A = 0.0;
            double t_enc_trlwe_B = 0.0;
            double t_decomp_A = 0.0;
            double t_decomp_B = 0.0;
            double t_add_ms = 0.0;
            double t_add_par_ms = 0.0;
            double t_pack_ms = 0.0;
            double t_save_pack = 0.0;
            double t_load_pack = 0.0;
            double t_verify_pack = 0.0;

            
            {
                
                int W = W_arg;
                if (W > static_cast<int>(P::n))
                {
                    if (show_epoch_detail)
                    {
                        std::cerr << "[Warn] given W=" << W
                                  << " > P::n=" << P::n
                                  << ", truncated.\n";
                    }
                    W = static_cast<int>(P::n);
                }

                
                std::vector<uint8_t> Apl(W), Bpl(W);
                std::uniform_int_distribution<int> bit01(0, 1);
                for (int i = 0; i < W; ++i)
                {
                    Apl[i] = bit01(eng);
                    Bpl[i] = bit01(eng);
                }

                // [MOD]
                // MEM_SNAPSHOT("before encrypt");

                
                BitsArray A_bits{}, B_bits{};
                {
                    auto A_arr = to_bit_array<P>(Apl);
                    auto B_arr = to_bit_array<P>(Bpl);
                    for (int i = 0; i < W; ++i)
                    {
                        A_bits[i] = A_arr[i];
                        B_bits[i] = B_arr[i];
                    }
                }

                
                TRLWE<P> ct_A, ct_B;
                {
                    AutoTimer t_encA(&t_enc_trlwe_A);
                    ct_A = PKE::encrypt_bits_trlwe_with_seeded_rlwe_pair<P>(spk2, A_bits, eng);
                }
                {
                    AutoTimer t_encB(&t_enc_trlwe_B);
                    ct_B = PKE::encrypt_bits_trlwe_with_seeded_rlwe_pair<P>(spk2, B_bits, eng);
                }


                const std::string ctfile_A = "cipher_A.trlwe";
                const std::string ctfile_B = "cipher_B.trlwe";
                TFHEpp::saveTRLWE<P>(ct_A, ctfile_A);
                TFHEpp::saveTRLWE<P>(ct_B, ctfile_B);
                if (show_epoch_detail)
                {
                    TFHEpp::printFileSize(ctfile_A);
                    TFHEpp::printFileSize(ctfile_B);
                }


                std::vector<TLWE<lvl2param>> A(W), B(W);
                {
                    AutoTimer t_deA(&t_decomp_A);
                    A = TRLWE_Decompress<P>(ct_A, W);
                }
                {
                    AutoTimer t_deB(&t_decomp_B);
                    B = TRLWE_Decompress<P>(ct_B, W);
                }



                HomAddResult res;           
                HomAddResult res_par_dummy; 

                {
                    AutoTimer timer(&t_add_ms);
                    res = HomAdd(A, B, W, ek, /*compute_overflow=*/true);
                }
                {
                    const int threads_hint =
                        static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
                    const std::vector<TLWE<lvl2param>> &A_cref = A;
                    const std::vector<TLWE<lvl2param>> &B_cref = B;
                    AutoTimer timer(&t_add_par_ms);
                    res_par_dummy = HomAdd_pal(A_cref, B_cref, W, ek,
                                               /*compute_overflow=*/true, threads_hint);
                }

                // [MOD]
                // MEM_SNAPSHOT("after homadd");

                
                auto dec_bit = [&sk](const TLWE<lvl2param> &ct) -> uint8_t
                {
                    auto v = TFHEpp::bootsSymDecrypt<lvl2param>(std::vector<TLWE<lvl2param>>{ct}, sk);
                    return v.at(0);
                };


                std::vector<uint8_t> Sum_bits(W);
                for (int i = 0; i < W; ++i)
                    Sum_bits[i] = dec_bit(res.SUM[i]);

                if (show_epoch_detail)
                {
                    auto print_bits = [](const char *tag, const std::vector<uint8_t> &v)
                    {
                        std::cout << tag;
                        for (uint8_t b : v)
                            std::cout << int(b) << ' ';
                        std::cout << '\n';
                    };

                    std::cout << "\n== Plaintext & Decrypted Result ==\n";
                    print_bits("A   (msb->lsb): ", Apl);
                    print_bits("B   (msb->lsb): ", Bpl);
                    print_bits("Sum (msb->lsb): ", Sum_bits);
                }

                cpp_int A_s = bits_to_int_signed(Apl, W);
                cpp_int B_s = bits_to_int_signed(Bpl, W);
                cpp_int S_s = bits_to_int_signed(Sum_bits, W);
                cpp_int plain_add = A_s + B_s;
                cpp_int plain_add_norm = normalize_twos_signed(plain_add, W);

                if (show_epoch_detail)
                {
                    std::cout << "-----------------------------\n";
                    std::cout << "A (signed)     = " << A_s << "\n";
                    std::cout << "B (signed)     = " << B_s << "\n";
                    std::cout << "A + B (signed) = " << plain_add << "\n";
                    std::cout << "Sum (signed)   = " << S_s << "\n";
                }

                int ovf_ser = -1;
                if (res.ovf_ct)
                {
                    ovf_ser = dec_bit(*res.ovf_ct);
                    if (show_epoch_detail)
                    {
                        std::cout << "Overflow bit   = " << ovf_ser << "\n";
                        if ((plain_add_norm != S_s) && ovf_ser == 0)
                            std::cout << "Warning❌: Result mismatch (no overflow flagged)!\n";
                    }
                }

                
                std::vector<uint8_t> Sum_bits_par(W);
                for (int i = 0; i < W; ++i)
                    Sum_bits_par[i] = dec_bit(res_par_dummy.SUM[i]);
                cpp_int S_s_par = bits_to_int_signed(Sum_bits_par, W);
                int ovf_par = -1;
                if (res_par_dummy.ovf_ct)
                    ovf_par = dec_bit(*res_par_dummy.ovf_ct);

                
                bool pass_serial_vs_plain = (S_s == plain_add_norm);

                
                bool pass_parallel_vs_serial = true;
                int first_mis_idx = -1;
                for (int i = 0; i < W; ++i)
                {
                    if (Sum_bits[i] != Sum_bits_par[i])
                    {
                        pass_parallel_vs_serial = false;
                        first_mis_idx = i;
                        break;
                    }
                }
                if (pass_parallel_vs_serial)
                {
                    if ((ovf_ser != -1 || ovf_par != -1) && (ovf_ser != ovf_par))
                        pass_parallel_vs_serial = false;
                }

                if (pass_serial_vs_plain)
                    ++ok_serial_vs_plain;
                if (pass_parallel_vs_serial)
                    ++ok_parallel_vs_serial;

                if (show_epoch_detail)
                {
                    std::cout << "\n== (Epoch " << (ep + 1) << ") Correctness Checks ==\n";
                    std::cout << "[Serial vs Plain 2^W]  : " << (pass_serial_vs_plain ? "PASS ✅" : "FAIL ❌") << "\n";
                    std::cout << "[Parallel vs Serial]    : " << (pass_parallel_vs_serial ? "PASS ✅" : "FAIL ❌") << "\n";
                    if (!pass_parallel_vs_serial)
                    {
                        std::cout << "  First mismatch index  : " << first_mis_idx << "\n";
                        std::cout << "  ovf(serial)=" << ovf_ser << ", ovf(par)=" << ovf_par << "\n";
                    }
                }

                // [MOD]
                // MEM_SNAPSHOT("before pack+io");

                
                AnnihilateKey<P> ahk;
                annihilatekeygen<P>(ahk, sk);

                TRLWE<P> packed_sum;
                {
                    AutoTimer t_pack_timer(&t_pack_ms);
                    HomPack::TLWE2TRLWEPacking<P>(packed_sum, res.SUM, ahk);
                }

                
                const std::string packed_file = "sum_packed.trlwe";
                {
                    AutoTimer t_sv(&t_save_pack);
                    saveTRLWE<P>(packed_sum, packed_file);
                }
                if (show_epoch_detail)
                    TFHEpp::printFileSize(packed_file);

                TRLWE<P> packed_loaded;
                {
                    AutoTimer t_ld(&t_load_pack);
                    loadTRLWE<P>(packed_loaded, packed_file);
                }
                {
                    AutoTimer t_vf(&t_verify_pack);
                    const auto dec_poly = trlweSymDecrypt<P>(packed_loaded, sk.key.get<P>());
                    bool ok = true;
                    for (int i = 0; i < W; ++i)
                    {
                        const uint8_t want = Sum_bits[i];
                        const uint8_t got = static_cast<uint8_t>(dec_poly[i]);
                        if (want != got)
                        {
                            ok = false;
                            if (show_epoch_detail)
                            {
                                std::cout << "[PACK-CHK] idx=" << i
                                          << " want=" << int(want)
                                          << " got=" << int(got) << " ❌\n";
                            }
                        }
                    }
                    if (show_epoch_detail)
                    {
                        std::cout << (ok ? "PACK-VERIFY: PASSED ✅\n" : "PACK-VERIFY: FAILED ❌\n");
                    }
                }

                // [MOD]
                // MEM_SNAPSHOT("after pack+io");

                
                std::remove(ctfile_A.c_str());
                std::remove(ctfile_B.c_str());
                std::remove(packed_file.c_str());

                
                if (show_epoch_detail)
                {
                    std::cout << "\n== (Epoch " << (ep + 1) << ") Encrypt+Pack Timing (ms) ==\n";
                    std::cout << "encrypt TRLWE(A):   " << t_enc_trlwe_A << " ms\n";
                    std::cout << "encrypt TRLWE(B):   " << t_enc_trlwe_B << " ms\n";
                    std::cout << "TRLWE_Decompress(A):" << t_decomp_A << " ms\n";
                    std::cout << "TRLWE_Decompress(B):" << t_decomp_B << " ms\n";

                    std::cout << "\n== (Epoch " << (ep + 1) << ") HomAdd Timing (ms) ==\n";
                    std::cout << "HomAdd (serial):    " << t_add_ms << " ms\n";
                    std::cout << "HomAdd (parallel):  " << t_add_par_ms << " ms\n";

                    std::cout << "\n== (Epoch " << (ep + 1) << ") Post-Compute Pack+IO Timing (ms) ==\n";
                    std::cout << "TLWE->TRLWE Pack:   " << t_pack_ms << " ms\n";
                    std::cout << "Save packed TRLWE:  " << t_save_pack << " ms\n";
                    std::cout << "Load packed TRLWE:  " << t_load_pack << " ms\n";
                    std::cout << "Verify packed:      " << t_verify_pack << " ms\n";
                }

                
                deep_free_vector(A);
                deep_free_vector(B);
                deep_free_vector(res.SUM);
                if (res.ovf_ct)
                    res.ovf_ct.reset();
                deep_free_vector(res_par_dummy.SUM);
                if (res_par_dummy.ovf_ct)
                    res_par_dummy.ovf_ct.reset();
                free_object(ct_A);
                free_object(ct_B);
                free_object(packed_sum);
                free_object(packed_loaded);
                free_object(ahk);

                
                // #if defined(__GLIBC__)
                // malloc_trim(0);
                // #endif
            } 

            
            t_encA.push_back(t_enc_trlwe_A);
            t_encB.push_back(t_enc_trlwe_B);
            t_deA.push_back(t_decomp_A);
            t_deB.push_back(t_decomp_B);

            t_add_ser.push_back(t_add_ms);
            t_add_par.push_back(t_add_par_ms); 

            t_pack.push_back(t_pack_ms);
            t_save.push_back(t_save_pack);
            t_load.push_back(t_load_pack);
            t_verify.push_back(t_verify_pack);
        }

        
        std::cout << "\n== KeyGen (once) Timing (ms) ==\n";
        std::cout << "  BK (emplacebkfft): " << t_keygen_bk << " ms\n";
        std::cout << "  IKS(emplaceiksk):  " << t_keygen_iks << " ms\n";
        std::cout << "  init_constants:    " << t_keygen_const << " ms\n";
        std::cout << "  RLWE pk pair:      " << t_keygen_pk << " ms\n";

        
        if (epochs > 1)
        {
            std::cout << "\n================= SUMMARY (average over " << epochs << " epochs) =================\n";

            std::cout << "\n== Encrypt+Pack Timing AVG (ms) ==\n";
            std::cout << "encrypt TRLWE(A):   " << avg(t_encA) << " ms\n";
            std::cout << "encrypt TRLWE(B):   " << avg(t_encB) << " ms\n";
            std::cout << "TRLWE_Decompress(A):" << avg(t_deA) << " ms\n";
            std::cout << "TRLWE_Decompress(B):" << avg(t_deB) << " ms\n";

            std::cout << "\n== HomAdd Timing AVG (ms) ==\n";
            std::cout << "HomAdd (serial):    " << avg(t_add_ser) << " ms\n";
            std::cout << "HomAdd (parallel):  " << avg(t_add_par) << " ms\n";

            std::cout << "\n== Post-Compute Pack+IO Timing AVG (ms) ==\n";
            std::cout << "TLWE->TRLWE Pack:   " << avg(t_pack) << " ms\n";
            std::cout << "Save packed TRLWE:  " << avg(t_save) << " ms\n";
            std::cout << "Load packed TRLWE:  " << avg(t_load) << " ms\n";
            std::cout << "Verify packed:      " << avg(t_verify) << " ms\n";

            
            auto percent = [](int num, int den) -> double
            {
                if (den <= 0)
                    return 0.0;
                return 100.0 * static_cast<double>(num) / static_cast<double>(den);
            };
            std::cout << "\n== Correctness Summary ==\n";
            std::cout << "Serial vs Plain (2^W): " << ok_serial_vs_plain << "/" << total_epochs
                      << " (" << percent(ok_serial_vs_plain, total_epochs) << "%)\n";
            std::cout << "Parallel vs Serial   : " << ok_parallel_vs_serial << "/" << total_epochs
                      << " (" << percent(ok_parallel_vs_serial, total_epochs) << "%)\n";
        }

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return -1;
    }
}
