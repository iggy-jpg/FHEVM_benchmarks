#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <string>
#include <array>

#include <tfhe++.hpp>
#include "HomPack.h"
#include "Autohog_Cop.hpp"
#include "SecureRandom.h"


struct ScopedTimer {
    std::chrono::high_resolution_clock::time_point t0;
    double* sink_ms;
    explicit ScopedTimer(double* out_ms) : t0(std::chrono::high_resolution_clock::now()), sink_ms(out_ms) {}
    ~ScopedTimer() {
        auto t1 = std::chrono::high_resolution_clock::now();
        *sink_ms += std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
    }
};

int main() try {
    using namespace TFHEpp;
    using P = lvl2param;

    autohog::SecureRng engine;
    std::bernoulli_distribution coin(0.5);

    // =========================
    
    // =========================
    double t_keygen_bk = 0.0;
    double t_keygen_iks = 0.0;
    double t_keygen_pk = 0.0;

    double t_encrypt_trlwe = 0.0;
    double t_save_trlwe = 0.0;
    double t_load_trlwe = 0.0;
    double t_decompress = 0.0;

    double t_pbs_total = 0.0;   
    double t_pack_avg = 0.0;    
    double t_save_pack = 0.0;
    double t_load_pack = 0.0;
    double t_verify_pack = 0.0;

    // =========================
    // KeyGen
    // =========================
    SecretKey sk;
    EvalKey ek;
    {
        ScopedTimer _(&t_keygen_bk);
        ek.emplacebkfft<lvl02param>(sk);
    }
    {
        ScopedTimer _(&t_keygen_iks);
        ek.emplaceiksk<lvl20param>(sk);
    }

    PKE::SeededRlwePkPair<P> spk;
    {
        ScopedTimer _(&t_keygen_pk);
        spk = PKE::generate_rlwe_pk_pair_seeded<P>(sk, engine);
    }

    
    const std::string pkfile = "rlwe_seeded_pk.bin";
    saveSeededRlwePublicKey<P>(spk, pkfile);
    printFileSize(pkfile);

    PKE::SeededRlwePkPair<P> spk2;
    loadSeededRlwePublicKey<P>(spk2, pkfile);

    // =========================
    
    // =========================
    constexpr uint32_t K = 10;
    static_assert(K <= P::n, "K must be <= P::n");

    std::array<uint8_t, K> m{};
    for (uint32_t i = 0; i < K; ++i) m[i] = static_cast<uint8_t>(coin(engine));

    std::array<bool, P::n> bits{};
    bits.fill(false);
    for (uint32_t i = 0; i < K; ++i) bits[i] = (m[i] != 0);

    TRLWE<P> ct_trlwe;
    {
        ScopedTimer _(&t_encrypt_trlwe);
        ct_trlwe = PKE::encrypt_bits_trlwe_with_seeded_rlwe_pair<P>(spk2, bits, engine);
    }

    const std::string ctfile = "cipher_seeded.trlwe";
    {
        ScopedTimer _(&t_save_trlwe);
        saveTRLWE<P>(ct_trlwe, ctfile);
    }
    printFileSize(ctfile);

    
    TRLWE<P> ct_loaded;
    {
        ScopedTimer _(&t_load_trlwe);
        loadTRLWE<P>(ct_loaded, ctfile);
    }

    std::vector<TLWE<P>> ct_vec;
    {
        ScopedTimer _(&t_decompress);
        ct_vec = TRLWE_Decompress<P>(ct_loaded, K);
    }

    // =========================
    
    // =========================
    std::vector<TLWE<P>> booted_vec(K);
    for (uint32_t i = 0; i < K; ++i) {
        auto t0 = std::chrono::high_resolution_clock::now();
        GateBootstrapping<lvl20param, lvl02param, P::μ>(booted_vec[i], ct_vec[i], ek);
        auto t1 = std::chrono::high_resolution_clock::now();
        t_pbs_total += std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;

        uint8_t dec_before = bootsSymDecrypt<P>({ct_vec[i]}, sk).at(0);
        uint8_t dec_after  = bootsSymDecrypt<P>({booted_vec[i]}, sk).at(0);

        std::cout << "PBS-CHK i=" << i
                  << " m=" << int(m[i])
                  << " dec(PKE)=" << int(dec_before)
                  << " dec(PBS)=" << int(dec_after)
                  << (m[i] != dec_before ? "  [WARN mismatch before PBS ❌]" : "")
                  << "\n";
    }

    // =========================
    
    // =========================
    AnnihilateKey<P> ahk;
    annihilatekeygen<P>(ahk, sk);

    constexpr uint32_t NUM_PACK_RUNS = 5;
    TRLWE<P> packed_after_pbs;
    {
        double accum_ms = 0.0;
        for (uint32_t r = 0; r < NUM_PACK_RUNS; ++r) {
            auto t0 = std::chrono::high_resolution_clock::now();
            HomPack::TLWE2TRLWEPacking<P>(packed_after_pbs, booted_vec, ahk);
            auto t1 = std::chrono::high_resolution_clock::now();
            accum_ms += std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
        }
        t_pack_avg = accum_ms / NUM_PACK_RUNS;
    }

    
    {
        ScopedTimer _(&t_verify_pack);
        const auto dec_poly = trlweSymDecrypt<P>(packed_after_pbs, sk.key.get<P>());
        bool ok = true;
        for (uint32_t i = 0; i < K; ++i) {
            uint8_t want = bootsSymDecrypt<P>({booted_vec[i]}, sk).at(0);
            uint8_t got  = static_cast<uint8_t>(dec_poly[i]);
            if (got != want) {
                ok = false;
                std::cout << "[PACK-CHK] idx=" << i
                          << " want=" << int(want)
                          << " got="  << int(got) << " ❌\n";
            }
        }
        std::cout << (ok ? "PACK-VERIFY (mem): PASSED ✅\n" : "PACK-VERIFY (mem): FAILED ❌\n");
    }

    
    const std::string packed_file = "cipher_after_pbs_packed.trlwe";
    {
        ScopedTimer _(&t_save_pack);
        saveTRLWE<P>(packed_after_pbs, packed_file);
    }
    printFileSize(packed_file);

    
    TRLWE<P> packed_loaded;
    {
        ScopedTimer _(&t_load_pack);
        loadTRLWE<P>(packed_loaded, packed_file);
    }
    {
        ScopedTimer _(&t_verify_pack); 
        const auto dec_poly = trlweSymDecrypt<P>(packed_loaded, sk.key.get<P>());
        bool ok = true;
        for (uint32_t i = 0; i < K; ++i) {
            uint8_t want = bootsSymDecrypt<P>({booted_vec[i]}, sk).at(0);
            uint8_t got  = static_cast<uint8_t>(dec_poly[i]);
            if (got != want) {
                ok = false;
                std::cout << "[PACK-CHK/RELOAD] idx=" << i
                          << " want=" << int(want)
                          << " got="  << int(got) << " ❌\n";
            }
        }
        std::cout << (ok ? "PACK-VERIFY (reload): PASSED ✅\n" : "PACK-VERIFY (reload): FAILED ❌\n");
    }

    auto print_bits = [](const char* tag, const std::array<uint8_t, K>& v) {
        std::cout << tag;
        for (uint8_t b : v) std::cout << int(b) << ' ';
        std::cout << '\n';
    };

    std::cout << "\n== Plaintext ==\n";
    print_bits("m (msb->lsb): ", m);

    std::cout << "\n== Timing (ms) ==\n";
    std::cout << "[KeyGen]\n";
    std::cout << "  BK (emplacebkfft): " << t_keygen_bk  << "\n";
    std::cout << "  IKS(emplaceiksk):  " << t_keygen_iks << "\n";
    std::cout << "  RLWE pk (seeded):  " << t_keygen_pk  << "\n";

    std::cout << "\n[Encrypt & IO]\n";
    std::cout << "  Encrypt TRLWE:     " << t_encrypt_trlwe << "\n";
    std::cout << "  Save TRLWE:        " << t_save_trlwe    << "\n";
    std::cout << "  Load TRLWE:        " << t_load_trlwe    << "\n";
    std::cout << "  Decompress (K):    " << t_decompress    << "\n";

    std::cout << "\n[PBS & Pack]\n";
    std::cout << "  PBS total (K):     " << t_pbs_total     << "\n";
    std::cout << "  TLWE->TRLWE avg:   " << t_pack_avg      << "\n";
    std::cout << "  Save packed:       " << t_save_pack     << "\n";
    std::cout << "  Load packed:       " << t_load_pack     << "\n";
    std::cout << "  Verify packed sum: " << t_verify_pack   << "  (mem+reload)\n";

    std::cout << "\nAll done.\n";
    return 0;

} catch (const std::exception& e) {
    std::cerr << "Fatal error: " << e.what() << "\n";
    return -1;
}
