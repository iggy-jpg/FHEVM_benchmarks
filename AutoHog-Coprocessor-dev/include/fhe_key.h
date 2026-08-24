#pragma once
#include <tfhe++.hpp>
#include "HomOP.h"
#include "HomUtil.h"
#include "SaveFile.h"

// Note: HomPKE.h defines the global PKE namespace
#include "HomPKE.h"         
#include "Autohog_Cop.hpp"
#include "HomPack.h"        // TRLWE_Decompress / TLWE2TRLWEPacking / annihilatekeygen
#include "HomPKE.h"         // RLWE Public Key Encryption Interface
#include "ZKPKE.h"
#include "HomDecompress.h"
#include <memory>
#include <thread>
#include <vector>
#include <random> 
#include <iostream> 
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <type_traits>

#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>

using namespace TFHEpp;
using std::vector;
using std::string;

#define USE_PAL 1

// ==========================================
// [FIX] Cereal Serialization Support
// ==========================================
// Must be injected into the ::PKE namespace where SeededRlwePkPair is defined
namespace PKE {
    template <class Archive, class P>
    void serialize(Archive & ar, SeededRlwePkPair<P> & spk)
    {
        ar(spk.seed, spk.b);
    }
}
// ==========================================

using tlwe_vec_t = std::vector<TLWE<lvl2param>>;
using tlwe_t = TLWE<lvl2param>;
using trlwe_t = TRLWE<lvl2param>;
using Bytes = vector<uint8_t>;
using BitsArray = std::array<bool, lvl2param::n>;

using trlwe_v2_t = TRLWE<lvl1param>;
using BitsArrayV2 = std::array<bool, lvl1param::n>;
using crs_lvl = lvl2param;

struct TLWEWrapper {
    tlwe_t data;
    TLWEWrapper(const tlwe_t& d) { data = d; }
    TLWEWrapper() = default; 
};

template<typename param = lvl2param>
struct TRLWEWrapperGeneric {
    TRLWE<param> data;
    size_t W;
    TRLWEWrapperGeneric(const TRLWE<param>& d, size_t w) {
        data = d;
        W = w;
    }
};

using TRLWEWrapper = TRLWEWrapperGeneric<lvl2param>;
using TRLWEV2Wrapper = TRLWEWrapperGeneric<lvl1param>;

struct TLWEVecWrapper {
    tlwe_vec_t data;
    size_t W;
    TLWEVecWrapper(const tlwe_vec_t& d) {
        data = d;
        W = d.size();
    }
};

using tlwe = TLWEWrapper;
using trlwe = TRLWEWrapper;
using trlwe_v2 = TRLWEV2Wrapper;
using tlwe_vec = TLWEVecWrapper;

using tlwe_ptr = std::unique_ptr<TLWEWrapper>;
using trlwe_ptr = std::unique_ptr<TRLWEWrapper>;
using trlwe_v2_ptr = std::unique_ptr<TRLWEV2Wrapper>;
using tlwe_vec_ptr = std::unique_ptr<TLWEVecWrapper>;

inline void require_valid_tlwe_vector(const tlwe_vec_ptr& value,
                                      const char* name) {
    if (!value || value->data.empty() || value->W != value->data.size() ||
        value->W > lvl2param::n) {
        throw std::invalid_argument(std::string("invalid ") + name);
    }
}

inline void require_compatible_tlwe_vectors(const tlwe_vec_ptr& lhs,
                                            const tlwe_vec_ptr& rhs,
                                            const char* operation) {
    require_valid_tlwe_vector(lhs, "left ciphertext operand");
    require_valid_tlwe_vector(rhs, "right ciphertext operand");
    if (lhs->W != rhs->W) {
        throw std::invalid_argument(std::string(operation) +
                                    " requires equal ciphertext widths");
    }
}

template<typename P>
inline void require_valid_trlwe_wrapper(
    const std::unique_ptr<TRLWEWrapperGeneric<P>>& value,
    const char* name) {
    if (!value || value->W == 0 || value->W > P::n) {
        throw std::invalid_argument(std::string("invalid ") + name);
    }
}

inline tlwe_ptr to_tlwe_ptr(const tlwe_vec_ptr& data) {
    require_valid_tlwe_vector(data, "TLWE vector");
    return tlwe_ptr(new tlwe(data->data[0]));
}
inline tlwe_vec_ptr to_tlwe_vec_ptr(const tlwe_ptr& data) {
    if (!data) throw std::invalid_argument("invalid TLWE ciphertext");
    tlwe_vec_t d = {data->data};
    return tlwe_vec_ptr(new tlwe_vec(d));
}

template<typename T>
Bytes to_vec(T t) {
    static_assert(std::is_integral_v<T>, "to_vec requires an integral type");
    using U = std::make_unsigned_t<T>;
    U value = static_cast<U>(t);
    Bytes result;
    for (size_t i = 0; i < sizeof(T); i++) {
        uint8_t byte = static_cast<uint8_t>(
            value >> ((sizeof(T) - 1 - i) * 8));
        result.push_back(byte);
    }
    return result;
}

template<typename T>
T from_vec(const Bytes& data) {
    static_assert(std::is_integral_v<T>, "from_vec requires an integral type");
    if (data.size() != sizeof(T)) {
        throw std::invalid_argument("integer byte vector has the wrong length");
    }
    using U = std::make_unsigned_t<T>;
    U t = 0;
    for (size_t i = 0; i < sizeof(T); i++) {
        t |= static_cast<U>(data[i]) << ((sizeof(T) - 1 - i) * 8);
    }
    return static_cast<T>(t);
}

Bytes serialize_ciphertext_with_proof(const trlwe_ptr& ciphertext, const ZKProof::Proof& proof);
bool deserialize_ciphertext_with_proof(trlwe_ptr& ciphertext, ZKProof::Proof& proof, const Bytes& data);

Bytes serialize_ciphertexts(const tlwe_vec_ptr& d);
Bytes serialize_compressed_ciphertexts(const trlwe_ptr& d);
Bytes serialize_compressed_ciphertexts_v2(const trlwe_v2_ptr& d);
void save_ciphertexts(const string& file, const tlwe_vec_ptr& d);
void save_compressed_ciphertexts(const string& file, const trlwe_ptr& d);
void save_compressed_ciphertexts_v2(const string& file, const trlwe_v2_ptr& d);

tlwe_vec_ptr deserialize_ciphertexts(const Bytes& bytes);
trlwe_ptr deserialize_compressed_ciphertexts(const Bytes& bytes);
trlwe_v2_ptr deserialize_compressed_ciphertexts_v2(const Bytes& bytes);
tlwe_vec_ptr load_ciphertexts(const string& file);
trlwe_ptr load_compressed_ciphertexts(const string& file);
trlwe_v2_ptr load_compressed_ciphertexts_v2(const string& file);

inline tlwe_vec_ptr decompress_ciphertexts(const trlwe_ptr& d) {
    require_valid_trlwe_wrapper<lvl2param>(d, "compressed ciphertext");
    tlwe_vec vec = TRLWE_Decompress<lvl2param>(d->data, d->W);
    return tlwe_vec_ptr(new tlwe_vec(vec));
}

typedef std::function<tlwe_t(uint8_t)> EncryptFn;
typedef std::function<uint8_t(tlwe_t)> DecryptFn;

// [FIX] Use ::PKE and TFHEpp:: explicitly to avoid ambiguity
void load_pk(::PKE::SeededRlwePkPair<TFHEpp::lvl2param> *spk, const char* pk, int len);
void load_pk_v2(::PKE::SeededRlwePkPair<TFHEpp::lvl1param> *spk, const char* pk, int len);
void load_crs(::ZKProof::PublicParams<crs_lvl> *zk_params, const char* crs, int len);
void load_crs(::ZKProof::PublicParams<crs_lvl> *zk_params, const char* crs, int len,
              const std::string& expected_sha256);

class FheKey {
private:
    SecretKey sk;
    EvalKey ek;
    AnnihilateKey<lvl2param> ahk2;
    
    // [FIX] Use ::PKE
    ::PKE::SeededRlwePkPair<TFHEpp::lvl1param> spk_lvl1;
    ::PKE::SeededRlwePkPair<TFHEpp::lvl2param> spk_lvl2;
    ::ZKProof::PublicParams<crs_lvl> zk_params;
    int cores = 1;

    // Helper: Serialize any object to string for comparison
    template <class T>
    bool is_equal_serialized(const T& a, const T& b) const {
        std::stringstream ss_a, ss_b;
        {
            cereal::BinaryOutputArchive oarchive_a(ss_a);
            oarchive_a(a);
        }
        {
            cereal::BinaryOutputArchive oarchive_b(ss_b);
            oarchive_b(b);
        }
        return ss_a.str() == ss_b.str();
    }

public:
    FheKey() {
        int cores = std::thread::hardware_concurrency();
        if (cores > 0) {
            this->cores = cores;
        }
    }

    bool has_zk_crs() const {
        return ZKPKE::has_loaded_crs(zk_params);
    }

    void require_zk_crs(const char* context) const {
        if (!has_zk_crs()) {
            throw std::runtime_error(std::string("missing CRS, cannot ") + context);
        }
    }

    // ==========================================
    // Equality Comparison Operator
    // ==========================================
    bool operator==(const FheKey& other) const {
        bool is_equal = true;

        // 1. Compare SecretKey
        if (!is_equal_serialized(sk, other.sk)) {
            std::cout << "❌ [Mismatch] SecretKey differs!" << std::endl;
            is_equal = false;
        }

        // 2. Compare EvalKey
        if (!is_equal_serialized(ek, other.ek)) {
            std::cout << "❌ [Mismatch] EvalKey differs!" << std::endl;
            is_equal = false;
        }

        // 3. Compare AnnihilateKey
        if (!is_equal_serialized(ahk2, other.ahk2)) {
            std::cout << "❌ [Mismatch] AnnihilateKey differs!" << std::endl;
            is_equal = false;
        }

        // 4. Compare Public Keys (Level 1)
        if (!is_equal_serialized(spk_lvl1, other.spk_lvl1)) {
            std::cout << "❌ [Mismatch] SeededRlwePkPair (Level 1) differs!" << std::endl;
            is_equal = false;
        }

        // 5. Compare Public Keys (Level 2)
        if (!is_equal_serialized(spk_lvl2, other.spk_lvl2)) {
            std::cout << "❌ [Mismatch] SeededRlwePkPair (Level 2) differs!" << std::endl;
            is_equal = false;
        }

        return is_equal;
    }

    bool operator!=(const FheKey& other) const {
        return !(*this == other);
    }

    void keygen(uint8_t ek_type) {
        auto gen_key = [&](){
            std::exception_ptr generation_error;
            std::mutex generation_error_mutex;
            auto launch_checked = [&](auto fn) {
                return std::thread([&, fn = std::move(fn)]() mutable {
                    try {
                        fn();
                    } catch (...) {
                        std::lock_guard<std::mutex> lock(generation_error_mutex);
                        if (!generation_error) generation_error = std::current_exception();
                    }
                });
            };

            auto gen_ek = [&]() {
                ek.emplacebkfft<lvl02param>(sk); 
                if (ek_type == 2 || ek_type == 4) {
                    ek.emplacebkfft<lvl12param>(sk); 
                }
                ek.emplaceiksk<lvl20param>(sk);
                if (ek_type == 1 || ek_type == 4) {
                    ek.emplaceiksk<lvl10param>(sk);  
                }
                if (ek_type == 3 || ek_type == 4) {
                    ek.emplaceiksk<lvl12param>(sk); 
                }
                std::cout << "gen ek ok" << std::endl;
            };
    
            auto gen_spk_lvl1 = [&]() {
                // [FIX] Use ::PKE
                autohog::SecureRng eng;
                spk_lvl1 = ::PKE::generate_rlwe_pk_pair_seeded<TFHEpp::lvl1param>(
                    sk, eng);
                std::cout << "gen spk lvl1 ok" << std::endl;
            };

            auto gen_spk_lvl2 = [&]() {
                autohog::SecureRng eng;
                spk_lvl2 = ::PKE::generate_rlwe_pk_pair_seeded<TFHEpp::lvl2param>(
                    sk, eng);
                init_constants(spk_lvl2);
                std::cout << "gen spk lvl2 ok" << std::endl;
            };

            auto gen_ahk = [&]() {
                annihilatekeygen<lvl2param>(ahk2, sk);
                std::cout << "gen ahk ok" << std::endl;
            };

            auto ek_thread = launch_checked(gen_ek);
            auto spk_lvl1_thread = launch_checked(gen_spk_lvl1);
            auto spk_lvl2_thread = launch_checked(gen_spk_lvl2);
            auto ahk_thread = launch_checked(gen_ahk);
            ek_thread.join();
            spk_lvl1_thread.join();
            spk_lvl2_thread.join();
            ahk_thread.join();
            if (generation_error) std::rethrow_exception(generation_error);
        };

        gen_key();

        // A proving CRS must not be created by the online key owner.  Clear any
        // previously loaded CRS and require callers to provision an externally
        // generated, authenticated ceremony artifact explicitly.
        zk_params = ::ZKProof::PublicParams<crs_lvl>();
        std::cout << "CRS not generated; load an authenticated external CRS before ZK use"
                  << std::endl;
    }

    void load(const string& eval_key_file, 
              const string& annihilate_key_file, 
              const string& public_key_file, 
              const string& public_key_v2_file,
              const string& crs_file,
              const string& crs_expected_sha256) {
        if (!crs_file.empty() && crs_expected_sha256.empty()) {
            throw std::invalid_argument(
                "A trusted CRS SHA-256 digest is required when loading a CRS file"
            );
        }

        std::exception_ptr load_error;
        std::mutex load_error_mutex;
        auto launch_checked = [&](auto fn) {
            return std::thread([&, fn = std::move(fn)]() mutable {
                try {
                    fn();
                } catch (...) {
                    std::lock_guard<std::mutex> lock(load_error_mutex);
                    if (!load_error) load_error = std::current_exception();
                }
            });
        };

        auto ek_thread = launch_checked([&]() {
            if (!loadEvalKey(ek, eval_key_file)) {
                throw std::runtime_error("Failed to load evaluation key");
            }
            std::cout << "load ek ok" << std::endl;
        });
        
        auto pk_lvl1_thread = launch_checked([&]() {
            if (!loadSeededRlwePublicKey<lvl1param>(spk_lvl1, public_key_v2_file)) {
                throw std::runtime_error("Failed to load lvl1 seeded RLWE public key");
            }
            std::cout << "load pk lvl1 ok" << std::endl;
        });
        
        auto pk_lvl2_thread = launch_checked([&]() {
            if (!loadSeededRlwePublicKey<lvl2param>(spk_lvl2, public_key_file)) {
                throw std::runtime_error("Failed to load lvl2 seeded RLWE public key");
            }
            std::cout << "load spk lvl2 ok" << std::endl;
        });
        //init_constants(spk_lvl2);

        auto ahk_thread = launch_checked([&]() {
            if (!loadAnnihilateKey<lvl2param>(ahk2, annihilate_key_file)) {
                throw std::runtime_error("Failed to load annihilation key");
            }
            std::cout << "load ahk ok" << std::endl;
        });

        auto crs_thread = launch_checked([&]() {
            if (!crs_file.empty()) {
                zk_params = ZKProof::load_public_params_by_expected_hash<crs_lvl>(
                    crs_file,
                    crs_expected_sha256
                );
            } else {
                zk_params = ::ZKProof::PublicParams<crs_lvl>();
            }
            std::cout << "load crs ok" << std::endl;
        });
        ek_thread.join();
        pk_lvl1_thread.join();
        pk_lvl2_thread.join();
        ahk_thread.join();
        crs_thread.join();
        if (load_error) std::rethrow_exception(load_error);
    }

    void load(const string& eval_key_file, 
              const string& annihilate_key_file, 
              const string& public_key_file, 
              const string& public_key_v2_file,
              const string& crs_file) {
        if (!crs_file.empty()) {
            throw std::invalid_argument(
                "Loading a CRS requires an independently trusted SHA-256 digest; "
                "use the six-argument FheKey::load overload"
            );
        }
        load(eval_key_file, annihilate_key_file, public_key_file,
             public_key_v2_file, "", "");
    }

    void load(const string& eval_key_file, 
              const string& annihilate_key_file, 
              const string& public_key_file, 
              const string& public_key_v2_file) {
        load(eval_key_file, annihilate_key_file, public_key_file, public_key_v2_file, "");
    }

    void load_sk(const string& secret_key_file) {
        // Secret Key is optional if we are a pure server and have all keys loaded
        if (!loadSecretKey(sk, secret_key_file)) {
            throw std::runtime_error("Failed to load secret key");
        }
    }

    void load_sk_from_memory(const vector<uint8_t>& data) {
        string data2(data.begin(), data.end());
        if (!loadSecretKeyFromMemory(sk, data2)) {
            throw std::runtime_error("Failed to load secret key from memory");
        }
    }

    void refresh(const string& secret_key_file, uint8_t ek_type) {
        if (!loadSecretKey(sk, secret_key_file)) {
            throw std::runtime_error("Failed to refresh from secret key");
        }
        keygen(ek_type);
    }

    void load_pk(const char* pk, int len) {
        ::load_pk(&spk_lvl2, pk, len);
        init_constants(spk_lvl2);
    }

    void load_pk_v2(const char* pk, int len) {
        ::load_pk_v2(&spk_lvl1, pk, len);
    }

    void load_crs(const char* crs, int len) {
        (void)crs;
        (void)len;
        throw std::invalid_argument(
            "Loading CRS bytes requires an independently trusted SHA-256 digest"
        );
    }

    void load_crs(const char* crs, int len, const string& expected_sha256) {
        if (crs == nullptr || len <= 0 || expected_sha256.empty()) {
            throw std::invalid_argument(
                "CRS bytes, positive length, and trusted SHA-256 digest are required"
            );
        }
        zk_params = ZKProof::load_public_params_from_memory_by_expected_hash<crs_lvl>(
            crs,
            len,
            expected_sha256
        );
    }

    Bytes save_eval_key() const {
        string str;
        if (!saveEvalKeyToMemory(ek, str)) {
            throw std::runtime_error("Failed to serialize evaluation key");
        }
        return Bytes(str.begin(), str.end());
    }

    Bytes save_secret_key() const {
        string str;
        if (!saveSecretKeyToMemory(sk, str)) {
            throw std::runtime_error("Failed to serialize secret key");
        }
        return Bytes(str.begin(), str.end());
    }

    Bytes save_annihilate_key() const {
        string str;
        if (!saveAnnihilateKeyToMemory<lvl2param>(ahk2, str)) {
            throw std::runtime_error("Failed to serialize annihilation key");
        }
        return Bytes(str.begin(), str.end());
    }

    Bytes save_public_key() const {
        string str;
        if (!saveSeededRlwePublicKeyToMemory<lvl2param>(spk_lvl2, str)) {
            throw std::runtime_error("Failed to serialize level-2 public key");
        }
        return Bytes(str.begin(), str.end());
    }

    Bytes save_public_key_v2() const {
        string str;
        if (!saveSeededRlwePublicKeyToMemory<lvl1param>(spk_lvl1, str)) {
            throw std::runtime_error("Failed to serialize level-1 public key");
        }
        return Bytes(str.begin(), str.end());
    }

    Bytes save_crs() const {
        require_zk_crs("save CRS to memory");
        string str;
        ZKProof::save_public_params_to_memory<crs_lvl>(zk_params, str);
        return Bytes(str.begin(), str.end());
    }
    void save(const string& secret_key_file, 
              const string& eval_key_file, 
              const string& annihilate_key_file, 
              const string& public_key_file, 
              const string& public_key_v2_file,
              const string& crs_file) const {
        std::exception_ptr save_error;
        std::mutex save_error_mutex;
        auto launch_checked = [&](auto fn) {
            return std::thread([&, fn = std::move(fn)]() mutable {
                try {
                    fn();
                } catch (...) {
                    std::lock_guard<std::mutex> lock(save_error_mutex);
                    if (!save_error) save_error = std::current_exception();
                }
            });
        };

        auto ek_thread = launch_checked([&]() {
            if (!saveEvalKey(ek, eval_key_file)) {
                throw std::runtime_error("Failed to save evaluation key");
            }
        });

        auto sk_thread = launch_checked([&]() {
            if (!saveSecretKey(sk, secret_key_file)) {
                throw std::runtime_error("Failed to save secret key");
            }
        });
        
        auto ahk_thread = launch_checked([&]() {
            // Save AnnihilateKey
            if (!annihilate_key_file.empty() &&
                !saveAnnihilateKey<lvl2param>(ahk2, annihilate_key_file)) {
                throw std::runtime_error("Failed to save annihilation key");
            }
        });
        
        auto pkv2_thread = launch_checked([&]() {
            if (!saveSeededRlwePublicKey<lvl1param>(spk_lvl1, public_key_v2_file)) {
                throw std::runtime_error("Failed to save level-1 public key");
            }
        });

        auto pk_thread = launch_checked([&]() {
            if (!saveSeededRlwePublicKey<lvl2param>(spk_lvl2, public_key_file)) {
                throw std::runtime_error("Failed to save level-2 public key");
            }
        });

        auto crs_thread = launch_checked([&]() {
            if (!crs_file.empty()) {
                require_zk_crs("save CRS");
                ZKProof::save_public_params<crs_lvl>(zk_params, crs_file);
            }
        });

        ek_thread.join();
        sk_thread.join();
        ahk_thread.join();
        pkv2_thread.join();
        pk_thread.join();
        crs_thread.join();
        if (save_error) std::rethrow_exception(save_error);
    }

    void save(const string& secret_key_file, 
              const string& eval_key_file, 
              const string& annihilate_key_file, 
              const string& public_key_file, 
              const string& public_key_v2_file) const {
        save(secret_key_file, eval_key_file, annihilate_key_file, public_key_file, public_key_v2_file, "");
    }

    // Overloaded save for backward compatibility
    void save(const string& secret_key_file, const string& eval_key_file, const string& public_key_file, const string& public_key_v2_file) const {
        if (!saveEvalKey(ek, eval_key_file) ||
            !saveSecretKey(sk, secret_key_file) ||
            !saveSeededRlwePublicKey<lvl1param>(spk_lvl1, public_key_v2_file) ||
            !saveSeededRlwePublicKey<lvl2param>(spk_lvl2, public_key_file)) {
            throw std::runtime_error("Failed to save one or more key files");
        }
        // saveAnnihilateKey skipped
    }

    trlwe_ptr compress(const tlwe_vec_ptr& ct) const {
        require_valid_tlwe_vector(ct, "TLWE vector to compress");
        trlwe_t packed_out;
        HomPack::TLWE2TRLWEPacking<lvl2param>(packed_out, ct->data, ahk2);
        return trlwe_ptr(new trlwe(packed_out, ct->W));
    }

    trlwe_ptr encrypt(uint8_t b) const {
        Bytes data = {b};
        return encryptData(data);
    }

    trlwe_ptr trivialEncrypt(uint8_t b) const {
        Bytes data = {b};
        return trivialEncryptData(data);
    }

    trlwe_v2_ptr encryptV2(uint8_t b) const {
        Bytes data = {b};
        return encryptDataV2(data);
    }

    Bytes toBitVec(const Bytes& data) const {
        if (data.size() > std::numeric_limits<size_t>::max() / 8) {
            throw std::length_error("plaintext length overflows bit count");
        }
        size_t W = data.size() * 8;
        Bytes A(W);
        for (size_t i = 0; i < W; i += 8) {
            for (size_t j = 0; j < 8; j++) {
                uint8_t tmp = (data[i / 8] >> (8 - 1 - j)) & 0x01;
                A[i + j] = tmp;
            }
        }
        return A;
    }

    trlwe_ptr encryptData(const Bytes& boolArray) const {
        if (boolArray.empty() || boolArray.size() > lvl2param::n) {
            throw std::invalid_argument("level-2 plaintext must contain 1..n bits");
        }
        autohog::SecureRng eng;
        BitsArray bits = to_bit_array<lvl2param>(boolArray);
        // [FIX] Use ::PKE
        trlwe_t ct = ::PKE::encrypt_bits_trlwe_with_seeded_rlwe_pair<lvl2param>(spk_lvl2, bits, eng);
        return trlwe_ptr(new trlwe(ct, boolArray.size()));
    }

    trlwe_ptr trivialEncryptData(const Bytes& boolArray) const {
        if (boolArray.empty() || boolArray.size() > lvl2param::n) {
            throw std::invalid_argument("level-2 plaintext must contain 1..n bits");
        }
        tlwe_vec_t vec;
        vec.resize(boolArray.size());
        for (size_t i = 0; i < boolArray.size(); i++) {
              tlwe_t res = {}; // Zero init
              res[lvl2param::n] = (boolArray[i] ? lvl2param::μ : -lvl2param::μ);
              vec[i] = res;
        }
        tlwe_vec_ptr ptr(new tlwe_vec(vec));
        return compress(ptr);
    }

    trlwe_ptr encryptBatch(const Bytes& data) const {
        Bytes boolArray = toBitVec(data);
        return encryptData(boolArray);
    }

    trlwe_ptr trivialEncryptBatch(const Bytes& data) const {
        Bytes boolArray = toBitVec(data);
        return trivialEncryptData(boolArray);
    }

    trlwe_v2_ptr encryptDataV2(const Bytes& boolArray) const {
        if (boolArray.empty() || boolArray.size() > lvl1param::n) {
            throw std::invalid_argument("level-1 plaintext must contain 1..n bits");
        }
        autohog::SecureRng eng;
        BitsArrayV2 bits = to_bit_array<lvl1param>(boolArray);
        // [FIX] Use ::PKE
        trlwe_v2_t ct = ::PKE::encrypt_bits_trlwe_with_seeded_rlwe_pair<lvl1param>(spk_lvl1, bits, eng);
        return trlwe_v2_ptr(new trlwe_v2(ct, boolArray.size()));
    }

    trlwe_v2_ptr encryptBatchV2(const Bytes& data) const {
        Bytes boolArray = toBitVec(data);
        return encryptDataV2(boolArray);
    }

    vector<uint8_t> encryptDataWithZk(const Bytes& boolArray, const Bytes& metadata) const {
        require_zk_crs("generate zk proof");
        if (boolArray.empty() || boolArray.size() > zk_params.k) {
            throw std::invalid_argument("zk plaintext exceeds the proved bit range");
        }
        autohog::SecureRng engine;
        auto result = ZKPKE::encrypt_bits_with_zk(boolArray, spk_lvl2, engine, zk_params, metadata);
        trlwe_ptr ct(new trlwe(result.ciphertext, boolArray.size()));
        return serialize_ciphertext_with_proof(ct, result.proof);
    }

    vector<uint8_t> encryptWithZk(uint8_t bit, const Bytes& metadata) const {
        Bytes boolArray = { bit };
        return encryptDataWithZk(boolArray, metadata);
    }

    vector<uint8_t> encryptWithZkBatch(const Bytes& data, const Bytes& metadata) const {
        Bytes boolArray = toBitVec(data);
        return encryptDataWithZk(boolArray, metadata);
    }

    bool verifyZkProof(trlwe_ptr& ct, const Bytes& data, const Bytes& metadata) {
        require_zk_crs("verify zk proof");
        ZKProof::Proof proof;
        bool is_ok = deserialize_ciphertext_with_proof(ct, proof, data);
        if (!is_ok) {
            return false;
        }
        auto public_commit = ZKProof::make_public_commit<crs_lvl>(spk_lvl2, ct->data, zk_params.k);
        is_ok = ZKProof::verify<crs_lvl>(proof, zk_params, public_commit, metadata);
        if (!is_ok) {
            return false;
        }
        return true;
    }

    trlwe_ptr decompressV1(const trlwe_v2_ptr& d) const {
        require_valid_trlwe_wrapper<lvl1param>(d, "level-1 compressed ciphertext");
        EvalKey* p_ek = (EvalKey*)&ek;
        auto ksk10_sp = p_ek->get<KeySwitchingKey<lvl10param>>();
        auto &ksk10 = *ksk10_sp;
        auto bkfft02_sp = p_ek->get<BootstrappingKeyFFT<lvl02param>>();
        auto &bkfft02 = *bkfft02_sp;

        trlwe_t packed_pbs_iks;
        HomDecompress::Decompress_PBSiks_lvl1_to_lvl2(
            packed_pbs_iks,
            d->data,
            ksk10,
            bkfft02,
            ahk2,
            d->W);
        return trlwe_ptr(new trlwe(packed_pbs_iks, d->W));
    }

    trlwe_ptr decompressV2(const trlwe_v2_ptr& d) const {
        require_valid_trlwe_wrapper<lvl1param>(d, "level-1 compressed ciphertext");
        EvalKey* p_ek = (EvalKey*)&ek;
        auto bkfft12_sp = p_ek->get<BootstrappingKeyFFT<lvl12param>>();
        auto &bkfft12 = *bkfft12_sp;

        trlwe_t packed_pbs_only;
        HomDecompress::Decompress_PBS_lvl1_to_lvl2(
              packed_pbs_only,
              d->data,
              bkfft12,
              ahk2,
              d->W);
        return trlwe_ptr(new trlwe(packed_pbs_only, d->W));
    }

    trlwe_ptr decompressV3(const trlwe_v2_ptr& d) const {
        require_valid_trlwe_wrapper<lvl1param>(d, "level-1 compressed ciphertext");
        EvalKey* p_ek = (EvalKey*)&ek;
        auto ksk12_sp = p_ek->get<KeySwitchingKey<lvl12param>>();
        auto &ksk12 = *ksk12_sp;

        trlwe_t packed_direct_iks;
        HomDecompress::Decompress_DirectIKS_lvl1_to_lvl2(
              packed_direct_iks,
              d->data,
              ksk12,
              ahk2,
              d->W);
        return trlwe_ptr(new trlwe(packed_direct_iks, d->W));
    }

    Bytes decryptBatch(const tlwe_vec_ptr& v) const {
        require_valid_tlwe_vector(v, "TLWE vector to decrypt");
        const tlwe_vec_t& vec = v->data;
        Bytes result;
        for (size_t i = 0; i < v->W; i += 8) {
            uint8_t k = 0;
            
            const size_t nbits = std::min<size_t>(8, v->W - i);
            for (size_t j = 0; j < nbits; j++) {
                uint8_t bit = bootsSymDecrypt<lvl2param>({vec[i + j]}, sk).at(0);
                if (bit) {
                    k |= (1 << (8 - 1 - j));
                }
            }
            result.push_back(k);
        }
        return result;
    }

    uint8_t decrypt(const tlwe_ptr& bit) const {
        if (!bit) throw std::invalid_argument("invalid TLWE ciphertext to decrypt");
        uint8_t b = bootsSymDecrypt<lvl2param>({bit->data}, sk).at(0);
        return b;
    }

    // ==========================================
    // Deterministic Arithmetic & Conversion
    // ==========================================

    tlwe_vec_ptr trivialAdd(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs) const {
        require_valid_tlwe_vector(lhs, "left trivial operand");
        require_valid_tlwe_vector(rhs, "right trivial operand");
        if (lhs->data.size() != rhs->data.size()) return nullptr;
        size_t size = lhs->data.size();
        tlwe_vec_t result_vec(size);
        bool carry = false;
        // Big-Endian traversal (LSB is at index size-1)
        for (size_t pos = size; pos-- > 0;) {
            const size_t i = pos;
            bool a = (lhs->data[i][lvl2param::n] == lvl2param::μ);
            bool b = (rhs->data[i][lvl2param::n] == lvl2param::μ);
            bool sum = a ^ b ^ carry;
            carry = (a && b) || (carry && (a ^ b));
            tlwe_t res = {};
            res[lvl2param::n] = sum ? lvl2param::μ : -lvl2param::μ;
            result_vec[i] = res;
        }
        return tlwe_vec_ptr(new tlwe_vec(result_vec));
    }

    tlwe_vec_ptr trivialSub(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs) const {
        require_valid_tlwe_vector(lhs, "left trivial operand");
        require_valid_tlwe_vector(rhs, "right trivial operand");
        if (lhs->data.size() != rhs->data.size()) return nullptr;
        size_t size = lhs->data.size();
        tlwe_vec_t result_vec(size);
        bool borrow = false;
        // Big-Endian traversal
        for (size_t pos = size; pos-- > 0;) {
            const size_t i = pos;
            bool a = (lhs->data[i][lvl2param::n] == lvl2param::μ);
            bool b = (rhs->data[i][lvl2param::n] == lvl2param::μ);
            bool diff = a ^ b ^ borrow;
            borrow = (!a && b) || (borrow && !(a ^ b));
            tlwe_t res = {};
            res[lvl2param::n] = diff ? lvl2param::μ : -lvl2param::μ;
            result_vec[i] = res;
        }
        return tlwe_vec_ptr(new tlwe_vec(result_vec));
    }

    tlwe_vec_ptr trivialToStandard(const tlwe_vec_ptr& triv_ct) const {
        require_valid_tlwe_vector(triv_ct, "trivial ciphertext vector");
        size_t bit_len = triv_ct->data.size();
        size_t byte_len = (bit_len + 7) / 8;
        vector<uint8_t> zeros(byte_len, 0);

        // 1. Generate noise (Mask + Error) using internal RNG
        auto zero_trlwe = encryptBatch(zeros); 
        auto zero_ct = decompress_ciphertexts(zero_trlwe);

        tlwe_vec_t res_vec(bit_len);
        for(size_t i = 0; i < bit_len; ++i) {
            res_vec[i] = zero_ct->data[i];
            // 2. Offset correction: Add μ to center the noise around 0
            res_vec[i][lvl2param::n] += lvl2param::μ;
            // 3. Add Trivial Payload
            res_vec[i][lvl2param::n] += triv_ct->data[i][lvl2param::n];
        }
        return tlwe_vec_ptr(new tlwe_vec(res_vec));
    }

    tlwe_vec_ptr standardToTrivial(const tlwe_vec_ptr& std_ct) const {
        require_valid_tlwe_vector(std_ct, "standard ciphertext vector");
        // 1. Decrypt to clean plaintext
        vector<uint8_t> pt = decryptBatch(std_ct);
        // 2. Re-encrypt as trivial (clean mask)
        auto triv_trlwe = trivialEncryptBatch(pt);
        // 3. Decompress back to vector
        return decompress_ciphertexts(triv_trlwe);
    }

    // =========================================================================
    
    // =========================================================================
    
    
    tlwe_vec_ptr add(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, bool compute_overflow = false) const {
        require_compatible_tlwe_vectors(lhs, rhs, "add");
        auto res = HomAdd(lhs->data, rhs->data, lhs->W, ek, compute_overflow);
        return tlwe_vec_ptr(new tlwe_vec(res.SUM));
    }

    tlwe_vec_ptr sub(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, bool compute_overflow = false) const {
        require_compatible_tlwe_vectors(lhs, rhs, "sub");
        auto res = HomSub(lhs->data, rhs->data, lhs->W, ek, compute_overflow);
        return tlwe_vec_ptr(new tlwe_vec(res.SUM));
    }

    
    tlwe_vec_ptr add_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, bool compute_overflow = false, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel add");
        auto res = HomAdd_pal(lhs->data, rhs->data, lhs->W, ek, compute_overflow, threads_hint);
        return tlwe_vec_ptr(new tlwe_vec(res.SUM));
    }

    tlwe_vec_ptr sub_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, bool compute_overflow = false, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel sub");
        auto res = HomSub_pal(lhs->data, rhs->data, lhs->W, ek, compute_overflow, threads_hint);
        return tlwe_vec_ptr(new tlwe_vec(res.SUM));
    }

    // =========================================================================
    
    // =========================================================================

    
    tlwe_vec_ptr add_unsigned(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, bool compute_carry = false) const {
        require_compatible_tlwe_vectors(lhs, rhs, "unsigned add");
        auto res = HomAdd_unsigned(lhs->data, rhs->data, lhs->W, ek, compute_carry);
        return tlwe_vec_ptr(new tlwe_vec(res.SUM));
    }

    tlwe_vec_ptr sub_unsigned(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, bool compute_borrow = false) const {
        require_compatible_tlwe_vectors(lhs, rhs, "unsigned sub");
        auto res = HomSub_unsigned(lhs->data, rhs->data, lhs->W, ek, compute_borrow);
        return tlwe_vec_ptr(new tlwe_vec(res.SUM));
    }

    
    tlwe_vec_ptr add_unsigned_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, bool compute_carry = false) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel unsigned add");
        auto res = HomAdd_unsigned_pal(lhs->data, rhs->data, lhs->W, ek, compute_carry);
        return tlwe_vec_ptr(new tlwe_vec(res.SUM));
    }

    tlwe_vec_ptr sub_unsigned_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, bool compute_borrow = false) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel unsigned sub");
        auto res = HomSub_unsigned_pal(lhs->data, rhs->data, lhs->W, ek, compute_borrow);
        return tlwe_vec_ptr(new tlwe_vec(res.SUM));
    }

    // ==========================================
    
    
    
    // ==========================================

    static std::vector<uint8_t> plain_bytes_to_bits(const Bytes& be, size_t W) {
        std::vector<uint8_t> bits(W, 0);
        
        for (size_t k = 0; k < W; ++k) {           
            const size_t byte_from_end = k / 8;
            if (byte_from_end >= be.size()) break;
            const uint8_t byte = be[be.size() - 1 - byte_from_end];
            bits[W - 1 - k] = (byte >> (k % 8)) & 1;
        }
        return bits;
    }

    
    tlwe_vec_ptr mul_low(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs) const {
        require_compatible_tlwe_vectors(lhs, rhs, "low multiplication");
        auto res = HomMulCompressorLow(lhs->data, rhs->data, lhs->W, ek);
        return tlwe_vec_ptr(new tlwe_vec(res));
    }

    tlwe_vec_ptr mul_low_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs,
                             int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel low multiplication");
        auto res = HomMulCompressorLow_pal(lhs->data, rhs->data, lhs->W, ek, threads_hint);
        return tlwe_vec_ptr(new tlwe_vec(res));
    }

    
    tlwe_vec_ptr mul_plain(const tlwe_vec_ptr& lhs, const Bytes& k) const {
        require_valid_tlwe_vector(lhs, "plain multiplication input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        return tlwe_vec_ptr(new tlwe_vec(HomMulPlain(lhs->data, kb, lhs->W, ek)));
    }
    tlwe_vec_ptr mul_plain_pal(const tlwe_vec_ptr& lhs, const Bytes& k,
                               int threads_hint = -1) const {
        require_valid_tlwe_vector(lhs, "parallel plain multiplication input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        return tlwe_vec_ptr(new tlwe_vec(HomMulPlain_pal(lhs->data, kb, lhs->W, ek, threads_hint)));
    }

    
    tlwe_vec_ptr mul_plain_low(const tlwe_vec_ptr& lhs, const Bytes& k) const {
        require_valid_tlwe_vector(lhs, "low plain multiplication input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        return tlwe_vec_ptr(new tlwe_vec(HomMulPlainLow(lhs->data, kb, lhs->W, ek)));
    }
    tlwe_vec_ptr mul_plain_low_pal(const tlwe_vec_ptr& lhs, const Bytes& k,
                                   int threads_hint = -1) const {
        require_valid_tlwe_vector(lhs, "parallel low plain multiplication input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        return tlwe_vec_ptr(new tlwe_vec(HomMulPlainLow_pal(lhs->data, kb, lhs->W, ek, threads_hint)));
    }

    
    
    
    tlwe_vec_ptr div_plain(const tlwe_vec_ptr& lhs, const Bytes& d) const {
        require_valid_tlwe_vector(lhs, "plain division input");
        const auto db = plain_bytes_to_bits(d, lhs->W);
        return tlwe_vec_ptr(new tlwe_vec(HomDivPlainQuotient(lhs->data, db, lhs->W, ek)));
    }
    tlwe_vec_ptr div_plain_pal(const tlwe_vec_ptr& lhs, const Bytes& d,
                               int threads_hint = -1) const {
        require_valid_tlwe_vector(lhs, "parallel plain division input");
        const auto db = plain_bytes_to_bits(d, lhs->W);
        return tlwe_vec_ptr(new tlwe_vec(
            HomDivPlainQuotient_pal(lhs->data, db, lhs->W, ek, threads_hint)));
    }

    
    std::pair<tlwe_vec_ptr, tlwe_vec_ptr> divmod_plain(
        const tlwe_vec_ptr& lhs, const Bytes& d, int threads_hint = -1) const {
        require_valid_tlwe_vector(lhs, "plain divmod input");
        const auto db = plain_bytes_to_bits(d, lhs->W);
        auto r = HomDivPlain_pal(lhs->data, db, lhs->W, ek, threads_hint);
        return {tlwe_vec_ptr(new tlwe_vec(r.quotient)),
                tlwe_vec_ptr(new tlwe_vec(r.remainder))};
    }

    
    tlwe_vec_ptr mod_plain(const tlwe_vec_ptr& lhs, const Bytes& d) const {
        require_valid_tlwe_vector(lhs, "plain modulo input");
        const auto db = plain_bytes_to_bits(d, lhs->W);
        return tlwe_vec_ptr(new tlwe_vec(HomDivPlain(lhs->data, db, lhs->W, ek).remainder));
    }
    tlwe_vec_ptr mod_plain_pal(const tlwe_vec_ptr& lhs, const Bytes& d,
                               int threads_hint = -1) const {
        require_valid_tlwe_vector(lhs, "parallel plain modulo input");
        const auto db = plain_bytes_to_bits(d, lhs->W);
        return tlwe_vec_ptr(new tlwe_vec(HomDivPlain_pal(lhs->data, db, lhs->W, ek, threads_hint).remainder));
    }

    
    std::vector<tlwe_vec_ptr> mul_plain_low_batch(
        const std::vector<tlwe_vec_ptr>& lhs, const Bytes& k,
        int threads_hint = -1) const {
        if (lhs.empty()) throw std::invalid_argument("plain multiplication batch is empty");
        const size_t W = lhs[0] ? lhs[0]->W : 0;
        std::vector<tlwe_vec_t> as;
        as.reserve(lhs.size());
        for (const auto& p : lhs) {
            require_valid_tlwe_vector(p, "plain multiplication batch item");
            if (p->W != W) throw std::invalid_argument("plain multiplication batch widths differ");
            as.push_back(p->data);
        }
        const auto kb = plain_bytes_to_bits(k, W);
        auto res = HomMulPlainLowBatch(as, kb, static_cast<int>(W), ek, threads_hint);
        std::vector<tlwe_vec_ptr> out;
        out.reserve(res.size());
        for (auto& r : res) out.emplace_back(new tlwe_vec(r));
        return out;
    }

    std::vector<tlwe_vec_ptr> div_plain_batch(
        const std::vector<tlwe_vec_ptr>& lhs, const Bytes& d,
        int threads_hint = -1) const {
        if (lhs.empty()) throw std::invalid_argument("plain division batch is empty");
        const size_t W = lhs[0] ? lhs[0]->W : 0;
        std::vector<tlwe_vec_t> as;
        as.reserve(lhs.size());
        for (const auto& p : lhs) {
            require_valid_tlwe_vector(p, "plain division batch item");
            if (p->W != W) throw std::invalid_argument("plain division batch widths differ");
            as.push_back(p->data);
        }
        const auto db = plain_bytes_to_bits(d, W);
        
        auto res = HomDivPlainQuotientBatch(as, db, static_cast<int>(W), ek, threads_hint);
        std::vector<tlwe_vec_ptr> out;
        out.reserve(res.size());
        for (auto& r : res) out.emplace_back(new tlwe_vec(r));
        return out;
    }

    
    tlwe_vec_ptr eq_plain(const tlwe_vec_ptr& lhs, const Bytes& k) const {
        require_valid_tlwe_vector(lhs, "plain equality input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        tlwe_vec_t r = {HomEQPlain(lhs->data, kb, lhs->W, ek)};
        return tlwe_vec_ptr(new tlwe_vec(r));
    }
    tlwe_vec_ptr gt_plain(const tlwe_vec_ptr& lhs, const Bytes& k) const {
        require_valid_tlwe_vector(lhs, "plain greater-than input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        tlwe_vec_t r = {HomGTPlain(lhs->data, kb, lhs->W, ek)};
        return tlwe_vec_ptr(new tlwe_vec(r));
    }
    tlwe_vec_ptr ge_plain(const tlwe_vec_ptr& lhs, const Bytes& k) const {
        require_valid_tlwe_vector(lhs, "plain greater-or-equal input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        tlwe_vec_t r = {HomGEPlain(lhs->data, kb, lhs->W, ek)};
        return tlwe_vec_ptr(new tlwe_vec(r));
    }
    tlwe_vec_ptr lt_plain(const tlwe_vec_ptr& lhs, const Bytes& k) const {
        require_valid_tlwe_vector(lhs, "plain less-than input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        tlwe_vec_t r = {HomLTPlain(lhs->data, kb, lhs->W, ek)};
        return tlwe_vec_ptr(new tlwe_vec(r));
    }
    tlwe_vec_ptr le_plain(const tlwe_vec_ptr& lhs, const Bytes& k) const {
        require_valid_tlwe_vector(lhs, "plain less-or-equal input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        tlwe_vec_t r = {HomLEPlain(lhs->data, kb, lhs->W, ek)};
        return tlwe_vec_ptr(new tlwe_vec(r));
    }

    
    tlwe_vec_ptr ge_plain_pal(const tlwe_vec_ptr& lhs, const Bytes& k, int threads_hint = -1) const {
        require_valid_tlwe_vector(lhs, "parallel plain greater-or-equal input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        tlwe_vec_t r = {HomGEPlain_pal(lhs->data, kb, lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(r));
    }
    tlwe_vec_ptr gt_plain_pal(const tlwe_vec_ptr& lhs, const Bytes& k, int threads_hint = -1) const {
        require_valid_tlwe_vector(lhs, "parallel plain greater-than input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        tlwe_vec_t r = {HomGTPlain_pal(lhs->data, kb, lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(r));
    }
    tlwe_vec_ptr lt_plain_pal(const tlwe_vec_ptr& lhs, const Bytes& k, int threads_hint = -1) const {
        require_valid_tlwe_vector(lhs, "parallel plain less-than input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        tlwe_vec_t r = {HomLTPlain_pal(lhs->data, kb, lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(r));
    }
    tlwe_vec_ptr le_plain_pal(const tlwe_vec_ptr& lhs, const Bytes& k, int threads_hint = -1) const {
        require_valid_tlwe_vector(lhs, "parallel plain less-or-equal input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        tlwe_vec_t r = {HomLEPlain_pal(lhs->data, kb, lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(r));
    }

    
private:
    tlwe_vec_t flip_msb_copy(const tlwe_vec_ptr& lhs) const {
        require_valid_tlwe_vector(lhs, "signed plain comparison input");
        tlwe_vec_t A = lhs->data;
        if (!A.empty()) HomNOT<lvl2param>(A[0], A[0]);
        return A;
    }
    std::vector<uint8_t> flip_msb_bits(const Bytes& k, size_t W) const {
        auto kb = plain_bytes_to_bits(k, W);
        if (!kb.empty()) kb[0] ^= 1;
        return kb;
    }
public:
    tlwe_vec_ptr ge_plain_signed(const tlwe_vec_ptr& lhs, const Bytes& k, int threads_hint = -1) const {
        tlwe_vec_t A = flip_msb_copy(lhs);
        tlwe_vec_t r = {HomGEPlain_pal(A, flip_msb_bits(k, lhs->W), lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(r));
    }
    tlwe_vec_ptr gt_plain_signed(const tlwe_vec_ptr& lhs, const Bytes& k, int threads_hint = -1) const {
        tlwe_vec_t A = flip_msb_copy(lhs);
        tlwe_vec_t r = {HomGTPlain_pal(A, flip_msb_bits(k, lhs->W), lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(r));
    }
    tlwe_vec_ptr lt_plain_signed(const tlwe_vec_ptr& lhs, const Bytes& k, int threads_hint = -1) const {
        tlwe_vec_t A = flip_msb_copy(lhs);
        tlwe_vec_t r = {HomLTPlain_pal(A, flip_msb_bits(k, lhs->W), lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(r));
    }
    tlwe_vec_ptr le_plain_signed(const tlwe_vec_ptr& lhs, const Bytes& k, int threads_hint = -1) const {
        tlwe_vec_t A = flip_msb_copy(lhs);
        tlwe_vec_t r = {HomLEPlain_pal(A, flip_msb_bits(k, lhs->W), lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(r));
    }

    
    tlwe_vec_ptr and_plain(const tlwe_vec_ptr& lhs, const Bytes& k) const {
        require_valid_tlwe_vector(lhs, "plain AND input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        return tlwe_vec_ptr(new tlwe_vec(HomANDPlain(lhs->data, kb, lhs->W)));
    }
    tlwe_vec_ptr or_plain(const tlwe_vec_ptr& lhs, const Bytes& k) const {
        require_valid_tlwe_vector(lhs, "plain OR input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        return tlwe_vec_ptr(new tlwe_vec(HomORPlain(lhs->data, kb, lhs->W)));
    }
    tlwe_vec_ptr xor_plain(const tlwe_vec_ptr& lhs, const Bytes& k) const {
        require_valid_tlwe_vector(lhs, "plain XOR input");
        const auto kb = plain_bytes_to_bits(k, lhs->W);
        return tlwe_vec_ptr(new tlwe_vec(HomXORPlain(lhs->data, kb, lhs->W)));
    }
    
    tlwe_vec_ptr select_plain(const Bytes& a, const Bytes& b,
                              const tlwe_vec_ptr& sel, size_t W) const {
        require_valid_tlwe_vector(sel, "plain select condition");
        if (sel->W != 1 || W == 0 || W > lvl2param::n) {
            throw std::invalid_argument("invalid plain select width or condition");
        }
        const auto ab = plain_bytes_to_bits(a, W);
        const auto bb = plain_bytes_to_bits(b, W);
        return tlwe_vec_ptr(new tlwe_vec(
            HomSelectPlain(ab, bb, sel->data[0], static_cast<int>(W))));
    }

    tlwe_vec_ptr mul(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs) const {
        require_compatible_tlwe_vectors(lhs, rhs, "multiplication");
        auto res = HomMul(lhs->data, rhs->data, lhs->W, ek);
        return tlwe_vec_ptr(new tlwe_vec(res));
    }

    tlwe_vec_ptr mul_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel multiplication");
        auto res = HomMul_pal(lhs->data, rhs->data, lhs->W, ek, threads_hint);
        return tlwe_vec_ptr(new tlwe_vec(res));
    }

    tlwe_vec_ptr div(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs) const {
        require_compatible_tlwe_vectors(lhs, rhs, "division");
        auto res = HomDiv(lhs->data, rhs->data, lhs->W, ek);
        return tlwe_vec_ptr(new tlwe_vec(res.quotient));
    }

    tlwe_vec_ptr div_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel division");
        auto res = HomDiv_pal(lhs->data, rhs->data, lhs->W, ek, threads_hint);
        return tlwe_vec_ptr(new tlwe_vec(res.quotient));
    }

    tlwe_vec_ptr cast(const tlwe_vec_ptr& in, int to_width, bool signed_extend = false) const {
        require_valid_tlwe_vector(in, "cast input");
        if (to_width <= 0 || static_cast<size_t>(to_width) > lvl2param::n) {
            throw std::invalid_argument("cast width is outside the supported range");
        }
        auto res = HomCast(in->data, in->W, to_width, signed_extend);
        return tlwe_vec_ptr(new tlwe_vec(res));
    }

    // =========================================================================
    
    // =========================================================================

    
    tlwe_vec_ptr gt(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs) const {
        require_compatible_tlwe_vectors(lhs, rhs, "greater-than comparison");
        tlwe_vec_t vec = {HomGT(lhs->data, rhs->data, lhs->W, ek)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr ge(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs) const {
        require_compatible_tlwe_vectors(lhs, rhs, "greater-or-equal comparison");
        tlwe_vec_t vec = {HomGE(lhs->data, rhs->data, lhs->W, ek)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr lt(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs) const {
        require_compatible_tlwe_vectors(lhs, rhs, "less-than comparison");
        tlwe_vec_t vec = {HomLT(lhs->data, rhs->data, lhs->W, ek)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr le(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs) const {
        require_compatible_tlwe_vectors(lhs, rhs, "less-or-equal comparison");
        tlwe_vec_t vec = {HomLE(lhs->data, rhs->data, lhs->W, ek)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr eq(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs) const {
        require_compatible_tlwe_vectors(lhs, rhs, "equality comparison");
        tlwe_vec_t vec = {HomEQ(lhs->data, rhs->data, lhs->W, ek)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr ne(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs) const {
        require_compatible_tlwe_vectors(lhs, rhs, "inequality comparison");
        tlwe_vec_t vec = {HomNE(lhs->data, rhs->data, lhs->W, ek)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }

    
    tlwe_vec_ptr gt_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel greater-than comparison");
        tlwe_vec_t vec = {HomGT_pal(lhs->data, rhs->data, lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr ge_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel greater-or-equal comparison");
        tlwe_vec_t vec = {HomGE_pal(lhs->data, rhs->data, lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr lt_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel less-than comparison");
        tlwe_vec_t vec = {HomLT_pal(lhs->data, rhs->data, lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr le_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel less-or-equal comparison");
        tlwe_vec_t vec = {HomLE_pal(lhs->data, rhs->data, lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr eq_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel equality comparison");
        tlwe_vec_t vec = {HomEQ_pal(lhs->data, rhs->data, lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr ne_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel inequality comparison");
        tlwe_vec_t vec = {HomNE_pal(lhs->data, rhs->data, lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }

    // =========================================================================
    
    // =========================================================================

    
    tlwe_vec_ptr gt_unsigned(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs) const {
        require_compatible_tlwe_vectors(lhs, rhs, "unsigned greater-than comparison");
        tlwe_vec_t vec = {HomGT_unsigned(lhs->data, rhs->data, lhs->W, ek)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr ge_unsigned(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs) const {
        require_compatible_tlwe_vectors(lhs, rhs, "unsigned greater-or-equal comparison");
        tlwe_vec_t vec = {HomGE_unsigned(lhs->data, rhs->data, lhs->W, ek)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr lt_unsigned(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs) const {
        require_compatible_tlwe_vectors(lhs, rhs, "unsigned less-than comparison");
        tlwe_vec_t vec = {HomLT_unsigned(lhs->data, rhs->data, lhs->W, ek)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr le_unsigned(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs) const {
        require_compatible_tlwe_vectors(lhs, rhs, "unsigned less-or-equal comparison");
        tlwe_vec_t vec = {HomLE_unsigned(lhs->data, rhs->data, lhs->W, ek)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }

    
    tlwe_vec_ptr gt_unsigned_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel unsigned greater-than comparison");
        tlwe_vec_t vec = {HomGT_unsigned_pal(lhs->data, rhs->data, lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr ge_unsigned_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel unsigned greater-or-equal comparison");
        tlwe_vec_t vec = {HomGE_unsigned_pal(lhs->data, rhs->data, lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr lt_unsigned_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel unsigned less-than comparison");
        tlwe_vec_t vec = {HomLT_unsigned_pal(lhs->data, rhs->data, lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }
    tlwe_vec_ptr le_unsigned_pal(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "parallel unsigned less-or-equal comparison");
        tlwe_vec_t vec = {HomLE_unsigned_pal(lhs->data, rhs->data, lhs->W, ek, threads_hint)};
        return tlwe_vec_ptr(new tlwe_vec(vec));
    }

    // =========================================================================
    
    // =========================================================================

    tlwe_vec_ptr bit_and(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "bitwise AND");
        auto res = HomAND(lhs->data, rhs->data, lhs->W, ek, threads_hint); 
        return tlwe_vec_ptr(new tlwe_vec(res));
    }

    tlwe_vec_ptr bit_or(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "bitwise OR");
        auto res = HomOR(lhs->data, rhs->data, lhs->W, ek, threads_hint); 
        return tlwe_vec_ptr(new tlwe_vec(res));
    }

    tlwe_vec_ptr bit_xor(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "bitwise XOR");
        auto res = HomXOR(lhs->data, rhs->data, lhs->W, ek, threads_hint); 
        return tlwe_vec_ptr(new tlwe_vec(res));
    }

    tlwe_vec_ptr bit_not(const tlwe_vec_ptr& in) const {
        require_valid_tlwe_vector(in, "bitwise NOT input");
        auto res = HomNOT(in->data, in->W, ek); 
        return tlwe_vec_ptr(new tlwe_vec(res));
    }

    tlwe_vec_ptr select(const tlwe_vec_ptr& lhs, const tlwe_vec_ptr& rhs, const tlwe_vec_ptr& sel, int threads_hint = -1) const {
        require_compatible_tlwe_vectors(lhs, rhs, "select");
        require_valid_tlwe_vector(sel, "select condition");
        if (sel->W != 1) throw std::invalid_argument("select condition must contain exactly one bit");
        auto res = HomSelect(lhs->data, rhs->data, sel->data[0], lhs->W, ek, threads_hint);
        return tlwe_vec_ptr(new tlwe_vec(res));
    }
};

inline std::unique_ptr<FheKey> create_fhe_key() {
    return std::unique_ptr<FheKey>(new FheKey());
}
