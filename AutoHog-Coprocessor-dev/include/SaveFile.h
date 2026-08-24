#pragma once
#include <memory>

#include <tfhe++.hpp>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>

#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "HomPKE.h"
#include "ZKProof.h"


namespace TFHEpp {

inline bool checkedSerializedSize(std::size_t count, std::size_t element_size,
                                  std::size_t prefix, std::size_t& total) {
    if (count > (std::numeric_limits<std::size_t>::max() - prefix) / element_size) {
        return false;
    }
    total = prefix + count * element_size;
    return true;
}

// ---------- File Size Utilities ----------

extern TLWE<lvl2param> CT_ZERO;
extern TLWE<lvl2param> CT_ONE;
void init_constants(const PKE::SeededRlwePkPair<lvl2param>& spk);

inline void print_size(std::uintmax_t bytes, std::ostream& os = std::cout) {
    double kib = bytes / 1024.0;
    double mib = bytes / (1024.0 * 1024.0);
    double gib = bytes / (1024.0 * 1024.0 * 1024.0);
    os << "Size: " << bytes << " bytes  (" 
       << std::fixed << std::setprecision(2)
       << kib << " KiB, " 
       << mib << " MiB, " 
       << gib << " GiB)\n";
}

/// Get file size
inline std::uintmax_t fileSize(const std::string& path) {
    try {
        return std::filesystem::file_size(path);
    } catch (...) {
        std::ifstream ifs(path, std::ios::binary | std::ios::ate);
        if (ifs) return static_cast<std::uintmax_t>(ifs.tellg());
    }
    return 0;
}

inline void printFileSize(const std::string& path, std::ostream& os = std::cout) {
    auto sz = fileSize(path);
    os << "[printFileSize] File: " << path << "\n";
    print_size(sz, os);
}

// ---------- EvalKey Save/Load ----------
inline bool saveEvalKey(const TFHEpp::EvalKey& ek, std::ostream& ofs) {
    try {
        {
            cereal::BinaryOutputArchive oarchive(ofs);
            oarchive(ek);
        }
        return static_cast<bool>(ofs);
    } catch (const std::exception& e) {
        std::cerr << "[saveEvalKey] Exception: " << e.what() << "\n";
        return false;
    }
}
inline bool saveEvalKey(const TFHEpp::EvalKey& ek, const std::string& path) {
    try {
        std::ofstream ofs(path, std::ios::binary);
        if (!ofs) {
            std::cerr << "[saveEvalKey] Cannot open file: " << path << "\n";
            return false;
        }
        if (!saveEvalKey(ek, ofs)) return false;
        ofs.flush();
        if (!ofs) return false;
        ofs.close();
        std::cout << "[saveEvalKey] Saved to " << path << "\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[saveEvalKey] Exception: " << e.what() << "\n";
        return false;
    }
}

inline bool saveEvalKeyToMemory(const TFHEpp::EvalKey& ek, std::string& data) {
    try {
        std::ostringstream ofs(std::ios::binary);
        if (!saveEvalKey(ek, ofs) || !ofs) return false;
        data = ofs.str();
        std::cout << "[saveEvalKey] Saved to memory "<< "\n";
        return static_cast<bool>(ofs);
    } catch (const std::exception& e) {
        std::cerr << "[saveEvalKey] Exception: " << e.what() << "\n";
        return false;
    }
}
inline bool loadEvalKey(TFHEpp::EvalKey& ek, const std::string& path) {
    try {
        std::ifstream ifs(path, std::ios::binary);
        if (!ifs) {
            std::cerr << "[loadEvalKey] Cannot open file: " << path << "\n";
            return false;
        }
        
        
        TFHEpp::EvalKey fresh_ek; 
        {
            
            cereal::BinaryInputArchive iarchive(ifs);
            iarchive(fresh_ek);
        }
        
        
        ek = std::move(fresh_ek);
        std::cout << "[loadEvalKey] Loaded from " << path << "\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[loadEvalKey] Exception: " << e.what() << "\n";
        return false;
    }
}

// ---------- AnnihilateKey Save/Load (Added) ----------
template<class P>
inline bool saveAnnihilateKey(const TFHEpp::AnnihilateKey<P>& ahk, std::ostream& ofs) {
    try {
        {
            
            cereal::BinaryOutputArchive oarchive(ofs);
            oarchive(ahk);
        }
        
        // Calculate and print size
        return static_cast<bool>(ofs);
    } catch (const std::exception& e) {
        std::cerr << "[saveAnnihilateKey] Exception: " << e.what() << "\n";
        return false;
    }
}

template<class P>
inline bool saveAnnihilateKey(const TFHEpp::AnnihilateKey<P>& ahk, const std::string& path) {
    try {
        std::ofstream ofs(path, std::ios::binary);
        if (!ofs) {
            std::cerr << "[saveAnnihilateKey] Cannot open file: " << path << "\n";
            return false;
        }
        if (!saveAnnihilateKey<P>(ahk, ofs)) return false;
        ofs.flush();
        if (!ofs) return false;
        ofs.close();
        
        // Calculate and print size
        auto sz = fileSize(path);
        std::cout << "[saveAnnihilateKey] Saved to " << path << " ";
        print_size(sz);
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[saveAnnihilateKey] Exception: " << e.what() << "\n";
        return false;
    }
}

template<class P>
inline bool saveAnnihilateKeyToMemory(const TFHEpp::AnnihilateKey<P>& ahk, std::string& data) {
    try {
        std::ostringstream ofs(std::ios::binary);
        if (!saveAnnihilateKey<P>(ahk, ofs) || !ofs) return false;
        data = ofs.str();
        
        // Calculate and print size
        std::cout << "[saveAnnihilateKey] Saved to memory" << " ";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[saveAnnihilateKey] Exception: " << e.what() << "\n";
        return false;
    }
}

template<class P>
inline bool loadAnnihilateKey(TFHEpp::AnnihilateKey<P>& ahk, const std::string& path) {
    try {
        std::ifstream ifs(path, std::ios::binary);
        if (!ifs) {
            std::cerr << "[loadAnnihilateKey] Cannot open file: " << path << "\n";
            return false;
        }
        
        
        
        auto fresh_ahk = std::make_unique<TFHEpp::AnnihilateKey<P>>();
        {
            cereal::BinaryInputArchive iarchive(ifs);
            iarchive(*fresh_ahk);
        }
        ahk = std::move(*fresh_ahk);
        std::cout << "[loadAnnihilateKey] Loaded from " << path << "\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[loadAnnihilateKey] Exception: " << e.what() << "\n";
        return false;
    }
}

// ---------- SecretKey Save/Load ----------
inline bool saveSecretKey(const TFHEpp::SecretKey& sk, std::ostream& ofs) {
    try {
        {
            
            cereal::BinaryOutputArchive oarchive(ofs);
            oarchive(sk);
        }

        return static_cast<bool>(ofs);
    } catch (const std::exception& e) {
        std::cerr << "[saveSecretKey] Exception: " << e.what() << "\n";
        return false;
    }
}
inline bool saveSecretKey(const TFHEpp::SecretKey& sk, const std::string& path) {
    try {
        std::ofstream ofs(path, std::ios::binary);
        if (!ofs) {
            std::cerr << "[saveSecretKey] Cannot open file: " << path << "\n";
            return false;
        }
        if (!saveSecretKey(sk, ofs)) return false;
        ofs.flush();
        if (!ofs) return false;
        ofs.close();

        std::cout << "[saveSecretKey] Saved to " << path << "\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[saveSecretKey] Exception: " << e.what() << "\n";
        return false;
    }
}

inline bool saveSecretKeyToMemory(const TFHEpp::SecretKey& sk, std::string& data) {
    try {
        std::ostringstream ofs(std::ios::binary);
        if (!saveSecretKey(sk, ofs) || !ofs) return false;
        data = ofs.str();

        std::cout << "[saveSecretKey] Saved to memory" << "\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[saveSecretKey] Exception: " << e.what() << "\n";
        return false;
    }
}

inline bool loadSecretKey(TFHEpp::SecretKey& sk, std::istream& ifs) {
    try {
        TFHEpp::SecretKey fresh_sk; 
        {
            cereal::BinaryInputArchive iarchive(ifs);
            iarchive(fresh_sk);
        }
        sk = std::move(fresh_sk);
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[loadSecretKey] Exception: " << e.what() << "\n";
        return false;
    }
}

inline bool loadSecretKey(TFHEpp::SecretKey& sk, const std::string& path) {
    try {
        std::ifstream ifs(path, std::ios::binary);
        if (!ifs) {
            std::cerr << "[loadSecretKey] Cannot open file: " << path << "\n";
            return false;
        }
        if (!loadSecretKey(sk, ifs)) return false;
        std::cout << "[loadSecretKey] Loaded from " << path << "\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[loadSecretKey] Exception: " << e.what() << "\n";
        return false;
    }
}

inline bool loadSecretKeyFromMemory(TFHEpp::SecretKey& sk, const std::string& data) {
    try {
        std::istringstream ifs(data, std::ios::binary);
        if (!ifs) {
            std::cerr << "[loadSecretKey] Cannot open file: " << "\n";
            return false;
        }
        if (!loadSecretKey(sk, ifs)) return false;
        std::cout << "[loadSecretKey] Loaded " << "\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[loadSecretKey] Exception: " << e.what() << "\n";
        return false;
    }
}
// ---------- TLWE Vector Save/Load ----------
template <typename P>
inline bool saveCipherVec(const std::vector<TFHEpp::TLWE<P>>& vec,
                          const std::string& path,
                          bool /*alsoPrintSize*/ = true) {
    try {
        std::ofstream ofs(path, std::ios::binary);
        if (!ofs) {
            std::cerr << "[saveCipherVec] Cannot open file: " << path << "\n";
            return false;
        }
        const std::size_t n = vec.size();
        if (n == 0 || n > P::n) return false;
        ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));
        if (n) {
            ofs.write(reinterpret_cast<const char*>(vec.data()),
                      static_cast<std::streamsize>(sizeof(TFHEpp::TLWE<P>) * n));
        }
        ofs.flush();
        if (!ofs) return false;
        ofs.close();
        std::cout << "[saveCipherVec] Saved " << n << " TLWE to " << path << "\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[saveCipherVec] Exception: " << e.what() << "\n";
        return false;
    }
}

template <typename P>
inline bool loadCipherVec(std::vector<TFHEpp::TLWE<P>>& vec,
                          const std::string& path) {
    try {
        std::ifstream ifs(path, std::ios::binary);
        if (!ifs) {
            std::cerr << "[loadCipherVec] Cannot open file: " << path << "\n";
            return false;
        }
        std::size_t n = 0;
        ifs.read(reinterpret_cast<char*>(&n), sizeof(n));
        std::size_t expected_size = 0;
        if (!ifs || n > P::n ||
            !checkedSerializedSize(n, sizeof(TFHEpp::TLWE<P>), sizeof(n), expected_size) ||
            fileSize(path) != expected_size) {
            std::cerr << "[loadCipherVec] Invalid element count or file size\n";
            return false;
        }
        std::vector<TFHEpp::TLWE<P>> loaded(n);
        if (n) {
            ifs.read(reinterpret_cast<char*>(loaded.data()),
                     static_cast<std::streamsize>(sizeof(TFHEpp::TLWE<P>) * n));
            if (!ifs) {
                std::cerr << "[loadCipherVec] Read truncated! expected "
                          << (sizeof(TFHEpp::TLWE<P>) * n) << " bytes\n";
                return false;
            }
        }
        vec = std::move(loaded);
        std::cout << "[loadCipherVec] Loaded " << n << " TLWE from " << path << "\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[loadCipherVec] Exception: " << e.what() << "\n";
        return false;
    }
}

// ---------- PublicKey Save/Load ----------
template <typename P>
inline bool savePublicKey(const PKE::LwePublicKey<P>& pk,
                          const std::string& path) {
    try {
        std::ofstream ofs(path, std::ios::binary);
        if (!ofs) {
            std::cerr << "[savePublicKey] Cannot open file: " << path << "\n";
            return false;
        }
        const std::size_t n = pk.zeros.size();
        if (n != PKE::zero_count_for_pk<P>()) return false;
        ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));
        if (n) {
            ofs.write(reinterpret_cast<const char*>(pk.zeros.data()),
                      static_cast<std::streamsize>(sizeof(TFHEpp::TLWE<P>) * n));
        }
        ofs.flush();
        if (!ofs) return false;
        ofs.close();
        std::cout << "[savePublicKey] Saved " << n << " TLWE to " << path << "\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[savePublicKey] Exception: " << e.what() << "\n";
        return false;
    }
}

template <typename P>
inline bool loadPublicKey(PKE::LwePublicKey<P>& pk,
                          const std::string& path) {
    try {
        std::ifstream ifs(path, std::ios::binary);
        if (!ifs) {
            std::cerr << "[loadPublicKey] Cannot open file: " << path << "\n";
            return false;
        }
        std::size_t n = 0;
        ifs.read(reinterpret_cast<char*>(&n), sizeof(n));
        const std::size_t required_n = PKE::zero_count_for_pk<P>();
        std::size_t expected_size = 0;
        if (!ifs || n != required_n ||
            !checkedSerializedSize(n, sizeof(TFHEpp::TLWE<P>), sizeof(n), expected_size) ||
            fileSize(path) != expected_size) {
            std::cerr << "[loadPublicKey] Invalid element count or file size\n";
            return false;
        }
        std::vector<TFHEpp::TLWE<P>> loaded(n);
        if (n) {
            ifs.read(reinterpret_cast<char*>(loaded.data()),
                     static_cast<std::streamsize>(sizeof(TFHEpp::TLWE<P>) * n));
            if (!ifs) {
                std::cerr << "[loadPublicKey] Read truncated! expected "
                          << (sizeof(TFHEpp::TLWE<P>) * n) << " bytes\n";
                return false;
            }
        }
        pk.zeros = std::move(loaded);
        std::cout << "[loadPublicKey] Loaded " << n << " TLWE from " << path << "\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[loadPublicKey] Exception: " << e.what() << "\n";
        return false;
    }
}

// ---------- CompressedPublicKey ----------
template <typename P>
inline bool saveCompressedPublicKey(const PKE::SeededLwePublicKey<P>& cpk,
                                    const std::string& path,
                                    bool /*alsoPrintSize*/ = true) {
    try {
        std::ofstream ofs(path, std::ios::binary);
        if (!ofs) {
            std::cerr << "[saveCompressedPublicKey] Cannot open file: " << path << "\n";
            return false;
        }
        ofs.write(reinterpret_cast<const char*>(&cpk.seed), sizeof(cpk.seed));
        std::size_t n = cpk.b_list.size();
        if (n != PKE::zero_count_for_pk<P>()) return false;
        ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));
        if (n) {
            ofs.write(reinterpret_cast<const char*>(cpk.b_list.data()),
                      static_cast<std::streamsize>(sizeof(typename P::T) * n));
        }
        ofs.flush();
        if (!ofs) return false;
        ofs.close();
        std::cout << "[saveCompressedPublicKey] Saved to " << path
                  << " (items=" << n << ")\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[saveCompressedPublicKey] Exception: " << e.what() << "\n";
        return false;
    }
}

template <typename P>
inline bool loadCompressedPublicKey(PKE::SeededLwePublicKey<P>& cpk,
                                    const std::string& path) {
    try {
        std::ifstream ifs(path, std::ios::binary);
        if (!ifs) {
            std::cerr << "[loadCompressedPublicKey] Cannot open file: " << path << "\n";
            return false;
        }
        decltype(cpk.seed) loaded_seed{};
        ifs.read(reinterpret_cast<char*>(&loaded_seed), sizeof(loaded_seed));
        if (!ifs) {
            std::cerr << "[loadCompressedPublicKey] Failed reading seed\n";
            return false;
        }
        std::size_t n = 0;
        ifs.read(reinterpret_cast<char*>(&n), sizeof(n));
        if (!ifs) {
            std::cerr << "[loadCompressedPublicKey] Failed reading size\n";
            return false;
        }
        const std::size_t required_n = PKE::zero_count_for_pk<P>();
        std::size_t expected_size = 0;
        const std::size_t prefix = sizeof(loaded_seed) + sizeof(n);
        if (n != required_n ||
            !checkedSerializedSize(n, sizeof(typename P::T), prefix, expected_size) ||
            fileSize(path) != expected_size) {
            std::cerr << "[loadCompressedPublicKey] Invalid element count or file size\n";
            return false;
        }
        std::vector<typename P::T> loaded_b(n);
        if (n) {
            ifs.read(reinterpret_cast<char*>(loaded_b.data()),
                     static_cast<std::streamsize>(sizeof(typename P::T) * n));
            if (!ifs) {
                std::cerr << "[loadCompressedPublicKey] Read truncated! expected "
                          << (sizeof(typename P::T) * n) << " bytes\n";
                return false;
            }
        }
        cpk.seed = loaded_seed;
        cpk.b_list = std::move(loaded_b);
        std::cout << "[loadCompressedPublicKey] Loaded from " << path
                  << " (items=" << n << ")\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[loadCompressedPublicKey] Exception: " << e.what() << "\n";
        return false;
    }
}

// ---------- RLWE PublicKey (OpenFHE-style) Save/Load ----------
template <typename P>
inline bool saveRlwePublicKey(const PKE::RlwePkPair<P>& pk,
                              const std::string& path) {
    try {
        std::ofstream ofs(path, std::ios::binary);
        if (!ofs) {
            std::cerr << "[saveRlwePublicKey] Cannot open file: " << path << "\n";
            return false;
        }
        for (size_t i = 0; i < P::k; ++i)
            ofs.write(reinterpret_cast<const char*>(pk.a[i].data()),
                      static_cast<std::streamsize>(sizeof(typename P::T) * P::n));
        ofs.write(reinterpret_cast<const char*>(pk.b.data()),
                  static_cast<std::streamsize>(sizeof(typename P::T) * P::n));
        ofs.flush();
        if (!ofs) return false;
        ofs.close();
        std::cout << "[saveRlwePublicKey] Saved to " << path << " (k=" << P::k << ", n=" << P::n << ")\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[saveRlwePublicKey] Exception: " << e.what() << "\n";
        return false;
    }
}

template <typename P>
inline bool loadRlwePublicKey(PKE::RlwePkPair<P>& pk,
                              const std::string& path) {
    try {
        constexpr std::uintmax_t expected_size =
            static_cast<std::uintmax_t>(P::k + 1) * sizeof(typename P::T) * P::n;
        if (fileSize(path) != expected_size) return false;
        std::ifstream ifs(path, std::ios::binary);
        if (!ifs) {
            std::cerr << "[loadRlwePublicKey] Cannot open file: " << path << "\n";
            return false;
        }
        PKE::RlwePkPair<P> loaded{};
        for (size_t i = 0; i < P::k; ++i) {
            ifs.read(reinterpret_cast<char*>(loaded.a[i].data()),
                     static_cast<std::streamsize>(sizeof(typename P::T) * P::n));
            if (!ifs) {
                std::cerr << "[loadRlwePublicKey] Read truncated on a[" << i << "]\n";
                return false;
            }
        }
        ifs.read(reinterpret_cast<char*>(loaded.b.data()),
                 static_cast<std::streamsize>(sizeof(typename P::T) * P::n));
        if (!ifs) {
            std::cerr << "[loadRlwePublicKey] Read truncated on b\n";
            return false;
        }
        pk = std::move(loaded);
        std::cout << "[loadRlwePublicKey] Loaded from " << path << " (k=" << P::k << ", n=" << P::n << ")\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[loadRlwePublicKey] Exception: " << e.what() << "\n";
        return false;
    }
}

template <typename P>
inline bool saveTLWE(const TFHEpp::TLWE<P>& ct, const std::string& path) {
    try {
        std::ofstream ofs(path, std::ios::binary);
        if (!ofs) return false;
        
        ofs.write(reinterpret_cast<const char*>(ct.data()),
                  static_cast<std::streamsize>(sizeof(typename P::T) * (P::n + 1)));
        ofs.flush();
        return static_cast<bool>(ofs);
    } catch (...) { return false; }
}

template <typename P>
inline bool loadTLWE(TFHEpp::TLWE<P>& ct, const std::string& path) {
    try {
        constexpr std::uintmax_t expected_size =
            sizeof(typename P::T) * (P::n + 1);
        if (fileSize(path) != expected_size) return false;
        std::ifstream ifs(path, std::ios::binary);
        if (!ifs) return false;
        TFHEpp::TLWE<P> loaded{};
        ifs.read(reinterpret_cast<char*>(loaded.data()),
                 static_cast<std::streamsize>(sizeof(typename P::T) * (P::n + 1)));
        if (!ifs) return false;
        ct = std::move(loaded);
        return true;
    } catch (...) { return false; }
}


// ---------- TRLWE Save/Load (Single) ----------
template <typename P>
inline bool saveTRLWE(const TFHEpp::TRLWE<P>& ct, const std::string& path) {
    try {
        std::ofstream ofs(path, std::ios::binary);
        if (!ofs) {
            std::cerr << "[saveTRLWE] Cannot open file: " << path << "\n";
            return false;
        }
        for (size_t i = 0; i < P::k + 1; ++i)
            ofs.write(reinterpret_cast<const char*>(ct[i].data()),
                      static_cast<std::streamsize>(sizeof(typename P::T) * P::n));
        ofs.flush();
        if (!ofs) return false;
        ofs.close();
        std::cout << "[saveTRLWE] Saved to " << path << " (k=" << P::k << ", n=" << P::n << ")\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[saveTRLWE] Exception: " << e.what() << "\n";
        return false;
    }
}

template <typename P>
inline bool loadTRLWE(TFHEpp::TRLWE<P>& ct, const std::string& path) {
    try {
        constexpr std::uintmax_t expected_size =
            static_cast<std::uintmax_t>(P::k + 1) * sizeof(typename P::T) * P::n;
        if (fileSize(path) != expected_size) return false;
        std::ifstream ifs(path, std::ios::binary);
        if (!ifs) {
            std::cerr << "[loadTRLWE] Cannot open file: " << path << "\n";
            return false;
        }
        TFHEpp::TRLWE<P> loaded{};
        for (size_t i = 0; i < P::k + 1; ++i) {
            ifs.read(reinterpret_cast<char*>(loaded[i].data()),
                     static_cast<std::streamsize>(sizeof(typename P::T) * P::n));
            if (!ifs) {
                std::cerr << "[loadTRLWE] Read truncated on component " << i << "\n";
                return false;
            }
        }
        ct = std::move(loaded);
        std::cout << "[loadTRLWE] Loaded from " << path << " (k=" << P::k << ", n=" << P::n << ")\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[loadTRLWE] Exception: " << e.what() << "\n";
        return false;
    }
}

// ==========================================
// [FIX] Extended SeededRlwePublicKey Save/Load
// ==========================================


/**
 * Saves a Seeded RLWE Public Key to a binary file.
 * * [Modification]: If the parameter set P is `lvl2param` (typically used for the Server Key),
 * this function automatically appends the global constants `CT_ZERO` and `CT_ONE` to the file.
 * This ensures that all servers loading this key will use bit-exact identical constants,
 * guaranteeing deterministic homomorphic operations.
 */
template <typename P>
inline bool saveSeededRlwePublicKey(const PKE::SeededRlwePkPair<P>& spk,
                                    std::ostream& ofs) {
    try {
        // 1. Save the Seed (used to regenerate the 'a' part of the key)
        ofs.write(reinterpret_cast<const char*>(&spk.seed), sizeof(spk.seed));
        
        // 2. Save the 'b' part of the key (AS + E)
        ofs.write(reinterpret_cast<const char*>(spk.b.data()),
                  static_cast<std::streamsize>(sizeof(typename P::T) * P::n));

        // 3. [NEW] Conditionally save constants for Level 2 keys
        if constexpr (std::is_same_v<P, lvl2param>) {
            // Write CT_ZERO
            ofs.write(reinterpret_cast<const char*>(TFHEpp::CT_ZERO.data()), sizeof(TFHEpp::CT_ZERO));
            // Write CT_ONE
            ofs.write(reinterpret_cast<const char*>(TFHEpp::CT_ONE.data()), sizeof(TFHEpp::CT_ONE));
            
            std::cout << "[saveSeededRlwePublicKey] info: Attached CT_ZERO and CT_ONE to the key file.\n";
        }

        return static_cast<bool>(ofs);
    } catch (const std::exception& e) {
        std::cerr << "[saveSeededRlwePublicKey] Exception: " << e.what() << "\n";
        return false;
    } catch (...) {
        std::cerr << "[saveSeededRlwePublicKey] Unknown exception occurred.\n";
        return false;
    }
}
template <typename P>
inline bool saveSeededRlwePublicKey(const PKE::SeededRlwePkPair<P>& spk,
                                    const std::string& path) {
    try {
        std::ofstream ofs(path, std::ios::binary);
        if (!ofs) {
            std::cerr << "[saveSeededRlwePublicKey] Error: Cannot open file: " << path << "\n";
            return false;
        }

        if (!saveSeededRlwePublicKey(spk, ofs)) return false;
        ofs.flush();
        if (!ofs) return false;
        ofs.close();
        
        // Optional: Print file size for verification
        auto sz = fileSize(path);
        std::cout << "[saveSeededRlwePublicKey] Success: Saved to " << path 
                  << " (Size: " << sz << " bytes)\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[saveSeededRlwePublicKey] Exception: " << e.what() << "\n";
        return false;
    } catch (...) {
        std::cerr << "[saveSeededRlwePublicKey] Unknown exception occurred.\n";
        return false;
    }
}

template <typename P>
inline bool saveSeededRlwePublicKeyToMemory(const PKE::SeededRlwePkPair<P>& spk,
                                            std::string& data) {
    try {
        std::ostringstream ofs(std::ios::binary);

        if (!saveSeededRlwePublicKey(spk, ofs) || !ofs) return false;
        data = ofs.str();
        
        // Optional: Print file size for verification
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[saveSeededRlwePublicKey] Exception: " << e.what() << "\n";
        return false;
    } catch (...) {
        std::cerr << "[saveSeededRlwePublicKey] Unknown exception occurred.\n";
        return false;
    }
}

/**
 * Loads a Seeded RLWE Public Key from a binary file.
 * For `lvl2param`, the legacy serialized constant fields are consumed only for
 * file-format compatibility. They are never trusted as process-global state;
 * deterministic trivial constants are rebuilt locally instead.
 */
template <typename P>
inline bool loadSeededRlwePublicKey(PKE::SeededRlwePkPair<P>& spk,
                                    const std::string& path) {
    try {
        constexpr std::uintmax_t base_size =
            sizeof(decltype(spk.seed)) + sizeof(typename P::T) * P::n;
        constexpr std::uintmax_t expected_size = base_size +
            (std::is_same_v<P, lvl2param>
                 ? 2 * sizeof(TFHEpp::TLWE<lvl2param>)
                 : 0);
        if (fileSize(path) != expected_size) {
            std::cerr << "[loadSeededRlwePublicKey] Error: Invalid file size.\n";
            return false;
        }

        std::ifstream ifs(path, std::ios::binary);
        if (!ifs) {
            std::cerr << "[loadSeededRlwePublicKey] Error: Cannot open file: " << path << "\n";
            return false;
        }

        PKE::SeededRlwePkPair<P> loaded_spk{};
        ifs.read(reinterpret_cast<char*>(&loaded_spk.seed), sizeof(loaded_spk.seed));
        if (!ifs) {
            std::cerr << "[loadSeededRlwePublicKey] Error: Failed to read seed.\n";
            return false;
        }

        ifs.read(reinterpret_cast<char*>(loaded_spk.b.data()),
                 static_cast<std::streamsize>(sizeof(typename P::T) * P::n));
        if (!ifs) {
            std::cerr << "[loadSeededRlwePublicKey] Error: Failed to read key body.\n";
            return false;
        }

        if constexpr (std::is_same_v<P, lvl2param>) {
            // Retain the serialized fields for format compatibility, but never
            // trust caller-controlled constants as process-global state.
            TFHEpp::TLWE<lvl2param> serialized_zero{};
            TFHEpp::TLWE<lvl2param> serialized_one{};
            ifs.read(reinterpret_cast<char*>(serialized_zero.data()), sizeof(serialized_zero));
            ifs.read(reinterpret_cast<char*>(serialized_one.data()), sizeof(serialized_one));
            if (!ifs) {
                std::cerr << "[loadSeededRlwePublicKey] Error: Truncated constants.\n";
                return false;
            }

            TFHEpp::init_constants(loaded_spk);
        }

        spk = std::move(loaded_spk);
        std::cout << "[loadSeededRlwePublicKey] Loaded key from " << path << "\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[loadSeededRlwePublicKey] Exception: " << e.what() << "\n";
        return false;
    } catch (...) {
        std::cerr << "[loadSeededRlwePublicKey] Unknown exception occurred.\n";
        return false;
    }
}


// ---------- SealedCiphertext (Ciphertext + ZK Proof + Metadata) Save/Load ----------

// File Format:
//   [4 bytes]  magic "ZKS1"
//   [TRLWE]    ciphertext (k+1) * n * sizeof(T)
//   [8 bytes]  proof_size
//   [variable] proof serialized data
//   [8 bytes]  metadata_size
//   [variable] metadata

template <typename P>
inline bool saveSealedCiphertext(const TFHEpp::TRLWE<P>& ct,
                                  const ZKProof::Proof& proof,
                                  const std::vector<uint8_t>& metadata,
                                  const std::string& path) {
    try {
        constexpr uint64_t kMaxProofBytes = 1ULL << 20;
        constexpr uint64_t kMaxMetadataBytes = 1ULL << 20;
        if (metadata.size() > kMaxMetadataBytes) {
            std::cerr << "[saveSealedCiphertext] Metadata is too large\n";
            return false;
        }
        auto proof_bytes = proof.serialize();
        if (proof_bytes.empty() || proof_bytes.size() > kMaxProofBytes) {
            std::cerr << "[saveSealedCiphertext] Proof size is invalid\n";
            return false;
        }
        std::ofstream ofs(path, std::ios::binary);
        if (!ofs) {
            std::cerr << "[saveSealedCiphertext] Cannot open file: " << path << "\n";
            return false;
        }

        // Magic
        const char magic[4] = {'Z', 'K', 'S', '1'};
        ofs.write(magic, 4);

        // Ciphertext
        for (size_t i = 0; i < P::k + 1; ++i) {
            ofs.write(reinterpret_cast<const char*>(ct[i].data()),
                      static_cast<std::streamsize>(sizeof(typename P::T) * P::n));
        }

        // Proof
        uint64_t proof_size = proof_bytes.size();
        ofs.write(reinterpret_cast<const char*>(&proof_size), sizeof(proof_size));
        ofs.write(reinterpret_cast<const char*>(proof_bytes.data()),
                  static_cast<std::streamsize>(proof_size));

        // Metadata
        uint64_t meta_size = metadata.size();
        ofs.write(reinterpret_cast<const char*>(&meta_size), sizeof(meta_size));
        if (meta_size > 0) {
            ofs.write(reinterpret_cast<const char*>(metadata.data()),
                      static_cast<std::streamsize>(meta_size));
        }

        ofs.flush();
        if (!ofs) return false;
        ofs.close();

        auto total = fileSize(path);
        std::cout << "[saveSealedCiphertext] Saved to " << path << "\n";
        print_size(total);
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[saveSealedCiphertext] Exception: " << e.what() << "\n";
        return false;
    }
}

template <typename P>
inline bool loadSealedCiphertext(TFHEpp::TRLWE<P>& ct,
                                  ZKProof::Proof& proof,
                                  std::vector<uint8_t>& metadata,
                                  const std::string& path) {
    try {
        constexpr uint64_t kMaxProofBytes = 1ULL << 20;
        constexpr uint64_t kMaxMetadataBytes = 1ULL << 20;
        const uint64_t total_size = fileSize(path);
        const uint64_t ciphertext_bytes =
            static_cast<uint64_t>(P::k + 1) * sizeof(typename P::T) * P::n;
        const uint64_t minimum_size = 4 + ciphertext_bytes + 2 * sizeof(uint64_t);
        if (total_size < minimum_size) {
            std::cerr << "[loadSealedCiphertext] File is too small\n";
            return false;
        }
        std::ifstream ifs(path, std::ios::binary);
        if (!ifs) {
            std::cerr << "[loadSealedCiphertext] Cannot open file: " << path << "\n";
            return false;
        }

        // Magic
        char magic[4];
        ifs.read(magic, 4);
        if (std::memcmp(magic, "ZKS1", 4) != 0) {
            std::cerr << "[loadSealedCiphertext] Invalid file format (bad magic)\n";
            return false;
        }

        TFHEpp::TRLWE<P> loaded_ct{};
        for (size_t i = 0; i < P::k + 1; ++i) {
            ifs.read(reinterpret_cast<char*>(loaded_ct[i].data()),
                     static_cast<std::streamsize>(sizeof(typename P::T) * P::n));
            if (!ifs) {
                std::cerr << "[loadSealedCiphertext] Read truncated on ct[" << i << "]\n";
                return false;
            }
        }

        // Proof
        uint64_t proof_size = 0;
        ifs.read(reinterpret_cast<char*>(&proof_size), sizeof(proof_size));
        if (!ifs || proof_size == 0 || proof_size > kMaxProofBytes) return false;
        const uint64_t after_proof_size = 4 + ciphertext_bytes + sizeof(uint64_t);
        if (proof_size > total_size - after_proof_size) return false;

        std::vector<uint8_t> proof_bytes(proof_size);
        ifs.read(reinterpret_cast<char*>(proof_bytes.data()),
                 static_cast<std::streamsize>(proof_size));
        if (!ifs) {
            std::cerr << "[loadSealedCiphertext] Read truncated on proof\n";
            return false;
        }
        ZKProof::Proof loaded_proof = ZKProof::Proof::deserialize(proof_bytes);

        // Metadata
        uint64_t meta_size = 0;
        ifs.read(reinterpret_cast<char*>(&meta_size), sizeof(meta_size));
        if (!ifs || meta_size > kMaxMetadataBytes) return false;
        const uint64_t expected_total = after_proof_size + proof_size +
                                        sizeof(uint64_t) + meta_size;
        if (expected_total != total_size) {
            std::cerr << "[loadSealedCiphertext] Invalid proof/metadata sizes\n";
            return false;
        }

        std::vector<uint8_t> loaded_metadata(static_cast<size_t>(meta_size));
        if (meta_size > 0) {
            ifs.read(reinterpret_cast<char*>(loaded_metadata.data()),
                     static_cast<std::streamsize>(meta_size));
            if (!ifs) {
                std::cerr << "[loadSealedCiphertext] Read truncated on metadata\n";
                return false;
            }
        }

        ct = std::move(loaded_ct);
        proof = std::move(loaded_proof);
        metadata = std::move(loaded_metadata);

        std::cout << "[loadSealedCiphertext] Loaded from " << path << "\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[loadSealedCiphertext] Exception: " << e.what() << "\n";
        return false;
    }
}


}
