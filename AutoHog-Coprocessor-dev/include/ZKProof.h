#pragma once

#include <cstddef>
#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "HomPKE.h"

namespace ZKProof {

constexpr size_t SEED_LEN = 32;

enum class ComputeLoad {
    Proof,
    Verify
};

class Proof {
public:
    void* c_hat_ptr;
    void* c_y_ptr;
    void* pi_ptr;

    Proof();
    ~Proof();
    Proof(const Proof&) = delete;
    Proof& operator=(const Proof&) = delete;
    Proof(Proof&& other) noexcept;
    Proof& operator=(Proof&& other) noexcept;

    std::vector<uint8_t> serialize() const;
    static Proof deserialize(const std::vector<uint8_t>& data);
};

template<typename P>
class PublicParams {
public:
    void* g_list_ptr;
    void* g_hat_list_ptr;

    uint64_t n;
    uint64_t d;
    uint64_t k;
    uint64_t b;
    uint64_t b_r;
    uint64_t q;
    uint64_t t;
    uint64_t msbs_zero_padding_bit_count;
    bool integrity_checked;

    std::array<uint8_t, 16> sid;
    std::array<uint8_t, 16> hash_ds;
    std::array<uint8_t, 16> hash_t_ds;
    std::array<uint8_t, 16> hash_agg_ds;
    std::array<uint8_t, 16> hash_lmap_ds;
    std::array<uint8_t, 16> hash_z_ds;
    std::array<uint8_t, 16> hash_w_ds;
    std::array<uint8_t, 16> hash_gamma_ds;

    PublicParams();
    ~PublicParams();
    PublicParams(const PublicParams&) = delete;
    PublicParams& operator=(const PublicParams&) = delete;
    PublicParams(PublicParams&& other) noexcept;
    PublicParams& operator=(PublicParams&& other) noexcept;
};

template<typename P>
struct PublicCommit {
    uint64_t d;
    uint64_t k;
    std::vector<int64_t> a;
    std::vector<int64_t> b;
    std::vector<int64_t> c1;
    std::vector<int64_t> c2;
};

template<typename P>
struct PrivateCommit {
    uint64_t d;
    uint64_t k;
    std::vector<int64_t> r;
    std::vector<int64_t> e1;
    std::vector<int64_t> e2;
    std::vector<int64_t> m;
};

template<typename P>
PublicParams<P> generate_public_params(
    uint64_t k,
    uint64_t b,
    uint64_t b_r,
    uint64_t q,
    uint64_t t,
    uint64_t msbs_zero_padding_bit_count,
    const std::array<uint8_t, SEED_LEN>& alpha_seed
);

template<typename P>
void save_public_params(const PublicParams<P>& params, std::ostream& ofs);
template<typename P>
void save_public_params(const PublicParams<P>& params, const std::string& filename);
template<typename P>
void save_public_params_to_memory(const PublicParams<P>& params, std::string& data);

template<typename P>
PublicParams<P> load_public_params(const std::string& filename);

template<typename P>
PublicParams<P> load_public_params_from_memory(const char* crs, int len);

template<typename P>
PublicParams<P> load_public_params_from_memory_by_expected_hash(
    const char* crs,
    int len,
    const std::string& expected_hash_hex
);
template<typename P>
PublicParams<P> load_public_params_by_hash(
    const std::string& filename,
    const std::string& hash_filename
);

template<typename P>
PublicParams<P> load_public_params_by_expected_hash(
    const std::string& filename,
    const std::string& expected_hash_hex
);

// Performs full CRS validation: parameter dimensions, group/subgroup checks,
// the required G1 hole, and batched pairing checks for the structured power
// chains. A successful result does not by itself prove that the trapdoor was
// destroyed; production callers must also authenticate the CRS provenance.
template<typename P>
bool validate_public_params(const PublicParams<P>& params);

std::string compute_file_sha256_hex(const std::string& filename);
void write_file_sha256_hex(const std::string& filename, const std::string& hash_filename);

template<typename P>
PublicCommit<P> make_public_commit(
    const PKE::SeededRlwePkPair<P>& spk,
    const TFHEpp::TRLWE<P>& ct,
    uint64_t num_messages
);

template<typename P>
PrivateCommit<P> make_private_commit(
    const std::vector<int64_t>& r,
    const std::vector<int64_t>& e1,
    const std::vector<int64_t>& e2,
    const std::vector<int64_t>& m
);

template<typename P>
Proof prove(
    const PublicParams<P>& params,
    const PublicCommit<P>& public_commit,
    const PrivateCommit<P>& private_commit,
    const std::vector<uint8_t>& metadata,
    ComputeLoad load,
    const std::array<uint8_t, SEED_LEN>& seed
);

template<typename P>
bool verify(
    const Proof& proof,
    const PublicParams<P>& params,
    const PublicCommit<P>& public_commit,
    const std::vector<uint8_t>& metadata
);

}
