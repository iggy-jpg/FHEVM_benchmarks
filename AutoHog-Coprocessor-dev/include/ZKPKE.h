#pragma once
#include "ZKProof.h"
#include "HomPKE.h"
#include "HomPack.h"
#include <limits>
#include <random>
#include <type_traits>
#include <vector>
#include <stdexcept>
#include <array>

namespace ZKPKE {

template<class P>
inline bool has_loaded_crs(const ZKProof::PublicParams<P>& params) {
    return params.g_list_ptr != nullptr &&
           params.g_hat_list_ptr != nullptr &&
           params.n != 0 &&
           params.d != 0 &&
           params.k != 0 &&
           params.integrity_checked;
}

template<class P>
struct EncryptionResult {
    TFHEpp::TRLWE<P> ciphertext;
    ZKProof::Proof proof;
    ZKProof::PrivateCommit<P> witness;
    PKE::SeededRlwePkPair<P> public_key;
};

template<class P, class URBG>
EncryptionResult<P> encrypt_with_zk_proof(
    const std::array<bool, P::n>& bits,
    uint64_t num_bits,
    const PKE::SeededRlwePkPair<P>& public_key,
    URBG& rng,
    const ZKProof::PublicParams<P>& params,
    const std::vector<uint8_t>& metadata
) {
    if (!has_loaded_crs(params)) {
        throw std::runtime_error("missing CRS, cannot generate zk proof");
    }
    const size_t statement_k = params.k;
    if (statement_k == 0 || statement_k > P::n) {
        throw std::runtime_error("invalid CRS statement dimension");
    }
    if (num_bits == 0 || num_bits > statement_k) {
        throw std::runtime_error("num_bits must be in [1, params.k]");
    }

    auto exact_negacyclic_mul_rev = [](const TFHEpp::Polynomial<P>& a_poly,
                                       const std::vector<int64_t>& r_small) {
        std::vector<__int128_t> acc(P::n, 0);
        for (size_t i = 0; i < P::n; ++i) {
            const __int128_t a = static_cast<__int128_t>(static_cast<int64_t>(a_poly[i]));
            for (size_t j = 0; j < P::n; ++j) {
                const size_t r_idx = P::n - 1 - j;
                const size_t idx = i + j;
                const __int128_t term = a * static_cast<__int128_t>(r_small[r_idx]);
                if (idx < P::n) acc[idx] += term;
                else acc[idx - P::n] -= term;
            }
        }

        return acc;
    };

    auto coeffs_to_poly = [](const std::vector<__int128_t>& coeffs) {
        TFHEpp::Polynomial<P> out{};
        for (size_t i = 0; i < P::n; ++i) {
            out[i] = static_cast<typename P::T>(coeffs[i]);
        }
        return out;
    };

    auto exact_divide_by_q = [](__int128_t value,
                                __uint128_t decoded_q,
                                const char* relation_name,
                                size_t index) {
        const __int128_t q = static_cast<__int128_t>(decoded_q);
        if (value % q != 0) {
            throw std::runtime_error(std::string(relation_name) +
                                     " residual is non-zero at index " +
                                     std::to_string(index));
        }
        return static_cast<int64_t>(value / q);
    };

    EncryptionResult<P> result;
    result.public_key = public_key;

    const auto e0 = PKE::sample_gaussian_poly<P>(rng, P::α);
    const auto e1 = PKE::sample_gaussian_poly<P>(rng, P::α);

    std::bernoulli_distribution bit_dist(0.5);
    std::vector<int64_t> r_vec(P::n);
    for (size_t i = 0; i < P::n; ++i) {
        r_vec[i] = bit_dist(rng) ? 1 : 0;
    }

    TFHEpp::Polynomial<P> m{};
    std::vector<int64_t> m_vec(P::n, 0);
    std::vector<int64_t> proof_m_vec(statement_k, 0);
    for (size_t j = 0; j < P::n; ++j) {
        const bool bit = bits[j];
        m_vec[j] = bit ? 1 : 0;
        m[j] = bit ? static_cast<typename P::T>(P::μ)
                   : (typename P::T(0) - static_cast<typename P::T>(P::μ));
    }
    for (size_t i = 0; i < statement_k; ++i) {
        proof_m_vec[i] = m_vec[statement_k - 1 - i];
    }

    TFHEpp::Polynomial<P> tmp{}, ca{}, cb{};
    TFHEpp::Polynomial<P> ai{};
    std::vector<int64_t> a_vec(P::n);
    std::vector<int64_t> b_vec(P::n);
    std::vector<int64_t> c1_vec(P::n);
    std::vector<int64_t> c2_vec(statement_k);

    PKE::uniform_poly_from_seed<P>(ai, result.public_key.seed, 0);
    for (size_t i = 0; i < P::n; ++i) {
        a_vec[i] = static_cast<int64_t>(ai[i]);
        b_vec[i] = static_cast<int64_t>(result.public_key.b[i]);
    }
    const auto ai_times_r_exact = exact_negacyclic_mul_rev(ai, r_vec);
    tmp = coeffs_to_poly(ai_times_r_exact);
    for (size_t j = 0; j < P::n; ++j) {
        ca[j] = tmp[j] + e1[j];
        c1_vec[j] = static_cast<int64_t>(ca[j]);
    }
    PKE::assign_poly_to_trlwe_comp<P>(result.ciphertext, 0, ca);

    const auto b_times_r_exact = exact_negacyclic_mul_rev(result.public_key.b, r_vec);
    cb = coeffs_to_poly(b_times_r_exact);
    for (size_t j = 0; j < P::n; ++j) {
        cb[j] += e0[j] + m[j];
    }
    PKE::assign_poly_to_trlwe_comp<P>(result.ciphertext, P::k, cb);

    std::vector<int64_t> e0_vec(statement_k);
    std::vector<int64_t> e1_vec(P::n);
    for (size_t i = 0; i < P::n; ++i) {
        e1_vec[i] = static_cast<int64_t>(e1[i]);
    }
    for (size_t i = 0; i < statement_k; ++i) {
        e0_vec[i] = static_cast<int64_t>(e0[statement_k - 1 - i]);
    }

    const __uint128_t decoded_q = (params.q == 0)
        ? (__uint128_t(1) << (8u * sizeof(typename P::T)))
        : static_cast<__uint128_t>(params.q);
    const __int128_t delta_enc = static_cast<__int128_t>(decoded_q / params.t);
    std::vector<int64_t> r1_vec(P::n);
    std::vector<int64_t> r2_vec(statement_k);
    for (size_t i = 0; i < P::n; ++i) {
        r1_vec[i] = exact_divide_by_q(
            ai_times_r_exact[i] +
            static_cast<__int128_t>(e1_vec[i]) -
            static_cast<__int128_t>(c1_vec[i]),
            decoded_q,
            "RLWE c1 relation",
            i
        );
    }
    for (size_t i = 0; i < statement_k; ++i) {
        c2_vec[i] = static_cast<int64_t>(cb[statement_k - 1 - i] + static_cast<typename P::T>(P::μ));
        r2_vec[i] = exact_divide_by_q(
            b_times_r_exact[statement_k - 1 - i] +
            static_cast<__int128_t>(e0_vec[i]) +
            delta_enc * static_cast<__int128_t>(proof_m_vec[i]) -
            static_cast<__int128_t>(c2_vec[i]),
            decoded_q,
            "RLWE c2 relation",
            i
        );
    }

    ZKProof::PublicCommit<P> public_commit;
    public_commit.d = P::n;
    public_commit.k = statement_k;
    public_commit.a = std::move(a_vec);
    public_commit.b = std::move(b_vec);
    public_commit.c1 = std::move(c1_vec);
    public_commit.c2 = std::move(c2_vec);

    result.witness = ZKProof::make_private_commit<P>(
        std::move(r_vec),
        std::move(e1_vec),
        std::move(e0_vec),
        std::move(proof_m_vec)
    );

    std::array<uint8_t, ZKProof::SEED_LEN> proof_seed{};
    autohog::secure_random_bytes(proof_seed.data(), proof_seed.size());

    result.proof = ZKProof::prove<P>(
        params,
        public_commit,
        result.witness,
        metadata,
        ZKProof::ComputeLoad::Verify,
        proof_seed
    );

    return result;
}

template<class P, class URBG>
EncryptionResult<P> encrypt_bits_with_zk(
    const std::vector<uint8_t>& bits,
    const PKE::SeededRlwePkPair<P>& public_key,
    URBG& rng,
    const ZKProof::PublicParams<P>& params,
    const std::vector<uint8_t>& metadata
) {
    if (bits.size() > params.k) {
        throw std::runtime_error("Too many bits for the configured ZK statement");
    }

    std::array<bool, P::n> bits_array{};
    bits_array.fill(false);
    for (size_t i = 0; i < bits.size(); ++i) {
        bits_array[i] = (bits[i] != 0);
    }

    return encrypt_with_zk_proof<P>(
        bits_array,
        bits.size(),
        public_key,
        rng,
        params,
        metadata
    );
}

template<class P, class URBG>
EncryptionResult<P> encrypt_bits_with_zk(
    const std::vector<uint8_t>& bits,
    const PKE::SeededRlwePkPair<P>& public_key,
    URBG& rng,
    const ZKProof::PublicParams<P>& params
) {
    std::vector<uint8_t> metadata = {'Z', 'K', 'P', 'K', 'E'};
    return encrypt_bits_with_zk(bits, public_key, rng, params, metadata);
}

template<class P>
bool verify_and_decompress(
    const ZKProof::Proof& proof,
    const ZKProof::PublicParams<P>& params,
    const PKE::SeededRlwePkPair<P>& public_key,
    const TFHEpp::TRLWE<P>& ciphertext,
    std::vector<TFHEpp::TLWE<P>>& output_tlwe,
    uint64_t num_bits,
    const std::vector<uint8_t>& metadata
) {
    // Never return slots that were not covered by the public statement.  The
    // proof commits to exactly params.k c2 coefficients.
    if (num_bits == 0 || num_bits > P::n || num_bits > params.k) {
        return false;
    }
    if (!has_loaded_crs(params)) {
        return false;
    }

    auto public_commit = ZKProof::make_public_commit<P>(public_key, ciphertext, params.k);
    const bool valid = ZKProof::verify<P>(proof, params, public_commit, metadata);
    if (!valid) return false;

    output_tlwe = TRLWE_Decompress<P>(ciphertext, static_cast<int>(num_bits));
    return true;
}

}
