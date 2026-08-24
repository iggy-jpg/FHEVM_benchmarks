#pragma once
#include <cmath>
#include <limits>
#include <vector>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include "tfhe++.hpp"
#include "SecureRandom.h"

namespace PKE {

// Convert a normalized real-valued noise sample to the integer Torus used by
// TFHE. Keep 2^BITS in the floating-point domain: for a 64-bit Torus,
// `1ULL << 64` is an invalid shift and 2^64 cannot be narrowed to uint64_t.
template<class Torus>
inline Torus torus_from_normalized_noise(long double noise)
{
    static_assert(std::is_integral_v<Torus> && !std::is_same_v<Torus, bool>,
                  "Torus must be a non-bool integral type");

    using Uns = std::make_unsigned_t<Torus>;
    constexpr int BITS = std::numeric_limits<Uns>::digits;
    const long double scaled = std::ldexp(noise, BITS);
    const long double min_value = static_cast<long double>(
        (std::numeric_limits<long long>::min)());
    const long double max_value = static_cast<long double>(
        (std::numeric_limits<long long>::max)());

    if (!std::isfinite(scaled) || scaled < min_value || scaled > max_value) {
        throw std::overflow_error("Gaussian Torus noise scaling overflow");
    }

    const long long rounded = std::llround(scaled);
    return static_cast<Torus>(static_cast<Uns>(rounded));
}

template<class P>
struct LwePublicKey {
    std::vector<TFHEpp::TLWE<P>> zeros;
};

template<class P>
inline void tlwe_add_inplace(TFHEpp::TLWE<P>& out, const TFHEpp::TLWE<P>& in) {
    constexpr size_t L = P::n + 1;
    for (size_t i = 0; i < L; ++i) out[i] += in[i];
}

template<class P>
inline void tlwe_add_selected_inplace(TFHEpp::TLWE<P>& out,
                                      const TFHEpp::TLWE<P>& in,
                                      uint32_t bit) {
    constexpr size_t L = P::n + 1;
    const typename P::T mask = typename P::T(0) - typename P::T(bit & 1u);
    for (size_t i = 0; i < L; ++i) out[i] += (in[i] & mask);
}

template<class P>
inline void tlwe_add_plain_inplace(TFHEpp::TLWE<P>& out, uint8_t m) {
    out[P::n] += (m ? P::μ : typename P::T(0) - P::μ);
}

template<class P>
inline size_t zero_count_for_pk() {
    constexpr size_t LOG2Q = 8 * sizeof(typename P::T);
    return (P::n + 1) * LOG2Q + 128;
}

template<class P>
inline LwePublicKey<P> generate_public_key(const TFHEpp::SecretKey& sk) {
    LwePublicKey<P> pk;
    const size_t M = zero_count_for_pk<P>();
    pk.zeros.resize(M);
    for (size_t i = 0; i < M; ++i) {
        pk.zeros[i] = TFHEpp::tlweSymEncrypt<P>(0, P::α, sk.key.get<P>());
    }
    return pk;
}

template<class P, class URBG>
inline TFHEpp::TLWE<P> encrypt_public_key_bit(const LwePublicKey<P>& pk,
                                              uint8_t message,
                                              URBG& rng) {
    std::bernoulli_distribution coin(0.5);
    TFHEpp::TLWE<P> out{};
    std::fill(out.begin(), out.end(), typename P::T(0));
    for (const auto& enc0 : pk.zeros) {
        const uint32_t r = static_cast<uint32_t>(coin(rng));
        tlwe_add_selected_inplace<P>(out, enc0, r);
    }
    tlwe_add_plain_inplace<P>(out, static_cast<uint8_t>(message & 1u));
    return out;
}



template<class P>
struct SeededLwePublicKey {
    autohog::PublicSeed seed;
    std::vector<typename P::T> b_list;
};

template<class P>
inline SeededLwePublicKey<P> generate_compressed_public_key(const TFHEpp::SecretKey& sk,
                                                            autohog::PublicSeed seed = autohog::secure_public_seed()) {
    const auto& s = sk.key.get<P>();
    const size_t M = zero_count_for_pk<P>();
    SeededLwePublicKey<P> cpk;
    cpk.seed = seed;
    cpk.b_list.resize(M);
    autohog::SecureRng noise_rng;
    autohog::PublicSeedStream public_stream(seed, 0x4c574541u);
    for (size_t i = 0; i < M; ++i) {
        unsigned __int128 acc = 0;
        for (size_t j = 0; j < P::n; ++j) {
            const typename P::T aij = public_stream.template next_uniform<typename P::T>();
            acc += static_cast<unsigned __int128>(aij) * static_cast<unsigned __int128>(s[j]);
        }
        typename P::T inner = static_cast<typename P::T>(acc);
        std::normal_distribution<double> gauss(0.0, P::α);
        const typename P::T e_torus = torus_from_normalized_noise<typename P::T>(
            static_cast<long double>(gauss(noise_rng)));
        cpk.b_list[i] = inner + e_torus;
    }
    return cpk;
}

template<class P, class URBG>
inline TFHEpp::TLWE<P> encrypt_with_compressed_public_key_bit(const SeededLwePublicKey<P>& cpk,
                                                              uint8_t message,
                                                              URBG& rng_coin) {
    const size_t M = cpk.b_list.size();
    TFHEpp::TLWE<P> out{};
    std::fill(out.begin(), out.end(), typename P::T(0));
    std::bernoulli_distribution coin(0.5);
    autohog::PublicSeedStream public_stream(cpk.seed, 0x4c574541u);
    for (size_t i = 0; i < M; ++i) {
        const typename P::T mask = typename P::T(0) - typename P::T(coin(rng_coin) ? 1 : 0);
        for (size_t j = 0; j < P::n; ++j) {
            const typename P::T aij = public_stream.template next_uniform<typename P::T>();
            out[j] += (aij & mask);
        }
        out[P::n] += (cpk.b_list[i] & mask);
    }
    out[P::n] += (message & 1u) ? P::μ : (typename P::T(0) - P::μ);
    return out;
}



template<class P>
struct RlwePkPair {
    std::array<TFHEpp::Polynomial<P>, P::k> a;
    TFHEpp::Polynomial<P> b;
};

template<class P, class URBG>
inline TFHEpp::Polynomial<P> sample_uniform_poly(URBG& rng) {
    TFHEpp::Polynomial<P> poly{};
    std::uniform_int_distribution<typename P::T> U;
    for (auto &co : poly) co = U(rng);
    return poly;
}

template<class P, class URBG>
inline TFHEpp::Polynomial<P> sample_gaussian_poly(URBG& rng, double stddev = P::α) {
    TFHEpp::Polynomial<P> poly{};
    std::normal_distribution<double> gauss(0.0, stddev);
    for (auto &co : poly) {
        co = torus_from_normalized_noise<typename P::T>(
            static_cast<long double>(gauss(rng)));
    }
    return poly;
}

template<class P>
inline void get_trlwe_secret_polys(const TFHEpp::SecretKey& sk,
                                   std::array<TFHEpp::Polynomial<P>, P::k>& s_out) {
    const auto& key = sk.key.get<P>();
    for (size_t t = 0; t < P::k; ++t)
        for (size_t i = 0; i < P::n; ++i)
            s_out[t][i] = key[t * P::n + i];
}

template<class P, class URBG>
inline RlwePkPair<P> generate_rlwe_pk_pair(const TFHEpp::SecretKey& sk, URBG& rng) {
    RlwePkPair<P> pk{};
    for (size_t i = 0; i < P::k; ++i) pk.a[i] = sample_uniform_poly<P>(rng);
    pk.b = sample_gaussian_poly<P>(rng, P::α);
    std::array<TFHEpp::Polynomial<P>, P::k> s{};
    get_trlwe_secret_polys<P>(sk, s);
    TFHEpp::Polynomial<P> tmp{};
    for (size_t i = 0; i < P::k; ++i) {
        TFHEpp::PolyMul<P>(tmp, pk.a[i], s[i]);
        for (size_t j = 0; j < P::n; ++j) pk.b[j] += tmp[j];
    }
    return pk;
}

template<class P>
inline void assign_poly_to_trlwe_comp(TFHEpp::TRLWE<P>& out, size_t i,
                                      const TFHEpp::Polynomial<P>& coeff_poly) {
    using CompT = std::remove_reference_t<decltype(out[i])>;
    if constexpr (std::is_same_v<CompT, TFHEpp::Polynomial<P>>) out[i] = coeff_poly;
    else { TFHEpp::PolynomialInFD<P> fd{}; TFHEpp::TwistIFFT<P>(fd, coeff_poly); out[i] = fd; }
}

template<class P, class URBG>
inline TFHEpp::TRLWE<P> encrypt_bits_trlwe_with_rlwe_pair(const RlwePkPair<P>& pk,
                                                          const std::array<bool, P::n>& bits,
                                                          URBG& rng) {
    TFHEpp::TRLWE<P> out{};
    const auto v  = sample_gaussian_poly<P>(rng, P::α);
    const auto e0 = sample_gaussian_poly<P>(rng, P::α);
    std::array<TFHEpp::Polynomial<P>, P::k> e1{};
    for (size_t i = 0; i < P::k; ++i) e1[i] = sample_gaussian_poly<P>(rng, P::α);

    TFHEpp::Polynomial<P> m{};
    for (size_t j = 0; j < P::n; ++j) m[j] = bits[j] ? P::μ : (typename P::T(0) - P::μ);

    TFHEpp::Polynomial<P> tmp{}, ca{}, cb{};

    for (size_t i = 0; i < P::k; ++i) {
        TFHEpp::PolyMul<P>(tmp, pk.a[i], v);
        for (size_t j = 0; j < P::n; ++j) ca[j] = tmp[j] + e1[i][j];
        assign_poly_to_trlwe_comp<P>(out, i, ca);
    }
    TFHEpp::PolyMul<P>(cb, pk.b, v);
    for (size_t j = 0; j < P::n; ++j) cb[j] += e0[j] + m[j];
    assign_poly_to_trlwe_comp<P>(out, P::k, cb);

    return out;
}

template<class P>
struct SeededRlwePkPair {
    autohog::PublicSeed seed;
    TFHEpp::Polynomial<P> b;
};

template<class P>
inline void uniform_poly_from_seed(TFHEpp::Polynomial<P>& out,
                                   const autohog::PublicSeed& seed,
                                   std::uint64_t idx) {
    const std::uint64_t byte_offset =
        autohog::public_object_byte_offset<typename P::T>(idx, P::n);
    autohog::PublicSeedStream stream(seed, 0x524c5745u, byte_offset);
    for (auto& c : out) c = stream.template next_uniform<typename P::T>();
}

template<class P, class URBG>
inline SeededRlwePkPair<P> generate_rlwe_pk_pair_seeded(
    const TFHEpp::SecretKey& sk,
    URBG& rng,
    autohog::PublicSeed seed = autohog::secure_public_seed()) {
    SeededRlwePkPair<P> spk{};
    spk.seed = seed;
    // Every in-tree caller supplies autohog::SecureRng here. The public seed
    // is generated independently above and is never used to derive this noise.
    spk.b = sample_gaussian_poly<P>(rng, P::α);
    std::array<TFHEpp::Polynomial<P>, P::k> s{};
    get_trlwe_secret_polys<P>(sk, s);
    TFHEpp::Polynomial<P> ai{}, tmp{};
    for (std::size_t i = 0; i < P::k; ++i) {
        uniform_poly_from_seed<P>(ai, spk.seed, i);
        TFHEpp::PolyMul<P>(tmp, ai, s[i]);
        for (std::size_t j = 0; j < P::n; ++j) spk.b[j] += tmp[j];
    }
    return spk;
}

template<class P, class URBG>
inline TFHEpp::TRLWE<P> encrypt_bits_trlwe_with_seeded_rlwe_pair(const SeededRlwePkPair<P>& spk,
                                                                 const std::array<bool, P::n>& bits,
                                                                 URBG& rng) {
    TFHEpp::TRLWE<P> out{};
    TFHEpp::Polynomial<P> v  = sample_gaussian_poly<P>(rng, P::α);
    TFHEpp::Polynomial<P> e0 = sample_gaussian_poly<P>(rng, P::α);
    std::array<TFHEpp::Polynomial<P>, P::k> e1{};
    for (std::size_t i = 0; i < P::k; ++i) e1[i] = sample_gaussian_poly<P>(rng, P::α);
    TFHEpp::Polynomial<P> m{};
    for (std::size_t j = 0; j < P::n; ++j) m[j] = bits[j] ? P::μ : (typename P::T(0) - P::μ);
    TFHEpp::Polynomial<P> tmp{}, ca{}, cb{};
    TFHEpp::Polynomial<P> ai{};
    for (std::size_t i = 0; i < P::k; ++i) {
        uniform_poly_from_seed<P>(ai, spk.seed, i);
        TFHEpp::PolyMul<P>(tmp, ai, v);
        for (std::size_t j = 0; j < P::n; ++j) ca[j] = tmp[j] + e1[i][j];
        assign_poly_to_trlwe_comp<P>(out, i, ca);
    }
    TFHEpp::PolyMul<P>(cb, spk.b, v);
    for (std::size_t j = 0; j < P::n; ++j) cb[j] += e0[j] + m[j];
    assign_poly_to_trlwe_comp<P>(out, P::k, cb);
    return out;
}

} // namespace PKE
