#pragma once

#include <array>
#include <vector>
#include <cstdint>
#include <type_traits>
#include <cmath>
#include <limits>
#include <algorithm>

namespace driftms {



template <class U>
struct Bits {
    static constexpr int value =
        std::numeric_limits<std::make_unsigned_t<U>>::digits;
};

// ===================== round-to-grid (nearest) & round_error =====================
//


//


template <class U>
inline U round_to_grid_nearest(U input, int log_mod)
{
    static_assert(std::is_unsigned<U>::value, "round_to_grid_nearest: U must be unsigned torus type");
    const int BITS = Bits<U>::value;
    const int sh   = BITS - log_mod;
    if (sh <= 0) return input;

    // round-to-nearest: add half grid before shifting
    const U half = static_cast<U>(U{1} << (sh - 1));
    const U q    = static_cast<U>((input + half) >> sh);
    return static_cast<U>(q << sh);
}

// rounded = input + round_error  (mod 2^BITS)
template <class U>
inline U round_error(U input, int log_mod)
{
    static_assert(std::is_unsigned<U>::value, "round_error: U must be unsigned torus type");
    const U rounded = round_to_grid_nearest<U>(input, log_mod);
    return static_cast<U>(rounded - input); // ring subtraction
}


template <class U>
inline double round_error_f64(U input, int log_mod)
{
    static_assert(std::is_unsigned<U>::value, "round_error_f64: U must be unsigned torus type");
    using S = std::make_signed_t<U>;
    const U err = round_error<U>(input, log_mod);
    return static_cast<double>(static_cast<S>(err));
}



struct NoiseEstimation {
    double expectancy; // μ
    double variance;   // σ_mask^2
};


/// e_i = round_error(a_i), e_b = round_error(b)
/// μ = e_b - (1/2) Σ e_i
/// σ^2 = (1/4) Σ e_i^2
template <class MaskIter, class U>
inline NoiseEstimation
measure_modulus_switch_noise_expectancy_variance_for_binary_key(
    MaskIter masks_begin,
    MaskIter masks_end,
    U        body,
    int      log_mod)
{
    long double sum_mask_errors   = 0.0L;
    long double sum_square_errors = 0.0L;

    for (auto it = masks_begin; it != masks_end; ++it) {
        const auto ai = static_cast<std::make_unsigned_t<U>>(*it);
        const long double e = static_cast<long double>(round_error_f64(ai, log_mod));
        sum_mask_errors   += e;
        sum_square_errors += e * e;
    }

    const long double eb =
        static_cast<long double>(round_error_f64(static_cast<std::make_unsigned_t<U>>(body), log_mod));

    const long double mu  = eb - 0.5L * sum_mask_errors;
    const long double var = 0.25L * sum_square_errors;

    return NoiseEstimation{static_cast<double>(mu), static_cast<double>(var)};
}

/// measure = |μ| + r_sigma * sqrt(σ_mask^2 + input_variance)
///


template <class MaskIter, class U>
inline double
measure_modulus_switch_noise_estimation_for_binary_key(
    double   r_sigma,
    double   input_variance,
    int      log_mod,
    MaskIter masks_begin,
    MaskIter masks_end,
    U        body)
{
    const auto ne =
        measure_modulus_switch_noise_expectancy_variance_for_binary_key<MaskIter, U>(
            masks_begin, masks_end, body, log_mod);

    const double sigma_total =
        std::sqrt(std::max(0.0, ne.variance) + std::max(0.0, input_variance));

    return std::abs(ne.expectancy) + r_sigma * sigma_total;
}



enum class CandidateKind {
    NoAddition = 0,
    AddZero    = 1,
};

struct Candidate {
    CandidateKind kind{CandidateKind::NoAddition};
    int           index{-1};
};

enum class CandidateResultKind {
    
    SatisfyingBound,
    BestNotSatisfyingBound
};

struct CandidateResult {
    CandidateResultKind kind{CandidateResultKind::BestNotSatisfyingBound};
    Candidate           candidate{};
    double              best_measure{std::numeric_limits<double>::infinity()}; 
};



template <class P>
struct DriftMsKey {
    using Torus    = typename P::T;
    static constexpr std::size_t LWE_DIM =
        static_cast<std::size_t>(P::k) * static_cast<std::size_t>(P::n);
    using TLWEType = std::array<Torus, LWE_DIM + 1>;

    std::vector<TLWEType> zeros; // encryptions of 0

    double r_sigma{0.0};
    double bound{0.0};          
    double input_variance{0.0};
    int    log_mod{0};
};


template <class P>
inline void tlwe_add_assign(typename DriftMsKey<P>::TLWEType &dst,
                            const typename DriftMsKey<P>::TLWEType &src)
{
    using Torus = typename P::T;
    constexpr std::size_t LWE_SIZE = DriftMsKey<P>::LWE_DIM + 1;
    for (std::size_t i = 0; i < LWE_SIZE; ++i) {
        dst[i] = static_cast<Torus>(dst[i] + src[i]);
    }
}


//




template <class P>
inline CandidateResult
choose_candidate_to_minimize_modulus_switch_noise_for_binary_key(
    const typename DriftMsKey<P>::TLWEType &lwe,
    const DriftMsKey<P>                    &dk)
{
    using Torus = typename P::T;
    constexpr std::size_t LWE_DIM = DriftMsKey<P>::LWE_DIM;

    
    struct MaskIterBase {
        const Torus *a;
        std::size_t  pos{0};
        bool operator!=(const MaskIterBase &o) const { return pos != o.pos; }
        void operator++() { ++pos; }
        Torus operator*() const { return a[pos]; }
    };

    MaskIterBase base_beg{lwe.data(), 0};
    MaskIterBase base_end{lwe.data(), LWE_DIM};
    const Torus body = lwe[LWE_DIM];

    const double base_measure =
        measure_modulus_switch_noise_estimation_for_binary_key<MaskIterBase, Torus>(
            dk.r_sigma, dk.input_variance, dk.log_mod, base_beg, base_end, body);

    CandidateResult res;
    res.kind = CandidateResultKind::BestNotSatisfyingBound;
    res.candidate = Candidate{CandidateKind::NoAddition, -1};
    res.best_measure = base_measure;

    
    if (dk.zeros.empty()) {
        return res;
    }

    
    for (std::size_t j = 0; j < dk.zeros.size(); ++j) {
        const auto &z = dk.zeros[j];

        struct MaskIterSum {
            const Torus *a;
            const Torus *z;
            std::size_t  pos{0};
            bool operator!=(const MaskIterSum &o) const { return pos != o.pos; }
            void operator++() { ++pos; }
            Torus operator*() const { return static_cast<Torus>(a[pos] + z[pos]); }
        };

        MaskIterSum sum_beg{lwe.data(), z.data(), 0};
        MaskIterSum sum_end{lwe.data(), z.data(), LWE_DIM};
        const Torus body_sum = static_cast<Torus>(lwe[LWE_DIM] + z[LWE_DIM]);

        const double m =
            measure_modulus_switch_noise_estimation_for_binary_key<MaskIterSum, Torus>(
                dk.r_sigma, dk.input_variance, dk.log_mod, sum_beg, sum_end, body_sum);

        if (m < res.best_measure) {
            res.best_measure = m;
            res.candidate.kind  = CandidateKind::AddZero;
            res.candidate.index = static_cast<int>(j);
        }
    }

    return res;
}



template <class P>
inline void
improve_lwe_ciphertext_modulus_switch_noise_for_binary_key(
    typename DriftMsKey<P>::TLWEType &lwe,
    const DriftMsKey<P>              &dk)
{
    const auto cr =
        choose_candidate_to_minimize_modulus_switch_noise_for_binary_key<P>(lwe, dk);

    const Candidate &cand = cr.candidate;
    if (cand.kind == CandidateKind::AddZero && cand.index >= 0 &&
        static_cast<std::size_t>(cand.index) < dk.zeros.size()) {
        tlwe_add_assign<P>(lwe, dk.zeros[static_cast<std::size_t>(cand.index)]);
    }
}

} // namespace driftms



//


//
//     body +=  sum_i ( round(a_i) - a_i ) / 2  -  half_case
//
//   half_case = 2^(BITS - log_mod - 1)
//





namespace centered {

template <class U>
using Signed = std::make_signed_t<U>;

template <class U>
inline Signed<U> to_signed(U x) {
    return static_cast<Signed<U>>(x);
}

template <class U>
inline U to_unsigned(Signed<U> x) {
    return static_cast<U>(x);
}

template <class U>
inline U round_to_grid_nearest(U input, int log_mod)
{
    static_assert(std::is_unsigned<U>::value);
    const int BITS = driftms::Bits<U>::value;
    const int sh = BITS - log_mod;
    if (sh <= 0) return input;

    
    const U half = static_cast<U>(U{1} << (sh - 1));
    const U q = static_cast<U>((input + half) >> sh);
    return static_cast<U>(q << sh);
}

template <class P>
inline void apply_centered_correction_to_body(
    typename driftms::DriftMsKey<P>::TLWEType &lwe,
    int log_mod)
{
    using Torus = typename P::T;
    static_assert(std::is_unsigned<Torus>::value,
                  "apply_centered_correction_to_body: Torus must be unsigned");

    constexpr std::size_t LWE_DIM = driftms::DriftMsKey<P>::LWE_DIM;
    const int BITS = driftms::Bits<Torus>::value;

    
    if (log_mod <= 0 || log_mod >= BITS) {
        return;
    }

    Torus sum_half_mask_round_errors = Torus{0};
    Signed<Torus> sum_halving_errors_doubled = Signed<Torus>{0};

    
    for (std::size_t i = 0; i < LWE_DIM; ++i) {
        const Torus a = lwe[i];

        // error = round(a) - a  (ring)
        const Torus rounded = round_to_grid_nearest<Torus>(a, log_mod);
        const Torus error_u = static_cast<Torus>(rounded - a);
        const Signed<Torus> signed_error = to_signed<Torus>(error_u);

        // half_error = signed_error / 2
        const Signed<Torus> half_error = signed_error / Signed<Torus>{2};

        // doubled halving error: 2*half_error - signed_error ∈ {-1,0,1}
        const Signed<Torus> halving_error_doubled = Signed<Torus>{2} * half_error - signed_error;

        sum_half_mask_round_errors = static_cast<Torus>(sum_half_mask_round_errors + to_unsigned<Torus>(half_error));

        sum_halving_errors_doubled += halving_error_doubled;
    }

    
    const Signed<Torus> sum_halving_errors_signed =
        sum_halving_errors_doubled / Signed<Torus>{2};
    const Torus sum_halving_errors_u =
        to_unsigned<Torus>(sum_halving_errors_signed);

    sum_half_mask_round_errors =
        static_cast<Torus>(sum_half_mask_round_errors - sum_halving_errors_u);

    // half_case = 2^(BITS - log_mod - 1)
    const Torus half_case =
        static_cast<Torus>(Torus{1} << (BITS - log_mod - 1));

    // body += sum(mask_round_error/2) - half_case
    lwe[LWE_DIM] =
        static_cast<Torus>(lwe[LWE_DIM] +
                           sum_half_mask_round_errors -
                           half_case);
}

} // namespace centered


namespace centered_tfhepp {



template <class U>
using Signed = std::make_signed_t<U>;

template <class U>
static inline Signed<U> to_signed(U x) {
    return static_cast<Signed<U>>(x);
}

template <class U>
static inline U to_unsigned(Signed<U> x) {
    return static_cast<U>(x);
}


//

//   roundoffset = 1 << (digits - 2 - targetP::nbit + bitwidth)

//


template <class Torus>
static inline Torus round_to_grid_nearest_br(Torus x, int sh)
{
    // nearest: (x + 2^{sh-1}) >> sh << sh
    
    const Torus half = Torus{1} << (sh - 1);
    const Torus q    = (x + half) >> sh;
    return static_cast<Torus>(q << sh);
}


//




//

//   e_i = round_nearest(a_i) - a_i    (signed)

//




//
template <class Pbk, uint32_t num_out = 1>
inline void apply_body_correction_like_br(
    TFHEpp::TLWE<typename Pbk::domainP> &tlwe)
{
    using Torus = typename Pbk::domainP::T;
    static_assert(std::is_unsigned<Torus>::value,
                  "apply_body_correction_like_br: Torus must be unsigned");

    constexpr uint32_t bitwidth = TFHEpp::bits_needed<num_out - 1>();
    constexpr int digits = std::numeric_limits<Torus>::digits;

    
    // sh = digits - 1 - targetP::nbit + bitwidth
    constexpr int sh = digits - 1 - (int)Pbk::targetP::nbit + (int)bitwidth;

    if constexpr (sh <= 0 || sh >= digits) {
        
        return;
    }

    constexpr int LWE_DIM = (int)(Pbk::domainP::k * Pbk::domainP::n);
    constexpr int BODY_IDX = LWE_DIM;

    
    __int128 sum_half_e = 0;

    
    Signed<Torus> sum_halving_errors_doubled = 0;

    for (int i = 0; i < LWE_DIM; ++i) {
        const Torus a = tlwe[i];

        // e = round(a) - a  (signed)
        const Torus rounded = round_to_grid_nearest_br<Torus>(a, sh);
        const Signed<Torus> e = to_signed<Torus>(static_cast<Torus>(rounded - a));

        // half_e = e / 2  (trunc toward 0)
        const Signed<Torus> half_e = e / Signed<Torus>{2};

        // 2*half_e - e ∈ {-1,0,1}
        const Signed<Torus> halving_err2 = Signed<Torus>{2} * half_e - e;

        sum_half_e += (__int128)half_e;
        sum_halving_errors_doubled += halving_err2;
    }

    
    sum_half_e -= (__int128)(sum_halving_errors_doubled / Signed<Torus>{2});

    // body += Σ e_i/2
    
    const Torus add_u = static_cast<Torus>(sum_half_e); // wrap into ring
    tlwe[BODY_IDX] = static_cast<Torus>(tlwe[BODY_IDX] + add_u);
}


//

//   centered_tfhepp::apply_before_blind_rotate<lvl02param>(tlwe0);
//   centered_tfhepp::apply_before_blind_rotate_many<lvl02param, 8>(tlwe0);
//

//
template <class Pbk>
inline void apply_before_blind_rotate(
    TFHEpp::TLWE<typename Pbk::domainP> &tlwe)
{
    apply_body_correction_like_br<Pbk, 1>(tlwe);
}

template <class Pbk, uint32_t num_out>
inline void apply_before_blind_rotate_many(
    TFHEpp::TLWE<typename Pbk::domainP> &tlwe)
{
    apply_body_correction_like_br<Pbk, num_out>(tlwe);
}

} // namespace centered_tfhepp