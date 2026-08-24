#include "ZKProof.h"
#include "HomPKE.h"
#include "zk_backend.h"
#include <chrono>
#include <cstdlib>
#include <openssl/evp.h>
#include <cctype>
#include <cstring>
#include <random>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <unordered_map>
#include <omp.h>

namespace ZKProof {

template<class T>
class ScrubbedValue {
public:
    ScrubbedValue() = default;
    ~ScrubbedValue() { autohog::secure_erase(&value_, sizeof(value_)); }

    ScrubbedValue(const ScrubbedValue&) = delete;
    ScrubbedValue& operator=(const ScrubbedValue&) = delete;

    T& get() { return value_; }
    const T& get() const { return value_; }

private:
    T value_{};
};

static size_t default_statement_k_lvl2() {
    return 320;
}

static constexpr size_t kG1CompressedBytes = autohog::zk::config::g1_serialized_size;
static constexpr size_t kG2CompressedBytes = autohog::zk::config::g2_serialized_size;
static constexpr size_t kProofSerializedBytes = autohog::zk::config::proof_serialized_size;
static constexpr char kCrsMagicCurrent[4] = {'Z', 'K', 'C', '3'};
static constexpr char kCrsMagicLegacy[4] = {'Z', 'K', 'C', '2'};

using ProfileClock = std::chrono::steady_clock;

static bool zk_profile_enabled() {
    static const bool enabled = [] {
        const char* env = std::getenv("AUTOHOG_ZK_PROFILE");
        if (!env || env[0] == '\0') return false;
        return std::strcmp(env, "0") != 0 &&
               std::strcmp(env, "false") != 0 &&
               std::strcmp(env, "FALSE") != 0;
    }();
    return enabled;
}

static double profile_elapsed_ms(ProfileClock::time_point start,
                                 ProfileClock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

static void emit_profile_summary(const char* label,
                                 const std::vector<std::pair<const char*, double>>& parts) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2);
    oss << "[" << label << " profile]";
    double total_ms = 0.0;
    for (const auto& part : parts) {
        total_ms += part.second;
        oss << " " << part.first << "=" << part.second << "ms";
    }
    oss << " total=" << total_ms << "ms";
    std::cerr << oss.str() << '\n';
}

#if !defined(AUTOHOG_ZK_BACKEND_mcl)
static const blst_fp12& pairing_identity_gt() {
    static const blst_fp12 identity = [] {
        const blst_p1_affine* g1_gen = blst_p1_affine_generator();
        const blst_p2_affine* g2_gen = blst_p2_affine_generator();

        blst_p1 neg_g1;
        blst_p1_from_affine(&neg_g1, g1_gen);
        blst_p1_cneg(&neg_g1, 1);
        blst_p1_affine neg_g1_aff;
        blst_p1_to_affine(&neg_g1_aff, &neg_g1);

        blst_fp12 ml_a, ml_b, identity_local;
        blst_miller_loop(&ml_a, g2_gen, g1_gen);
        blst_miller_loop(&ml_b, g2_gen, &neg_g1_aff);
        blst_fp12_mul(&ml_a, &ml_a, &ml_b);
        blst_final_exp(&identity_local, &ml_a);
        return identity_local;
    }();
    return identity;
}
#endif

static void shake256_hash(
    const uint8_t* input,
    size_t input_len,
    uint8_t* output,
    size_t output_len
) {
    if ((input == nullptr && input_len != 0) ||
        (output == nullptr && output_len != 0)) {
        throw std::invalid_argument("SHAKE256 received an invalid buffer");
    }
    std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> ctx(
        EVP_MD_CTX_new(), &EVP_MD_CTX_free);
    if (!ctx) throw std::runtime_error("Failed to create EVP context");
    if (EVP_DigestInit_ex(ctx.get(), EVP_shake256(), nullptr) != 1 ||
        (input_len != 0 && EVP_DigestUpdate(ctx.get(), input, input_len) != 1) ||
        EVP_DigestFinalXOF(ctx.get(), output, output_len) != 1) {
        throw std::runtime_error("SHAKE256 operation failed");
    }
}

static __uint128_t decode_q(uint64_t q) {
    return (q == 0) ? ((__uint128_t)1 << 64) : (__uint128_t)q;
}

static blst_fr fr_zero() {
    blst_scalar zero_scalar;
    std::memset(&zero_scalar, 0, sizeof(zero_scalar));
    blst_fr r;
    blst_fr_from_scalar(&r, &zero_scalar);
    return r;
}


static blst_fr fr_add(const blst_fr& a, const blst_fr& b) {
    blst_fr r; blst_fr_add(&r, &a, &b); return r;
}
static blst_fr fr_sub(const blst_fr& a, const blst_fr& b) {
    blst_fr r; blst_fr_sub(&r, &a, &b); return r;
}
static blst_fr fr_mul(const blst_fr& a, const blst_fr& b) {
    blst_fr r; blst_fr_mul(&r, &a, &b); return r;
}
static blst_fr fr_neg(const blst_fr& a) {
    blst_fr r; blst_fr_cneg(&r, &a, true); return r;
}

static blst_scalar fr_to_scalar(const blst_fr& a) {
    blst_scalar s; blst_scalar_from_fr(&s, &a); return s;
}

static bool fr_is_zero(const blst_fr& a) {
    blst_scalar s = fr_to_scalar(a);
    for (int i = 0; i < 32; i++) if (s.b[i]) return false;
    return true;
}

static bool scalar_is_zero(const blst_scalar& s) {
    for (int i = 0; i < 32; ++i) {
        if (s.b[i] != 0) return false;
    }
    return true;
}

static blst_fr fr_from_u64(uint64_t v) {
    blst_scalar s; memset(&s, 0, sizeof(s));
    memcpy(s.b, &v, 8);
    blst_fr r; blst_fr_from_scalar(&r, &s); return r;
}

blst_fr fr_from_i64(int64_t val) {
    blst_fr r;
    uint64_t abs_val = (val < 0)
        ? (uint64_t(0) - static_cast<uint64_t>(val))
        : static_cast<uint64_t>(val);
    uint64_t vals[4] = {abs_val, 0, 0, 0};
    blst_fr_from_uint64(&r, vals);
    if (val < 0) {
        return fr_neg(r); 
    }
    return r;
}
static blst_fr fr_from_u128(__uint128_t v) {
    blst_scalar s; memset(&s, 0, sizeof(s));
    uint64_t lo = static_cast<uint64_t>(v);
    uint64_t hi = static_cast<uint64_t>(v >> 64);
    memcpy(s.b,     &lo, 8);
    memcpy(s.b + 8, &hi, 8);
    blst_fr r; blst_fr_from_scalar(&r, &s); return r;
}

static blst_fr fr_inv(const blst_fr& a) {
    blst_fr r;
    blst_fr_inverse(&r, &a);
    return r;
}

static blst_fr fr_pow_u64(blst_fr base, uint64_t exp) {
    blst_fr result = fr_from_u64(1);
    while (exp > 0) {
        if (exp & 1) {
            result = fr_mul(result, base);
        }
        exp >>= 1;
        if (exp == 0) {
            break;
        }
        base = fr_mul(base, base);
    }
    return result;
}

static std::vector<uint8_t> fr_to_bytes(const blst_fr& a) {
    blst_scalar s = fr_to_scalar(a);
    return std::vector<uint8_t>(s.b, s.b + 32);
}

static std::vector<uint8_t> u64_to_le_bytes(uint64_t v) {
    std::vector<uint8_t> b(8); memcpy(b.data(), &v, 8); return b;
}

static std::vector<uint8_t> i64_vec_to_bytes(const std::vector<int64_t>& v) {
    std::vector<uint8_t> out(v.size() * 8);
    for (size_t i = 0; i < v.size(); i++) memcpy(out.data() + i * 8, &v[i], 8);
    return out;
}

static std::vector<uint8_t> fr_vec_to_bytes(const std::vector<blst_fr>& v) {
    std::vector<uint8_t> out;
    for (const auto& f : v) {
        auto b = fr_to_bytes(f);
        out.insert(out.end(), b.begin(), b.end());
    }
    return out;
}

static std::vector<uint8_t> p1_to_bytes(const blst_p1& p) {
    blst_p1_affine aff;
    blst_p1_to_affine(&aff, &p);
    std::vector<uint8_t> c(kG1CompressedBytes);
    blst_p1_affine_compress(c.data(), &aff);
    return c;
}

static std::vector<uint8_t> p2_to_bytes(const blst_p2& p) {
    blst_p2_affine aff;
    blst_p2_to_affine(&aff, &p);
    std::vector<uint8_t> c(kG2CompressedBytes);
    blst_p2_affine_compress(c.data(), &aff);
    return c;
}

static bool g1_affine_usable(const blst_p1_affine* p) {
    return p && !blst_p1_affine_is_inf(p) &&
           blst_p1_affine_on_curve(p) &&
           blst_p1_affine_in_g1(p);
}

static bool g2_affine_usable(const blst_p2_affine* p) {
    return p && !blst_p2_affine_is_inf(p) &&
           blst_p2_affine_on_curve(p) &&
           blst_p2_affine_in_g2(p);
}

static bool proof_points_usable(const Proof& proof) {
    return g2_affine_usable(static_cast<const blst_p2_affine*>(proof.c_hat_ptr)) &&
           g1_affine_usable(static_cast<const blst_p1_affine*>(proof.c_y_ptr)) &&
           g1_affine_usable(static_cast<const blst_p1_affine*>(proof.pi_ptr));
}

static blst_p1 p1_zero() {
    blst_scalar zero_scalar;
    std::memset(&zero_scalar, 0, sizeof(zero_scalar));
    blst_p1 zero;
    blst_p1_from_affine(&zero, blst_p1_affine_generator());
    blst_p1_mult(&zero, &zero, zero_scalar.b, 256);
    return zero;
}

static blst_p2 p2_zero() {
    blst_scalar zero_scalar;
    std::memset(&zero_scalar, 0, sizeof(zero_scalar));
    blst_p2 zero;
    blst_p2_from_affine(&zero, blst_p2_affine_generator());
    blst_p2_mult(&zero, &zero, zero_scalar.b, 256);
    return zero;
}

static void add_p1_term(blst_p1& acc, bool& has_acc, const blst_p1& term) {
    if (has_acc) {
        blst_p1_add(&acc, &acc, &term);
    } else {
        acc = term;
        has_acc = true;
    }
}

static void add_p2_term(blst_p2& acc, bool& has_acc, const blst_p2& term) {
    if (has_acc) {
        blst_p2_add(&acc, &acc, &term);
    } else {
        acc = term;
        has_acc = true;
    }
}

static blst_p1 combine_p1_partials(const std::vector<blst_p1>& partials,
                                   const std::vector<uint8_t>& has_partial) {
    blst_p1 result = p1_zero();
    bool has_result = false;
    for (size_t i = 0; i < partials.size(); ++i) {
        if (!has_partial[i]) continue;
        add_p1_term(result, has_result, partials[i]);
    }
    return has_result ? result : p1_zero();
}

static blst_p2 combine_p2_partials(const std::vector<blst_p2>& partials,
                                   const std::vector<uint8_t>& has_partial) {
    blst_p2 result = p2_zero();
    bool has_result = false;
    for (size_t i = 0; i < partials.size(); ++i) {
        if (!has_partial[i]) continue;
        add_p2_term(result, has_result, partials[i]);
    }
    return has_result ? result : p2_zero();
}

static size_t msm_parallel_chunks(size_t npoints);

static std::vector<limb_t>& p1_pippenger_scratch(size_t point_count) {
    static thread_local std::vector<limb_t> scratch;
    const size_t scratch_bytes = blst_p1s_mult_pippenger_scratch_sizeof(point_count);
    const size_t scratch_limbs = (scratch_bytes + sizeof(limb_t) - 1) / sizeof(limb_t);
    if (scratch.size() < scratch_limbs) {
        scratch.resize(scratch_limbs);
    }
    return scratch;
}

static std::vector<limb_t>& p2_pippenger_scratch(size_t point_count) {
    static thread_local std::vector<limb_t> scratch;
    const size_t scratch_bytes = blst_p2s_mult_pippenger_scratch_sizeof(point_count);
    const size_t scratch_limbs = (scratch_bytes + sizeof(limb_t) - 1) / sizeof(limb_t);
    if (scratch.size() < scratch_limbs) {
        scratch.resize(scratch_limbs);
    }
    return scratch;
}

static blst_p2 p2_add_affines(const std::vector<const blst_p2_affine*>& points) {
    if (points.empty()) {
        return p2_zero();
    }

    const size_t chunks = msm_parallel_chunks(points.size());
    if (chunks == 1) {
        blst_p2 result;
        blst_p2s_add(&result, points.data(), points.size());
        return result;
    }

    std::vector<blst_p2> partials(chunks);
    std::vector<uint8_t> has_partial(chunks, 0);
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t chunk = 0; chunk < static_cast<ptrdiff_t>(chunks); ++chunk) {
        const size_t begin = points.size() * static_cast<size_t>(chunk) / chunks;
        const size_t end = points.size() * static_cast<size_t>(chunk + 1) / chunks;
        if (begin >= end) continue;

        blst_p2s_add(
            &partials[static_cast<size_t>(chunk)],
            points.data() + begin,
            end - begin
        );
        has_partial[static_cast<size_t>(chunk)] = 1;
    }

    return combine_p2_partials(partials, has_partial);
}

static size_t msm_parallel_chunks(size_t npoints) {
    if (npoints == 0 || omp_in_parallel()) {
        return 1;
    }

    const size_t max_threads = static_cast<size_t>(std::max(1, omp_get_max_threads()));
    const size_t target_chunk_size = 2048;
    size_t chunks = (npoints + target_chunk_size - 1) / target_chunk_size;
    chunks = std::max<size_t>(1, std::min(chunks, max_threads));
    return chunks;
}

static size_t power_chain_parallel_chunks(size_t count) {
    if (count == 0 || omp_in_parallel()) {
        return 1;
    }
    const size_t max_threads = static_cast<size_t>(std::max(1, omp_get_max_threads()));
    const size_t target_chunk_size = 2048;
    size_t chunks = (count + target_chunk_size - 1) / target_chunk_size;
    chunks = std::max<size_t>(1, std::min(chunks, max_threads));
    return chunks;
}

static blst_p1 p1_mult_pippenger(
    const std::vector<const blst_p1_affine*>& points,
    const std::vector<blst_scalar>& scalars
) {
    if (points.empty()) {
        return p1_zero();
    }

    std::vector<const unsigned char*> scalar_ptrs(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        scalar_ptrs[i] = scalars[i].b;
    }

    const size_t chunks = msm_parallel_chunks(points.size());
    if (chunks == 1) {
        auto& scratch = p1_pippenger_scratch(points.size());
        blst_p1 result;
        blst_p1s_mult_pippenger(
            &result,
            points.data(),
            points.size(),
            scalar_ptrs.data(),
            256,
            scratch.data()
        );
        return result;
    }

    std::vector<blst_p1> partials(chunks);
    std::vector<uint8_t> has_partial(chunks, 0);
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t chunk = 0; chunk < static_cast<ptrdiff_t>(chunks); ++chunk) {
        const size_t begin = points.size() * static_cast<size_t>(chunk) / chunks;
        const size_t end = points.size() * static_cast<size_t>(chunk + 1) / chunks;
        if (begin >= end) continue;

        const size_t chunk_count = end - begin;
        auto& scratch = p1_pippenger_scratch(chunk_count);
        blst_p1s_mult_pippenger(
            &partials[static_cast<size_t>(chunk)],
            points.data() + begin,
            chunk_count,
            scalar_ptrs.data() + begin,
            256,
            scratch.data()
        );
        has_partial[static_cast<size_t>(chunk)] = 1;
    }

    return combine_p1_partials(partials, has_partial);
}

static blst_p2 p2_mult_pippenger(
    const std::vector<const blst_p2_affine*>& points,
    const std::vector<blst_scalar>& scalars
) {
    if (points.empty()) {
        return p2_zero();
    }

    std::vector<const unsigned char*> scalar_ptrs(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        scalar_ptrs[i] = scalars[i].b;
    }

    const size_t chunks = msm_parallel_chunks(points.size());
    if (chunks == 1) {
        auto& scratch = p2_pippenger_scratch(points.size());
        blst_p2 result;
        blst_p2s_mult_pippenger(
            &result,
            points.data(),
            points.size(),
            scalar_ptrs.data(),
            256,
            scratch.data()
        );
        return result;
    }

    std::vector<blst_p2> partials(chunks);
    std::vector<uint8_t> has_partial(chunks, 0);
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t chunk = 0; chunk < static_cast<ptrdiff_t>(chunks); ++chunk) {
        const size_t begin = points.size() * static_cast<size_t>(chunk) / chunks;
        const size_t end = points.size() * static_cast<size_t>(chunk + 1) / chunks;
        if (begin >= end) continue;

        const size_t chunk_count = end - begin;
        auto& scratch = p2_pippenger_scratch(chunk_count);
        blst_p2s_mult_pippenger(
            &partials[static_cast<size_t>(chunk)],
            points.data() + begin,
            chunk_count,
            scalar_ptrs.data() + begin,
            256,
            scratch.data()
        );
        has_partial[static_cast<size_t>(chunk)] = 1;
    }

    return combine_p2_partials(partials, has_partial);
}

static blst_p1 p1_msm_direct(const blst_p1_affine* bases,
                             const blst_fr* scalars,
                             size_t count) {
    std::vector<const blst_p1_affine*> points(count);
    std::vector<blst_scalar> scalar_storage(count);

    #pragma omp parallel for schedule(static) if (count >= 4096)
    for (size_t i = 0; i < count; ++i) {
        points[i] = &bases[i];
        scalar_storage[i] = fr_to_scalar(scalars[i]);
    }

    return p1_mult_pippenger(points, scalar_storage);
}

static blst_p1 p1_msm_reverse_y_w(const blst_p1_affine* g_list,
                                  size_t n,
                                  const std::vector<blst_fr>& y,
                                  const std::vector<bool>& w,
                                  size_t count) {
    std::vector<const blst_p1_affine*> points;
    std::vector<blst_scalar> scalar_storage;
    points.reserve(count);
    scalar_storage.reserve(count);

    for (size_t j = 0; j < count; ++j) {
        if (!w[j]) continue;
        points.push_back(&g_list[n - 1 - j]);
        scalar_storage.push_back(fr_to_scalar(y[j]));
    }

    return p1_mult_pippenger(points, scalar_storage);
}

static blst_p1 p1_msm_verify_term0(const blst_p1_affine* g_list,
                                   size_t n,
                                   size_t big_d,
                                   const std::vector<blst_fr>& y,
                                   const std::vector<blst_fr>& t_vec,
                                   const std::vector<blst_fr>& a_theta,
                                   const blst_fr& delta_eq,
                                   const blst_fr& delta_y,
                                   const blst_fr& delta_theta) {
    std::vector<const blst_p1_affine*> points;
    std::vector<blst_scalar> scalar_storage;
    points.reserve(n);
    scalar_storage.reserve(n);

    for (size_t i = 1; i <= n; ++i) {
        blst_fr yi = y[i - 1];
        blst_fr ti = t_vec[i - 1];
        blst_fr factor = fr_mul(fr_sub(fr_mul(delta_eq, ti), delta_y), yi);
        if (i <= big_d) {
            factor = fr_add(factor, fr_mul(delta_theta, a_theta[i - 1]));
        }
        blst_scalar scalar = fr_to_scalar(factor);
        if (scalar_is_zero(scalar)) continue;
        points.push_back(&g_list[n - i]);
        scalar_storage.push_back(scalar);
    }

    return p1_mult_pippenger(points, scalar_storage);
}

static blst_p2 p2_msm_t_vec(const blst_p2_affine* g_hat_list,
                            size_t n,
                            const std::vector<blst_fr>& t_vec,
                            const blst_fr& delta_eq) {
    std::vector<const blst_p2_affine*> points;
    std::vector<blst_scalar> scalar_storage;
    points.reserve(n);
    scalar_storage.reserve(n);

    for (size_t i = 1; i <= n; ++i) {
        blst_fr coef = fr_mul(delta_eq, t_vec[i - 1]);
        blst_scalar scalar = fr_to_scalar(coef);
        if (scalar_is_zero(scalar)) continue;
        points.push_back(&g_hat_list[i - 1]);
        scalar_storage.push_back(scalar);
    }

    return p2_mult_pippenger(points, scalar_storage);
}

template<typename P>
static bool public_params_dimensions_usable(const PublicParams<P>& params) {
    if (params.d != P::n || params.k == 0 || params.k > params.d ||
        params.b == 0 || params.b_r == 0 || params.t == 0 ||
        params.msbs_zero_padding_bit_count >= 64) {
        return false;
    }

    const uint64_t effective_t = params.t >> params.msbs_zero_padding_bit_count;
    if (effective_t == 0) {
        return false;
    }

    const auto ilog2_floor = [](uint64_t value) {
        size_t result = 0;
        while (value > 1) {
            value >>= 1;
            ++result;
        }
        return result;
    };
    const auto checked_add = [](size_t a, size_t b, size_t& out) {
        if (a > (std::numeric_limits<size_t>::max)() - b) return false;
        out = a + b;
        return true;
    };
    const auto checked_mul = [](size_t a, size_t b, size_t& out) {
        if (a != 0 && b > (std::numeric_limits<size_t>::max)() / a) return false;
        out = a * b;
        return true;
    };

    const size_t d = static_cast<size_t>(params.d);
    const size_t k = static_cast<size_t>(params.k);
    const size_t step_t = ilog2_floor(effective_t);
    const size_t step_b = 1 + ilog2_floor(params.b);
    const size_t step_br = 1 + ilog2_floor(params.b_r);

    size_t k_t = 0;
    size_t d_plus_k = 0;
    size_t bound_steps = 0;
    size_t bounds_part = 0;
    size_t expected_n = 0;
    if (!checked_mul(k, step_t, k_t) ||
        !checked_add(d, k, d_plus_k) ||
        !checked_add(step_b, step_br, bound_steps) ||
        !checked_mul(d_plus_k, bound_steps, bounds_part) ||
        !checked_add(d, k_t, expected_n) ||
        !checked_add(expected_n, bounds_part, expected_n) ||
        !checked_add(expected_n, 1, expected_n)) {
        return false;
    }

    return params.n == expected_n && params.n >= 3 &&
           params.n <= (std::numeric_limits<size_t>::max)() / 2;
}

template<typename P>
static bool public_params_basic_usable(const PublicParams<P>& params) {
    return params.g_list_ptr && params.g_hat_list_ptr &&
           public_params_dimensions_usable(params);
}

static bool pairing_equal(const blst_p1_affine* p_left,
                          const blst_p2_affine* q_left,
                          const blst_p1_affine* p_right,
                          const blst_p2_affine* q_right) {
    blst_fp12 ml_left, ml_right, gt_left, gt_right;
    blst_miller_loop(&ml_left, q_left, p_left);
    blst_miller_loop(&ml_right, q_right, p_right);
    blst_final_exp(&gt_left, &ml_left);
    blst_final_exp(&gt_right, &ml_right);
    return blst_fp12_is_equal(&gt_left, &gt_right);
}

template<typename P>
static std::vector<uint8_t> crs_validation_transcript(const PublicParams<P>& params) {
    const auto* g_list = static_cast<const blst_p1_affine*>(params.g_list_ptr);
    const auto* g_hat_list = static_cast<const blst_p2_affine*>(params.g_hat_list_ptr);
    static constexpr char label[] = "AUTOHOG-CRS-POWER-CHAIN-CHECK-V1";

    std::vector<uint8_t> transcript;
    transcript.reserve(
        sizeof(label) - 1 + 8 * 8 + 8 * 16 +
        static_cast<size_t>(2 * params.n) * kG1CompressedBytes +
        static_cast<size_t>(params.n) * kG2CompressedBytes);
    transcript.insert(transcript.end(), label, label + sizeof(label) - 1);

    const auto append_u64 = [&](uint64_t value) {
        std::array<uint8_t, 8> bytes{};
        autohog::store_u64_be(bytes.data(), value);
        transcript.insert(transcript.end(), bytes.begin(), bytes.end());
    };
    append_u64(params.n);
    append_u64(params.d);
    append_u64(params.k);
    append_u64(params.b);
    append_u64(params.b_r);
    append_u64(params.q);
    append_u64(params.t);
    append_u64(params.msbs_zero_padding_bit_count);

    const auto append_array = [&](const auto& bytes) {
        transcript.insert(transcript.end(), bytes.begin(), bytes.end());
    };
    append_array(params.sid);
    append_array(params.hash_ds);
    append_array(params.hash_t_ds);
    append_array(params.hash_agg_ds);
    append_array(params.hash_lmap_ds);
    append_array(params.hash_z_ds);
    append_array(params.hash_w_ds);
    append_array(params.hash_gamma_ds);

    std::vector<uint8_t> encoded_g1(kG1CompressedBytes);
    for (size_t i = 0; i < 2 * params.n; ++i) {
        if (i == params.n) {
            std::fill(encoded_g1.begin(), encoded_g1.end(), 0);
        } else {
            blst_p1_affine_compress(encoded_g1.data(), &g_list[i]);
        }
        transcript.insert(transcript.end(), encoded_g1.begin(), encoded_g1.end());
    }

    std::vector<uint8_t> encoded_g2(kG2CompressedBytes);
    for (size_t i = 0; i < params.n; ++i) {
        blst_p2_affine_compress(encoded_g2.data(), &g_hat_list[i]);
        transcript.insert(transcript.end(), encoded_g2.begin(), encoded_g2.end());
    }
    return transcript;
}

template<typename P>
static bool public_params_power_chain_usable(const PublicParams<P>& params) {
    const auto* g_list = static_cast<const blst_p1_affine*>(params.g_list_ptr);
    const auto* g_hat_list = static_cast<const blst_p2_affine*>(params.g_hat_list_ptr);
    const size_t n = static_cast<size_t>(params.n);

    // Anchor both groups to the same first exponent.
    if (!pairing_equal(&g_list[0], blst_p2_affine_generator(),
                       blst_p1_affine_generator(), &g_hat_list[0])) {
        return false;
    }

    const size_t g1_relation_count = 2 * n - 3;
    const size_t g2_relation_count = n - 1;
    const size_t total_relation_count = g1_relation_count + g2_relation_count;
    auto transcript = crs_validation_transcript(params);
    std::vector<uint8_t> challenge_bytes(total_relation_count * 16);
    shake256_hash(transcript.data(), transcript.size(),
                  challenge_bytes.data(), challenge_bytes.size());

    const auto challenge_scalar = [&](size_t index) {
        blst_scalar scalar{};
        std::memcpy(scalar.b, challenge_bytes.data() + index * 16, 16);
        if (scalar_is_zero(scalar)) scalar.b[0] = 1;
        return scalar;
    };

    std::vector<const blst_p1_affine*> g1_left_points;
    std::vector<const blst_p1_affine*> g1_right_points;
    std::vector<blst_scalar> g1_scalars;
    g1_left_points.reserve(g1_relation_count);
    g1_right_points.reserve(g1_relation_count);
    g1_scalars.reserve(g1_relation_count);
    size_t challenge_index = 0;
    for (size_t i = 1; i < 2 * n; ++i) {
        if (i == n || i == n + 1) continue;
        g1_left_points.push_back(&g_list[i]);
        g1_right_points.push_back(&g_list[i - 1]);
        g1_scalars.push_back(challenge_scalar(challenge_index++));
    }
    if (g1_scalars.size() != g1_relation_count) return false;

    const blst_p1 g1_left = p1_mult_pippenger(g1_left_points, g1_scalars);
    const blst_p1 g1_right = p1_mult_pippenger(g1_right_points, g1_scalars);
    blst_p1_affine g1_left_aff, g1_right_aff;
    blst_p1_to_affine(&g1_left_aff, &g1_left);
    blst_p1_to_affine(&g1_right_aff, &g1_right);
    if (!pairing_equal(&g1_left_aff, &g_hat_list[0],
                       &g1_right_aff, &g_hat_list[1])) {
        return false;
    }

    std::vector<const blst_p2_affine*> g2_left_points;
    std::vector<const blst_p2_affine*> g2_right_points;
    std::vector<blst_scalar> g2_scalars;
    g2_left_points.reserve(g2_relation_count);
    g2_right_points.reserve(g2_relation_count);
    g2_scalars.reserve(g2_relation_count);
    for (size_t i = 1; i < n; ++i) {
        g2_left_points.push_back(&g_hat_list[i]);
        g2_right_points.push_back(&g_hat_list[i - 1]);
        g2_scalars.push_back(challenge_scalar(challenge_index++));
    }
    if (g2_scalars.size() != g2_relation_count ||
        challenge_index != total_relation_count) {
        return false;
    }

    const blst_p2 g2_left = p2_mult_pippenger(g2_left_points, g2_scalars);
    const blst_p2 g2_right = p2_mult_pippenger(g2_right_points, g2_scalars);
    blst_p2_affine g2_left_aff, g2_right_aff;
    blst_p2_to_affine(&g2_left_aff, &g2_left);
    blst_p2_to_affine(&g2_right_aff, &g2_right);
    if (!pairing_equal(&g_list[0], &g2_left_aff,
                       &g_list[1], &g2_right_aff)) {
        return false;
    }

    // g_list[n] intentionally omits alpha^(n+1). Bridge the gap with G2[2].
    return pairing_equal(&g_list[n + 1], &g_hat_list[0],
                         &g_list[n - 1], &g_hat_list[2]);
}

template<typename P>
static bool public_params_usable(const PublicParams<P>& params) {
    if (!public_params_basic_usable(params)) {
        return false;
    }

    const auto* g_list = static_cast<const blst_p1_affine*>(params.g_list_ptr);
    const auto* g_hat_list = static_cast<const blst_p2_affine*>(params.g_hat_list_ptr);

    if (!blst_p1_affine_is_inf(&g_list[params.n])) {
        return false;
    }

    int ok = 1;
    #pragma omp parallel for schedule(static) reduction(&:ok)
    for (size_t i = 0; i < 2 * params.n; i++) {
        if (i == params.n) {
            continue;
        }
        ok &= g1_affine_usable(&g_list[i]) ? 1 : 0;
    }

    #pragma omp parallel for schedule(static) reduction(&:ok)
    for (size_t i = 0; i < params.n; i++) {
        ok &= g2_affine_usable(&g_hat_list[i]) ? 1 : 0;
    }

    return ok != 0 && public_params_power_chain_usable(params);
}

template<typename P>
bool validate_public_params(const PublicParams<P>& params) {
    return public_params_usable(params);
}


static std::vector<blst_fr> hash_to_fr_vec(
    const std::vector<uint8_t>& data, size_t count
) {
    std::vector<uint8_t> out(count * 32);
    shake256_hash(data.data(), data.size(), out.data(), out.size());
    std::vector<blst_fr> result(count);
    #pragma omp parallel for schedule(static) if (count >= 1024)
    for (size_t i = 0; i < count; i++) {
        blst_scalar s; blst_scalar_from_lendian(&s, out.data() + i * 32);
        blst_fr_from_scalar(&result[i], &s);
    }
    return result;
}


static std::vector<blst_fr> hash_to_fr_vec_128bit(
    const std::vector<uint8_t>& data, size_t count
) {
    std::vector<uint8_t> out(count * 16);
    shake256_hash(data.data(), data.size(), out.data(), out.size());
    std::vector<blst_fr> result(count);
    #pragma omp parallel for schedule(static) if (count >= 1024)
    for (size_t i = 0; i < count; i++) {
        blst_scalar s; memset(&s, 0, sizeof(s));
        memcpy(s.b, out.data() + i * 16, 16);
        blst_fr_from_scalar(&result[i], &s);
    }
    return result;
}

static std::string trim_ascii_whitespace(std::string value) {
    const auto not_space = [](unsigned char ch) {
        return !std::isspace(ch);
    };
    value.erase(value.begin(), std::find_if(value.begin(), value.end(), not_space));
    value.erase(std::find_if(value.rbegin(), value.rend(), not_space).base(), value.end());
    return value;
}

static bool is_lower_hex_string(const std::string& value, size_t expected_len) {
    if (value.size() != expected_len) return false;
    return std::all_of(value.begin(), value.end(), [](unsigned char ch) {
        return (ch >= '0' && ch <= '9') || (ch >= 'a' && ch <= 'f');
    });
}

static std::string normalize_sha256_hex(std::string value, const std::string& source_desc) {
    value = trim_ascii_whitespace(std::move(value));
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    if (!is_lower_hex_string(value, 64)) {
        throw std::runtime_error("Invalid SHA-256 text in " + source_desc);
    }
    return value;
}

static size_t next_pow2(size_t n) {
    if (n == 0) return 1;
    size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

struct NttPlan {
    size_t n = 0;
    size_t log2_n = 0;
    std::vector<size_t> bit_reverse;
    std::vector<std::vector<blst_fr>> forward_twiddles;
    std::vector<std::vector<blst_fr>> inverse_twiddles;
};

static blst_fr get_primitive_root_of_unity(size_t log2_n) {
    if (log2_n > autohog::zk::config::fr_two_adicity) {
        std::ostringstream oss;
        oss << "Current ZK curve " << autohog::zk::config::curve_name
            << " has Fr 2-adicity " << autohog::zk::config::fr_two_adicity
            << ", but this proof requires a 2^" << log2_n
            << "-th root of unity for polynomial multiplication";
        throw std::runtime_error(oss.str());
    }
#if defined(AUTOHOG_ZK_CURVE_bn_snark1)
  return get_primitive_root_of_unity_bn_snark1(log2_n);
#endif
static const uint8_t exp_bytes[32] = {
    0xff, 0xff, 0xff, 0xff, 0xfe, 0x5b, 0xfe, 0xff,
    0x02, 0xa4, 0xbd, 0x53, 0x05, 0xd8, 0xa1, 0x09,
    0x08, 0xd8, 0x39, 0x33, 0x48, 0x7d, 0x9d, 0x29,
    0x53, 0xa7, 0xed, 0x73, 0x00, 0x00, 0x00, 0x00
    };

    blst_fr result = fr_from_u64(1);
    blst_fr base = fr_from_u64(7);
    for (int i = 0; i < 256; i++) {
        int byte_i = i / 8, bit_i = i % 8;
        if ((exp_bytes[byte_i] >> bit_i) & 1) {
            blst_fr tmp; blst_fr_mul(&tmp, &result, &base);
            result = tmp;
        }
        blst_fr sq; blst_fr_mul(&sq, &base, &base);
        base = sq;
    }
    blst_fr root = result;
    for (size_t i = 0; i < 32 - log2_n; i++) {
        blst_fr sq; blst_fr_mul(&sq, &root, &root);
        root = sq;
    }
    return root;
}

static NttPlan build_ntt_plan(size_t n) {
    if (n == 0 || (n & (n - 1)) != 0) {
        throw std::runtime_error("NTT plan size must be a non-zero power of two");
    }

    NttPlan plan;
    plan.n = n;
    while ((size_t(1) << plan.log2_n) < n) {
        ++plan.log2_n;
    }

    if (plan.log2_n > autohog::zk::config::fr_two_adicity) {
        std::ostringstream oss;
        oss << "Current ZK curve " << autohog::zk::config::curve_name
            << " has Fr 2-adicity " << autohog::zk::config::fr_two_adicity
            << ", but this proof needs an NTT of size 2^" << plan.log2_n;
        throw std::runtime_error(oss.str());
    }

    plan.bit_reverse.resize(n);
    if (n > 1) {
        plan.bit_reverse[0] = 0;
        for (size_t i = 1; i < n; ++i) {
            plan.bit_reverse[i] =
                (plan.bit_reverse[i >> 1] >> 1) |
                ((i & 1) << (plan.log2_n - 1));
        }
    }

    plan.forward_twiddles.reserve(plan.log2_n);
    plan.inverse_twiddles.reserve(plan.log2_n);
    for (size_t stage = 0; stage < plan.log2_n; ++stage) {
        const size_t len = size_t(1) << (stage + 1);
        const size_t half = len >> 1;
        const blst_fr root = get_primitive_root_of_unity(stage + 1);
        const blst_fr inv_root = fr_inv(root);

        std::vector<blst_fr> forward(half);
        std::vector<blst_fr> inverse(half);
        forward[0] = fr_from_u64(1);
        inverse[0] = fr_from_u64(1);
        for (size_t j = 1; j < half; ++j) {
            forward[j] = fr_mul(forward[j - 1], root);
            inverse[j] = fr_mul(inverse[j - 1], inv_root);
        }
        plan.forward_twiddles.push_back(std::move(forward));
        plan.inverse_twiddles.push_back(std::move(inverse));
    }

    return plan;
}

static const NttPlan& get_ntt_plan(size_t n) {
    static std::mutex cache_mutex;
    static std::unordered_map<size_t, std::shared_ptr<NttPlan>> cache;

    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        auto it = cache.find(n);
        if (it != cache.end()) {
            return *it->second;
        }
    }

    auto plan = std::make_shared<NttPlan>(build_ntt_plan(n));

    std::lock_guard<std::mutex> lock(cache_mutex);
    auto [it, inserted] = cache.emplace(n, std::move(plan));
    return *it->second;
}

static void bit_reverse_permute(std::vector<blst_fr>& a, const NttPlan& plan) {
    for (size_t i = 0; i < plan.n; ++i) {
        const size_t j = plan.bit_reverse[i];
        if (i < j) std::swap(a[i], a[j]);
    }
}

static void ntt(std::vector<blst_fr>& a, bool inverse) {
    const size_t n = a.size();
    if (n == 1) return;

    const NttPlan& plan = get_ntt_plan(n);
    bit_reverse_permute(a, plan);

    for (size_t stage = 0, len = 2; len <= n; ++stage, len <<= 1) {
        size_t half = len >> 1;
        const auto& twiddles = inverse
            ? plan.inverse_twiddles[stage]
            : plan.forward_twiddles[stage];

        #pragma omp parallel for schedule(static) if (n >= 4096)
        for (ptrdiff_t block = 0; block < static_cast<ptrdiff_t>(n); block += static_cast<ptrdiff_t>(len)) {
            const size_t i = static_cast<size_t>(block);
            for (size_t j = 0; j < half; j++) {
                const blst_fr& w = twiddles[j];
                blst_fr u = a[i + j];
                blst_fr t; blst_fr_mul(&t, &a[i + j + half], &w);
                blst_fr_add(&a[i + j], &u, &t);
                blst_fr_sub(&a[i + j + half], &u, &t);
            }
        }
    }

    if (inverse) {
        blst_fr n_fr = fr_from_u64(static_cast<uint64_t>(n));
        blst_fr n_inv = fr_inv(n_fr);
        #pragma omp parallel for schedule(static) if (n >= 4096)
        for (ptrdiff_t idx = 0; idx < static_cast<ptrdiff_t>(n); ++idx) {
            const size_t i = static_cast<size_t>(idx);
            blst_fr tmp; blst_fr_mul(&tmp, &a[i], &n_inv);
            a[i] = tmp;
        }
    }
}


static std::vector<blst_fr> poly_mul_fr(
    const std::vector<blst_fr>& a,
    const std::vector<blst_fr>& b
) {
    if (a.empty() || b.empty()) return {};

    size_t result_len = a.size() + b.size() - 1;
    size_t ntt_size   = next_pow2(result_len);

    if (ntt_size > (1u << 30)) {
        throw std::runtime_error("poly_mul_fr: polynomial too large for NTT (>2^30)");
    }

    std::vector<blst_fr> fa(ntt_size, fr_zero());
    std::vector<blst_fr> fb(ntt_size, fr_zero());
    for (size_t i = 0; i < a.size(); i++) fa[i] = a[i];
    for (size_t i = 0; i < b.size(); i++) fb[i] = b[i];

    ntt(fa, false);
    ntt(fb, false);

    #pragma omp parallel for schedule(static) if (ntt_size >= 4096)
    for (size_t i = 0; i < ntt_size; i++) {
        blst_fr tmp; blst_fr_mul(&tmp, &fa[i], &fb[i]);
        fa[i] = tmp;
    }

    ntt(fa, true);

    fa.resize(result_len);
    return fa;
}

static std::vector<blst_fr> poly_sub_fr(
    const std::vector<blst_fr>& a, const std::vector<blst_fr>& b
) {
    std::vector<blst_fr> result(std::max(a.size(), b.size()), fr_zero());
    #pragma omp parallel for schedule(static) if (result.size() >= 4096)
    for (size_t i = 0; i < result.size(); i++) {
        const blst_fr ai = (i < a.size()) ? a[i] : fr_zero();
        const blst_fr bi = (i < b.size()) ? b[i] : fr_zero();
        result[i] = fr_sub(ai, bi);
    }
    return result;
}

static size_t ilog2_u64(uint64_t v) {
    if (v <= 1) return 0;
    size_t r = 0; uint64_t t = v;
    while (t > 1) { t >>= 1; r++; }
    return r;
}

static int64_t exact_divide_by_q(
    __int128_t value,
    __uint128_t decoded_q,
    const char* relation_name,
    size_t index
) {
    const __int128_t q = static_cast<__int128_t>(decoded_q);
    if (value % q != 0) {
        std::ostringstream oss;
        oss << relation_name << " residual is non-zero at index " << index;
        throw std::runtime_error(oss.str());
    }

    const __int128_t quotient = value / q;
    if (quotient < static_cast<__int128_t>(std::numeric_limits<int64_t>::min()) ||
        quotient > static_cast<__int128_t>(std::numeric_limits<int64_t>::max())) {
        std::ostringstream oss;
        oss << relation_name << " quotient is out of int64 range at index " << index;
        throw std::runtime_error(oss.str());
    }

    return static_cast<int64_t>(quotient);
}

static std::vector<int64_t> compute_r1(
    const std::vector<int64_t>& e1,
    const std::vector<int64_t>& c1,
    const std::vector<int64_t>& a,
    const std::vector<int64_t>& r,
    size_t d,
    __uint128_t decoded_q
) {
    std::vector<__int128_t> acc(d, 0);
    for (size_t i = 0; i < d; i++) {
        acc[i] = static_cast<__int128_t>(e1[i]) - static_cast<__int128_t>(c1[i]);
    }

    const int num_threads = omp_in_parallel() ? 1 : omp_get_max_threads();
    std::vector<std::vector<__int128_t>> partials(
        num_threads, std::vector<__int128_t>(d, 0)
    );

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        auto& local = partials[tid];

        #pragma omp for schedule(static)
        for (ptrdiff_t ii = 0; ii < static_cast<ptrdiff_t>(d); ii++) {
            const size_t i = static_cast<size_t>(ii);
            for (size_t j = 0; j < d; j++) {
                size_t idx = i + j;
                const __int128_t term =
                    static_cast<__int128_t>(a[i]) * static_cast<__int128_t>(r[d - j - 1]);
                if (idx < d) {
                    local[idx] += term;
                } else {
                    local[idx - d] -= term;
                }
            }
        }
    }

    for (const auto& local : partials) {
        for (size_t i = 0; i < d; i++) {
            acc[i] += local[i];
        }
    }

    std::vector<int64_t> r1(d, 0);
    #pragma omp parallel for schedule(static) if (d >= 512)
    for (size_t i = 0; i < d; i++) {
        __int128_t val = acc[i];
        r1[i] = exact_divide_by_q(val, decoded_q, "RLWE c1 relation", i);
    }
    return r1;
}

static std::vector<int64_t> compute_r2(
    const std::vector<int64_t>& e2,
    const std::vector<int64_t>& c2,
    const std::vector<int64_t>& m,
    const std::vector<int64_t>& b,
    const std::vector<int64_t>& r,
    size_t d,
    size_t row_offset,
    uint64_t delta,
    __uint128_t decoded_q
) {
    std::vector<int64_t> r2(m.size(), 0);
    #pragma omp parallel for schedule(static) if (m.size() >= 64)
    for (size_t i = 0; i < m.size(); i++) {
        const size_t global_i = row_offset + i;
        __int128_t val = (__int128_t)delta * m[i] + e2[i] - c2[i];
        __int128_t dot = 0;
        for (size_t j = 0; j < d; j++) {
            const __int128_t coeff_b =
                (global_i + j < d)
                    ? static_cast<__int128_t>(b[d - j - global_i - 1])
                    : -static_cast<__int128_t>(b[2 * d - j - global_i - 1]);
            dot += static_cast<__int128_t>(r[d - j - 1]) * coeff_b;
        }
        val += dot;
        r2[i] = exact_divide_by_q(val, decoded_q, "RLWE c2 relation", i);
    }
    return r2;
}

template<typename P>
static void check_rlwe_plaintext_relation(
    const PublicParams<P>& params,
    const PublicCommit<P>& public_commit,
    const PrivateCommit<P>& private_commit,
    uint64_t delta,
    __uint128_t decoded_q,
    std::vector<int64_t>& r1,
    std::vector<int64_t>& r2
) {
    for (size_t i = 0; i < private_commit.m.size(); i++) {
        if (private_commit.m[i] != 0 && private_commit.m[i] != 1) {
            std::ostringstream oss;
            oss << "Plaintext witness must be binary at index " << i;
            throw std::runtime_error(oss.str());
        }
    }

    const size_t row_offset = params.d - public_commit.c2.size();
    r1 = compute_r1(
        private_commit.e1,
        public_commit.c1,
        public_commit.a,
        private_commit.r,
        params.d,
        decoded_q
    );
    r2 = compute_r2(
        private_commit.e2,
        public_commit.c2,
        private_commit.m,
        public_commit.b,
        private_commit.r,
        params.d,
        row_offset,
        delta,
        decoded_q
    );
}

template<typename P>
static std::vector<uint8_t> build_statement_x_bytes(
    const PublicParams<P>& params,
    const PublicCommit<P>& public_commit
) {
    std::vector<uint8_t> x_bytes;
    auto app = [&](uint64_t v) {
        auto b = u64_to_le_bytes(v);
        x_bytes.insert(x_bytes.end(), b.begin(), b.end());
    };

    app(params.q);
    app(static_cast<uint64_t>(params.d));
    app(params.b);
    app(params.t);
    app(params.msbs_zero_padding_bit_count);
    app(params.b_r);
    app(static_cast<uint64_t>(public_commit.c2.size()));

    auto ai  = i64_vec_to_bytes(public_commit.a);
    auto bi  = i64_vec_to_bytes(public_commit.b);
    auto c1i = i64_vec_to_bytes(public_commit.c1);
    auto c2i = i64_vec_to_bytes(public_commit.c2);
    x_bytes.insert(x_bytes.end(), ai.begin(),  ai.end());
    x_bytes.insert(x_bytes.end(), bi.begin(),  bi.end());
    x_bytes.insert(x_bytes.end(), c1i.begin(), c1i.end());
    x_bytes.insert(x_bytes.end(), c2i.begin(), c2i.end());
    return x_bytes;
}

static void compute_a_theta_fn(
    const std::vector<blst_fr>& theta0,
    const std::vector<blst_fr>& m_bit_theta,
    size_t d, size_t k,
    const std::vector<int64_t>& a_int,
    const std::vector<int64_t>& b_int,
    std::vector<blst_fr>& a_theta,
    size_t row_offset,
    uint64_t effective_t,
    uint64_t delta,
    uint64_t b_i,
    uint64_t b_r,
    __uint128_t decoded_q
) {
    blst_fr q_fr     = fr_from_u128(decoded_q);
    blst_fr delta_fr = fr_from_u64(delta);

    std::vector<blst_fr> theta1(theta0.begin(), theta0.begin() + d);
    std::vector<blst_fr> theta2(theta0.begin() + d, theta0.begin() + d + k);

    size_t offset = 0;

    #pragma omp parallel for schedule(static) if (d >= 256)
    for (size_t i = 0; i < d; i++) {
        blst_fr dot = fr_zero();
        for (size_t j = 0; j < d; j++) {
            blst_fr coeff = (i <= j)
                ? fr_from_i64(a_int[j - i])
                : fr_neg(fr_from_i64(a_int[d + j - i]));
            dot = fr_add(dot, fr_mul(coeff, theta1[j]));
        }
        for (size_t j = 0; j < k; j++) {
            const size_t global_j = row_offset + j;
            blst_fr coeff = (i + global_j < d)
                ? fr_from_i64(b_int[d - i - global_j - 1])
                : fr_neg(fr_from_i64(b_int[2 * d - i - global_j - 1]));
            dot = fr_add(dot, fr_mul(coeff, theta2[j]));
        }
        a_theta[offset + i] = dot;
    }
    offset += d;

    size_t step_t = ilog2_u64(effective_t);
    #pragma omp parallel for schedule(static) if (k >= 64)
    for (size_t i = 0; i < k; i++) {
        blst_fr cur = fr_mul(delta_fr, theta2[i]);
        for (size_t j = 0; j < step_t; j++) {
            if (j > 0) cur = fr_add(cur, cur);
            blst_fr value = cur;
            if (j > 0) {
                value = fr_add(value, m_bit_theta[i * (step_t - 1) + (j - 1)]);
            }
            a_theta[offset + step_t * i + j] = value;
        }
    }
    offset += k * step_t;

    size_t step_b = 1 + ilog2_u64(b_i);
    #pragma omp parallel for schedule(static) if (d >= 256)
    for (size_t i = 0; i < d; i++) {
        blst_fr cur = theta1[i];
        for (size_t j = 0; j < step_b; j++) {
            if (j > 0) cur = fr_add(cur, cur);
            a_theta[offset + step_b * i + j] = (j == step_b - 1) ? fr_neg(cur) : cur;
        }
    }
    offset += d * step_b;

    #pragma omp parallel for schedule(static) if (k >= 64)
    for (size_t i = 0; i < k; i++) {
        blst_fr cur = theta2[i];
        for (size_t j = 0; j < step_b; j++) {
            if (j > 0) cur = fr_add(cur, cur);
            a_theta[offset + step_b * i + j] = (j == step_b - 1) ? fr_neg(cur) : cur;
        }
    }
    offset += k * step_b;

    size_t step_br = 1 + ilog2_u64(b_r);
    #pragma omp parallel for schedule(static) if (d >= 256)
    for (size_t i = 0; i < d; i++) {
        blst_fr cur = fr_mul(fr_neg(q_fr), theta1[i]);
        for (size_t j = 0; j < step_br; j++) {
            if (j > 0) cur = fr_add(cur, cur);
            a_theta[offset + step_br * i + j] = (j == step_br - 1) ? fr_neg(cur) : cur;
        }
    }
    offset += d * step_br;

    #pragma omp parallel for schedule(static) if (k >= 64)
    for (size_t i = 0; i < k; i++) {
        blst_fr cur = fr_mul(fr_neg(q_fr), theta2[i]);
        for (size_t j = 0; j < step_br; j++) {
            if (j > 0) cur = fr_add(cur, cur);
            a_theta[offset + step_br * i + j] = (j == step_br - 1) ? fr_neg(cur) : cur;
        }
    }
}


Proof::Proof() {
    c_hat_ptr = new blst_p2_affine;
    c_y_ptr = new blst_p1_affine;
    pi_ptr = new blst_p1_affine;
    blst_scalar zero_scalar;
    std::memset(&zero_scalar, 0, sizeof(zero_scalar));
    blst_p1 zero_g1;
    blst_p1_from_affine(&zero_g1, blst_p1_affine_generator());
    blst_p1_mult(&zero_g1, &zero_g1, zero_scalar.b, 256);
    blst_p2 zero_g2;
    blst_p2_from_affine(&zero_g2, blst_p2_affine_generator());
    blst_p2_mult(&zero_g2, &zero_g2, zero_scalar.b, 256);
    blst_p2_to_affine(static_cast<blst_p2_affine*>(c_hat_ptr), &zero_g2);
    blst_p1_to_affine(static_cast<blst_p1_affine*>(c_y_ptr), &zero_g1);
    blst_p1_to_affine(static_cast<blst_p1_affine*>(pi_ptr), &zero_g1);
}

Proof::~Proof() {
    delete static_cast<blst_p2_affine*>(c_hat_ptr);
    delete static_cast<blst_p1_affine*>(c_y_ptr);
    delete static_cast<blst_p1_affine*>(pi_ptr);
}

Proof::Proof(Proof&& other) noexcept {
    c_hat_ptr = other.c_hat_ptr;
    c_y_ptr = other.c_y_ptr;
    pi_ptr = other.pi_ptr;
    other.c_hat_ptr = nullptr;
    other.c_y_ptr = nullptr;
    other.pi_ptr = nullptr;
}

Proof& Proof::operator=(Proof&& other) noexcept {
    if (this != &other) {
        delete static_cast<blst_p2_affine*>(c_hat_ptr);
        delete static_cast<blst_p1_affine*>(c_y_ptr);
        delete static_cast<blst_p1_affine*>(pi_ptr);
        c_hat_ptr = other.c_hat_ptr;
        c_y_ptr = other.c_y_ptr;
        pi_ptr = other.pi_ptr;
        other.c_hat_ptr = nullptr;
        other.c_y_ptr = nullptr;
        other.pi_ptr = nullptr;
    }
    return *this;
}

std::vector<uint8_t> Proof::serialize() const {
    std::vector<uint8_t> buffer;
    buffer.reserve(kProofSerializedBytes);

    std::vector<uint8_t> c_hat_bytes(kG2CompressedBytes);
    blst_p2_affine_compress(c_hat_bytes.data(), static_cast<const blst_p2_affine*>(c_hat_ptr));
    buffer.insert(buffer.end(), c_hat_bytes.begin(), c_hat_bytes.end());

    std::vector<uint8_t> c_y_bytes(kG1CompressedBytes);
    blst_p1_affine_compress(c_y_bytes.data(), static_cast<const blst_p1_affine*>(c_y_ptr));
    buffer.insert(buffer.end(), c_y_bytes.begin(), c_y_bytes.end());

    std::vector<uint8_t> pi_bytes(kG1CompressedBytes);
    blst_p1_affine_compress(pi_bytes.data(), static_cast<const blst_p1_affine*>(pi_ptr));
    buffer.insert(buffer.end(), pi_bytes.begin(), pi_bytes.end());

    return buffer;
}

Proof Proof::deserialize(const std::vector<uint8_t>& data) {
    if (data.size() != kProofSerializedBytes) {
        throw std::runtime_error("Invalid proof size");
    }
    Proof proof;
    const size_t c_hat_offset = 0;
    const size_t c_y_offset = kG2CompressedBytes;
    const size_t pi_offset = kG2CompressedBytes + kG1CompressedBytes;
    if (blst_p2_uncompress(static_cast<blst_p2_affine*>(proof.c_hat_ptr), data.data() + c_hat_offset) != BLST_SUCCESS ||
        blst_p1_uncompress(static_cast<blst_p1_affine*>(proof.c_y_ptr), data.data() + c_y_offset) != BLST_SUCCESS ||
        blst_p1_uncompress(static_cast<blst_p1_affine*>(proof.pi_ptr), data.data() + pi_offset) != BLST_SUCCESS) {
        throw std::runtime_error("Invalid compressed proof point");
    }
    return proof;
}

template<typename P>
PublicParams<P>::PublicParams() : g_list_ptr(nullptr), g_hat_list_ptr(nullptr), 
    n(0), d(0), k(0), b(0), b_r(0), q(0), t(0),
    msbs_zero_padding_bit_count(0), integrity_checked(false) {}

template<typename P>
PublicParams<P>::~PublicParams() {
    if (g_list_ptr) delete[] static_cast<blst_p1_affine*>(g_list_ptr);
    if (g_hat_list_ptr) delete[] static_cast<blst_p2_affine*>(g_hat_list_ptr);
}

template<typename P>
PublicParams<P>::PublicParams(PublicParams&& other) noexcept {
    g_list_ptr = other.g_list_ptr;
    g_hat_list_ptr = other.g_hat_list_ptr;
    n = other.n; d = other.d; k = other.k; 
    b = other.b; b_r = other.b_r; q = other.q; t = other.t;
    msbs_zero_padding_bit_count = other.msbs_zero_padding_bit_count;
    integrity_checked = other.integrity_checked;
    sid = other.sid;
    hash_ds = other.hash_ds; hash_t_ds = other.hash_t_ds;
    hash_agg_ds = other.hash_agg_ds; hash_lmap_ds = other.hash_lmap_ds;
    hash_z_ds = other.hash_z_ds; hash_w_ds = other.hash_w_ds;
    hash_gamma_ds = other.hash_gamma_ds;
    other.g_list_ptr = nullptr;
    other.g_hat_list_ptr = nullptr;
    other.integrity_checked = false;
}

template<typename P>
PublicParams<P>& PublicParams<P>::operator=(PublicParams&& other) noexcept {
    if (this != &other) {
        if (g_list_ptr) delete[] static_cast<blst_p1_affine*>(g_list_ptr);
        if (g_hat_list_ptr) delete[] static_cast<blst_p2_affine*>(g_hat_list_ptr);
        g_list_ptr = other.g_list_ptr;
        g_hat_list_ptr = other.g_hat_list_ptr;
        n = other.n; d = other.d; k = other.k;
        b = other.b; b_r = other.b_r; q = other.q; t = other.t;
        msbs_zero_padding_bit_count = other.msbs_zero_padding_bit_count;
        integrity_checked = other.integrity_checked;
        sid = other.sid;
        hash_ds = other.hash_ds; hash_t_ds = other.hash_t_ds;
        hash_agg_ds = other.hash_agg_ds; hash_lmap_ds = other.hash_lmap_ds;
        hash_z_ds = other.hash_z_ds; hash_w_ds = other.hash_w_ds;
        hash_gamma_ds = other.hash_gamma_ds;
        other.g_list_ptr = nullptr;
        other.g_hat_list_ptr = nullptr;
        other.integrity_checked = false;
    }
    return *this;
}

template class PublicParams<TFHEpp::lvl2param>;
template struct PublicCommit<TFHEpp::lvl2param>;
template struct PrivateCommit<TFHEpp::lvl2param>;

template<typename P>
PublicParams<P> generate_public_params(
    uint64_t k,
    uint64_t b, uint64_t b_r, uint64_t q, uint64_t t,
    uint64_t msbs_zero_padding_bit_count,
    const std::array<uint8_t, SEED_LEN>& alpha_seed
) {
    PublicParams<P> params;
    params.d = P::n;
    const size_t effective_k = (k == 0) ? default_statement_k_lvl2() : static_cast<size_t>(k);
    if (effective_k == 0 || effective_k > P::n) {
        throw std::runtime_error("CRS message dimension k must be in [1, P::n]");
    }
    params.k = effective_k;
    params.b = b;
    params.b_r = b_r;
    params.q = q;
    params.t = t;
    params.msbs_zero_padding_bit_count = msbs_zero_padding_bit_count;
    
    const size_t d = params.d;
    const size_t k_max = params.k;

    uint64_t effective_t = t >> msbs_zero_padding_bit_count;
    size_t step_t  = ilog2_u64(effective_t);
    size_t step_b  = 1 + ilog2_u64(b);
    size_t step_br = 1 + ilog2_u64(b_r);

    size_t big_d = d + k_max * step_t + (d + k_max) * (step_b + step_br);
    params.n = big_d + 1;
    
    // Derive every public domain separator from all 256 seed bits. The old
    // construction copied only the first 64 bits into mt19937_64, reducing the
    // effective input space even though alpha_seed itself is 32 bytes.
    static constexpr std::array<uint8_t, 30> ds_label{
        'A', 'U', 'T', 'O', 'H', 'O', 'G', '-', 'Z', 'K', '-', 'D', 'O', 'M', 'A', 'I', 'N',
        '-', 'S', 'E', 'P', 'A', 'R', 'A', 'T', 'O', 'R', '-', 'V', '1'
    };
    std::vector<uint8_t> ds_input;
    ds_input.reserve(ds_label.size() + alpha_seed.size());
    ds_input.insert(ds_input.end(), ds_label.begin(), ds_label.end());
    ds_input.insert(ds_input.end(), alpha_seed.begin(), alpha_seed.end());

    const size_t ds_output_size =
        params.hash_ds.size() + params.hash_t_ds.size() + params.hash_agg_ds.size() +
        params.hash_lmap_ds.size() + params.hash_z_ds.size() + params.hash_w_ds.size() +
        params.hash_gamma_ds.size() + params.sid.size();
    std::vector<uint8_t> ds_output(ds_output_size);
    shake256_hash(ds_input.data(), ds_input.size(), ds_output.data(), ds_output.size());
    autohog::secure_erase(ds_input.data(), ds_input.size());

    size_t ds_offset = 0;
    const auto copy_derived_bytes = [&](auto& destination) {
        std::copy_n(ds_output.data() + ds_offset, destination.size(), destination.data());
        ds_offset += destination.size();
    };
    copy_derived_bytes(params.hash_ds);
    copy_derived_bytes(params.hash_t_ds);
    copy_derived_bytes(params.hash_agg_ds);
    copy_derived_bytes(params.hash_lmap_ds);
    copy_derived_bytes(params.hash_z_ds);
    copy_derived_bytes(params.hash_w_ds);
    copy_derived_bytes(params.hash_gamma_ds);
    copy_derived_bytes(params.sid);
    
    ScrubbedValue<blst_scalar> alpha;
    blst_scalar_from_lendian(&alpha.get(), alpha_seed.data());
    ScrubbedValue<blst_fr> alpha_fr;
    blst_fr_from_scalar(&alpha_fr.get(), &alpha.get());
    if (fr_is_zero(alpha_fr.get())) {
        throw std::runtime_error("Invalid CRS seed: alpha maps to zero");
    }
    
    params.g_list_ptr    = new blst_p1_affine[2 * params.n];
    params.g_hat_list_ptr = new blst_p2_affine[params.n];
    
    blst_p1_affine* g_list = static_cast<blst_p1_affine*>(params.g_list_ptr);
    blst_p2_affine* g_hat_list = static_cast<blst_p2_affine*>(params.g_hat_list_ptr);
    
    const blst_p1_affine* g1_gen = blst_p1_affine_generator();
    const blst_p2_affine* g2_gen = blst_p2_affine_generator();
    
    blst_scalar zero_scalar;
    std::memset(&zero_scalar, 0, sizeof(zero_scalar));
    blst_p1 zero_g1;
    blst_p1_from_affine(&zero_g1, g1_gen);
    blst_p1_mult(&zero_g1, &zero_g1, zero_scalar.b, 256);

    blst_p1_affine zero_g1_aff;
    blst_p1_to_affine(&zero_g1_aff, &zero_g1);

    const size_t g1_count = 2 * params.n;
    const size_t g1_chunks = power_chain_parallel_chunks(g1_count);
    #pragma omp parallel for schedule(static) if (g1_chunks > 1)
    for (int chunk = 0; chunk < static_cast<int>(g1_chunks); ++chunk) {
        const size_t begin = g1_count * static_cast<size_t>(chunk) / g1_chunks;
        const size_t end = g1_count * static_cast<size_t>(chunk + 1) / g1_chunks;
        if (begin >= end) {
            continue;
        }

        ScrubbedValue<blst_fr> start_power;
        start_power.get() = fr_pow_u64(alpha_fr.get(), static_cast<uint64_t>(begin + 1));
        ScrubbedValue<blst_scalar> start_scalar;
        start_scalar.get() = fr_to_scalar(start_power.get());

        blst_p1 current_g1;
        blst_p1_from_affine(&current_g1, g1_gen);
        blst_p1_mult(&current_g1, &current_g1, start_scalar.get().b, 256);

        for (size_t i = begin; i < end; ++i) {
            if (i == params.n) {
                g_list[i] = zero_g1_aff;
            } else {
                blst_p1_to_affine(&g_list[i], &current_g1);
            }
            if (i + 1 < end) {
                blst_p1 next_g1;
                blst_p1_mult(&next_g1, &current_g1, alpha.get().b, 256);
                current_g1 = next_g1;
            }
        }
    }

    const size_t g2_count = params.n;
    const size_t g2_chunks = power_chain_parallel_chunks(g2_count);
    #pragma omp parallel for schedule(static) if (g2_chunks > 1)
    for (int chunk = 0; chunk < static_cast<int>(g2_chunks); ++chunk) {
        const size_t begin = g2_count * static_cast<size_t>(chunk) / g2_chunks;
        const size_t end = g2_count * static_cast<size_t>(chunk + 1) / g2_chunks;
        if (begin >= end) {
            continue;
        }

        ScrubbedValue<blst_fr> start_power;
        start_power.get() = fr_pow_u64(alpha_fr.get(), static_cast<uint64_t>(begin + 1));
        ScrubbedValue<blst_scalar> start_scalar;
        start_scalar.get() = fr_to_scalar(start_power.get());

        blst_p2 current_g2;
        blst_p2_from_affine(&current_g2, g2_gen);
        blst_p2_mult(&current_g2, &current_g2, start_scalar.get().b, 256);

        for (size_t i = begin; i < end; ++i) {
            blst_p2_to_affine(&g_hat_list[i], &current_g2);
            if (i + 1 < end) {
                blst_p2 next_g2;
                blst_p2_mult(&next_g2, &current_g2, alpha.get().b, 256);
                current_g2 = next_g2;
            }
        }
    }
    if (!public_params_usable(params)) {
        throw std::runtime_error("Generated invalid CRS");
    }
    params.integrity_checked = true;
    return params;
}

std::string compute_file_sha256_hex(const std::string& filename) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) throw std::runtime_error("Cannot open file for SHA-256: " + filename);

    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    if (!ctx) throw std::runtime_error("Failed to create SHA-256 EVP context");
    if (EVP_DigestInit_ex(ctx, EVP_sha256(), nullptr) != 1) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("Failed to initialize SHA-256 digest");
    }

    std::array<char, 1 << 16> buffer{};
    while (ifs) {
        ifs.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const std::streamsize read_count = ifs.gcount();
        if (read_count > 0 &&
            EVP_DigestUpdate(ctx, buffer.data(), static_cast<size_t>(read_count)) != 1) {
            EVP_MD_CTX_free(ctx);
            throw std::runtime_error("Failed to update SHA-256 digest");
        }
    }

    if (!ifs.eof()) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("Failed while reading file for SHA-256: " + filename);
    }

    unsigned char md[EVP_MAX_MD_SIZE];
    unsigned int md_len = 0;
    if (EVP_DigestFinal_ex(ctx, md, &md_len) != 1) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("Failed to finalize SHA-256 digest");
    }
    EVP_MD_CTX_free(ctx);

    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (unsigned int i = 0; i < md_len; ++i) {
        oss << std::setw(2) << static_cast<unsigned int>(md[i]);
    }
    return oss.str();
}

static std::string compute_memory_sha256_hex(const void* data, size_t len) {
    if (data == nullptr && len != 0) {
        throw std::invalid_argument("Cannot hash a null CRS buffer");
    }

    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    if (!ctx) throw std::runtime_error("Failed to create SHA-256 EVP context");
    if (EVP_DigestInit_ex(ctx, EVP_sha256(), nullptr) != 1 ||
        (len != 0 && EVP_DigestUpdate(ctx, data, len) != 1)) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("Failed to hash CRS memory buffer");
    }

    unsigned char md[EVP_MAX_MD_SIZE];
    unsigned int md_len = 0;
    if (EVP_DigestFinal_ex(ctx, md, &md_len) != 1) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("Failed to finalize CRS SHA-256 digest");
    }
    EVP_MD_CTX_free(ctx);

    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (unsigned int i = 0; i < md_len; ++i) {
        oss << std::setw(2) << static_cast<unsigned int>(md[i]);
    }
    return oss.str();
}

void write_file_sha256_hex(const std::string& filename, const std::string& hash_filename) {
    const std::string digest = compute_file_sha256_hex(filename);
    std::ofstream ofs(hash_filename, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Cannot open hash file for writing: " + hash_filename);
    }
    ofs << digest << '\n';
    if (!ofs) {
        throw std::runtime_error("Failed writing hash file: " + hash_filename);
    }
}

template<typename P>
PublicCommit<P> make_public_commit(
    const PKE::SeededRlwePkPair<P>& spk,
    const TFHEpp::TRLWE<P>& ct,
    uint64_t num_messages
) {
    if (num_messages == 0 || num_messages > P::n) {
        throw std::runtime_error("PublicCommit message dimension must be in [1, P::n]");
    }

    PublicCommit<P> commit;
    commit.d = P::n;
    commit.k = num_messages;
    
    commit.a.resize(P::n);
    TFHEpp::Polynomial<P> a_poly;
    PKE::uniform_poly_from_seed<P>(a_poly, spk.seed, 0);
    for (size_t i = 0; i < P::n; i++) {
        commit.a[i] = static_cast<int64_t>(a_poly[i]);
    }
    
    commit.b.resize(P::n);
    for (size_t i = 0; i < P::n; i++) {
        commit.b[i] = static_cast<int64_t>(spk.b[i]);
    }
    
    commit.c1.resize(P::n);
    for (size_t i = 0; i < P::n; i++) {
        commit.c1[i] = static_cast<int64_t>(ct[0][i]);
    }
    
    commit.c2.resize(num_messages);
    for (size_t i = 0; i < num_messages; i++) {
        commit.c2[i] = static_cast<int64_t>(ct[P::k][num_messages - 1 - i] + static_cast<typename P::T>(P::μ));
    }
    
    return commit;
}

template<typename P>
PrivateCommit<P> make_private_commit(
    const std::vector<int64_t>& r,
    const std::vector<int64_t>& e1,
    const std::vector<int64_t>& e2,
    const std::vector<int64_t>& m
) {
    if (r.size() != P::n || e1.size() != P::n) {
        throw std::runtime_error("PrivateCommit r/e1 vectors must have length P::n");
    }
    if (m.empty() || m.size() > P::n || e2.size() != m.size()) {
        throw std::runtime_error("PrivateCommit m/e2 vectors must have matching length in [1, P::n]");
    }
    
    PrivateCommit<P> commit;
    commit.d = P::n;
    commit.k = m.size();
    commit.r = r;
    commit.e1 = e1;
    commit.e2 = e2;
    commit.m = m;
    return commit;
}

template<typename P>
Proof prove(
    const PublicParams<P>& params,
    const PublicCommit<P>& public_commit,
    const PrivateCommit<P>& private_commit,
    const std::vector<uint8_t>& metadata,
    ComputeLoad load,
    const std::array<uint8_t, SEED_LEN>& seed
) {
    (void)load; 
    const bool profile_enabled = zk_profile_enabled();
    const auto prove_start_total = profile_enabled ? ProfileClock::now() : ProfileClock::time_point{};
    auto prove_stage_start = prove_start_total;
    std::vector<std::pair<const char*, double>> prove_profile_parts;

    if (!public_params_basic_usable(params)) {
        throw std::runtime_error("Invalid CRS structure for pke-v1 proof");
    }
    if (!params.integrity_checked && !public_params_usable(params)) {
        throw std::runtime_error("Unusable CRS for pke-v1 proof");
    }

    const size_t d = params.d;
    const size_t k = public_commit.c2.size();

    if (k == 0 || k > params.k ||
        public_commit.a.size() != d || public_commit.b.size() != d ||
        public_commit.c1.size() != d || public_commit.c2.size() != k ||
        private_commit.r.size() != d || private_commit.e1.size() != d ||
        private_commit.e2.size() != k || private_commit.m.size() != k) {
        throw std::runtime_error("Invalid commit dimensions for pke-v1 proof");
    }

    __uint128_t decoded_q  = decode_q(params.q);
    uint64_t    delta_enc  = static_cast<uint64_t>(decoded_q / params.t);
    uint64_t    effective_t = params.t >> params.msbs_zero_padding_bit_count;

    size_t step_t  = ilog2_u64(effective_t);
    size_t step_b  = 1 + ilog2_u64(params.b);
    size_t step_br = 1 + ilog2_u64(params.b_r);
    size_t m_bit_constraints = k * (step_t - 1);
    size_t big_d   = d + k * step_t + (d + k) * (step_b + step_br);
    size_t n       = params.n;

    if (big_d + 1 > n) {
        throw std::runtime_error("CRS dimensions are incompatible with proof statement");
    }

    const size_t row_offset = d - k;
    std::vector<int64_t> r1;
    std::vector<int64_t> r2;
    // This proof variant is specifically for boolean plaintext RLWE statements:
    // c1 = a*r + e1 mod q and c2 = b*r + e2 + Delta*m mod q with m in {0, 1}.
    check_rlwe_plaintext_relation(
        params,
        public_commit,
        private_commit,
        delta_enc,
        decoded_q,
        r1,
        r2
    );
    if (profile_enabled) {
        const auto now = ProfileClock::now();
        prove_profile_parts.emplace_back("relation", profile_elapsed_ms(prove_stage_start, now));
        prove_stage_start = now;
    }

    std::vector<bool> w(n, false);

    auto u64cast = [](int64_t x) { return static_cast<uint64_t>(x); };

    size_t pos = 0;
    for (int i = static_cast<int>(d) - 1; i >= 0; i--)
        w[pos++] = (u64cast(private_commit.r[i]) & 1) != 0;
    for (size_t i = 0; i < k; i++)
        for (size_t j = 0; j < step_t; j++)
            w[pos++] = ((u64cast(private_commit.m[i]) >> j) & 1) != 0;
    for (size_t i = 0; i < d; i++)
        for (size_t j = 0; j < step_b; j++)
            w[pos++] = ((u64cast(private_commit.e1[i]) >> j) & 1) != 0;
    for (size_t i = 0; i < k; i++)
        for (size_t j = 0; j < step_b; j++)
            w[pos++] = ((u64cast(private_commit.e2[i]) >> j) & 1) != 0;
    for (size_t i = 0; i < d; i++)
        for (size_t j = 0; j < step_br; j++)
            w[pos++] = ((u64cast(r1[i]) >> j) & 1) != 0;
    for (size_t i = 0; i < r2.size(); i++)
        for (size_t j = 0; j < step_br; j++)
            w[pos++] = ((u64cast(r2[i]) >> j) & 1) != 0;
    if (profile_enabled) {
        const auto now = ProfileClock::now();
        prove_profile_parts.emplace_back("witness_bits", profile_elapsed_ms(prove_stage_start, now));
        prove_stage_start = now;
    }

    std::vector<uint8_t> gamma_input;
    gamma_input.insert(gamma_input.end(),
        params.hash_gamma_ds.begin(), params.hash_gamma_ds.end());
    gamma_input.insert(gamma_input.end(), seed.begin(), seed.end());
    auto gamma_vec = hash_to_fr_vec(gamma_input, 2);
    blst_fr gamma_fr   = gamma_vec[0];
    blst_fr gamma_y_fr = gamma_vec[1];

    blst_p1_affine* g_list    = static_cast<blst_p1_affine*>(params.g_list_ptr);
    blst_p2_affine* g_hat_list = static_cast<blst_p2_affine*>(params.g_hat_list_ptr);
    const blst_p1_affine* g1_gen = blst_p1_affine_generator();
    const blst_p2_affine* g2_gen = blst_p2_affine_generator();

    blst_p2 c_hat;
    blst_p2_from_affine(&c_hat, g2_gen);
    { blst_scalar s = fr_to_scalar(gamma_fr); blst_p2_mult(&c_hat, &c_hat, s.b, 256); }
    std::vector<const blst_p2_affine*> c_hat_points;
    c_hat_points.reserve(big_d);
    for (size_t j = 0; j < big_d; j++) {
        if (!w[j]) continue;
        c_hat_points.push_back(&g_hat_list[j]);
    }
    if (!c_hat_points.empty()) {
        blst_p2 c_hat_terms = p2_add_affines(c_hat_points);
        blst_p2_add(&c_hat, &c_hat, &c_hat_terms);
    }
    Proof proof;
    blst_p2_to_affine(static_cast<blst_p2_affine*>(proof.c_hat_ptr), &c_hat);
    auto c_hat_bytes = p2_to_bytes(c_hat);
    if (profile_enabled) {
        const auto now = ProfileClock::now();
        prove_profile_parts.emplace_back("c_hat", profile_elapsed_ms(prove_stage_start, now));
        prove_stage_start = now;
    }

    std::vector<uint8_t> x_bytes = build_statement_x_bytes(params, public_commit);

    std::vector<uint8_t> y_input;
    y_input.insert(y_input.end(), params.hash_ds.begin(), params.hash_ds.end());
    y_input.insert(y_input.end(), params.sid.begin(), params.sid.end());
    y_input.insert(y_input.end(), metadata.begin(), metadata.end());
    y_input.insert(y_input.end(), x_bytes.begin(), x_bytes.end());
    y_input.insert(y_input.end(), c_hat_bytes.begin(), c_hat_bytes.end());
    auto y = hash_to_fr_vec(y_input, n);

    blst_p1 c_y;
    blst_p1_from_affine(&c_y, g1_gen);
    { blst_scalar s = fr_to_scalar(gamma_y_fr); blst_p1_mult(&c_y, &c_y, s.b, 256); }
    blst_p1 c_y_terms = p1_msm_reverse_y_w(g_list, n, y, w, big_d);
    blst_p1_add(&c_y, &c_y, &c_y_terms);
    blst_p1_to_affine(static_cast<blst_p1_affine*>(proof.c_y_ptr), &c_y);
    auto c_y_bytes = p1_to_bytes(c_y);
    if (profile_enabled) {
        const auto now = ProfileClock::now();
        prove_profile_parts.emplace_back("hash_y_and_c_y", profile_elapsed_ms(prove_stage_start, now));
        prove_stage_start = now;
    }

    std::vector<uint8_t> theta_input;
    theta_input.insert(theta_input.end(), params.hash_lmap_ds.begin(), params.hash_lmap_ds.end());
    theta_input.insert(theta_input.end(), params.sid.begin(), params.sid.end());
    theta_input.insert(theta_input.end(), metadata.begin(), metadata.end());
    theta_input.insert(theta_input.end(), x_bytes.begin(), x_bytes.end());
    theta_input.insert(theta_input.end(), c_hat_bytes.begin(), c_hat_bytes.end());
    theta_input.insert(theta_input.end(), c_y_bytes.begin(), c_y_bytes.end());
    auto theta_all = hash_to_fr_vec(theta_input, d + k + m_bit_constraints + 1);
    std::vector<blst_fr> theta0(theta_all.begin(), theta_all.begin() + d + k);
    std::vector<blst_fr> m_bit_theta(
        theta_all.begin() + d + k,
        theta_all.begin() + d + k + m_bit_constraints
    );
    blst_fr delta_theta = theta_all[d + k + m_bit_constraints];
    if (fr_is_zero(delta_theta)) {
        throw std::runtime_error("Degenerate Fiat-Shamir challenge: delta_theta = 0");
    }

    std::vector<blst_fr> a_theta(big_d, fr_zero());
    compute_a_theta_fn(theta0, m_bit_theta, d, k, public_commit.a, public_commit.b,
                       a_theta, row_offset, effective_t, delta_enc,
                       params.b, params.b_r, decoded_q);

    blst_fr t_theta = fr_zero();
    for (size_t i = 0; i < d; i++)
        t_theta = fr_add(t_theta, fr_mul(theta0[i], fr_from_i64(public_commit.c1[i])));
    for (size_t i = 0; i < k; i++)
        t_theta = fr_add(t_theta, fr_mul(theta0[d + i], fr_from_i64(public_commit.c2[i])));
    if (profile_enabled) {
        const auto now = ProfileClock::now();
        prove_profile_parts.emplace_back("theta_and_a_theta", profile_elapsed_ms(prove_stage_start, now));
        prove_stage_start = now;
    }

    auto y_bytes = fr_vec_to_bytes(y);
    std::vector<uint8_t> t_input;
    t_input.insert(t_input.end(), params.hash_t_ds.begin(), params.hash_t_ds.end());
    t_input.insert(t_input.end(), params.sid.begin(), params.sid.end());
    t_input.insert(t_input.end(), metadata.begin(), metadata.end());
    t_input.insert(t_input.end(), y_bytes.begin(), y_bytes.end());
    t_input.insert(t_input.end(), x_bytes.begin(), x_bytes.end());
    t_input.insert(t_input.end(), c_hat_bytes.begin(), c_hat_bytes.end());
    t_input.insert(t_input.end(), c_y_bytes.begin(), c_y_bytes.end());
    auto t_vec = hash_to_fr_vec_128bit(t_input, n);

    std::vector<uint8_t> agg_input;
    agg_input.insert(agg_input.end(), params.hash_agg_ds.begin(), params.hash_agg_ds.end());
    agg_input.insert(agg_input.end(), params.sid.begin(), params.sid.end());
    agg_input.insert(agg_input.end(), metadata.begin(), metadata.end());
    agg_input.insert(agg_input.end(), x_bytes.begin(), x_bytes.end());
    agg_input.insert(agg_input.end(), c_hat_bytes.begin(), c_hat_bytes.end());
    agg_input.insert(agg_input.end(), c_y_bytes.begin(), c_y_bytes.end());
    auto delta_agg = hash_to_fr_vec(agg_input, 2);
    blst_fr delta_eq = delta_agg[0];
    blst_fr delta_y  = delta_agg[1];
    if (fr_is_zero(delta_eq) || fr_is_zero(delta_y)) {
        throw std::runtime_error("Degenerate Fiat-Shamir challenge: delta_eq/delta_y = 0");
    }
    if (profile_enabled) {
        const auto now = ProfileClock::now();
        prove_profile_parts.emplace_back("hash_t_and_agg", profile_elapsed_ms(prove_stage_start, now));
        prove_stage_start = now;
    }

    std::vector<blst_fr> poly_1(big_d + 1, fr_zero());
    poly_1[0] = gamma_fr;
    for (size_t j = 0; j < big_d; j++)
        poly_1[j + 1] = fr_from_u64(w[j] ? 1 : 0);

    std::vector<blst_fr> poly_0(n + 1, fr_zero());
    poly_0[0] = fr_mul(delta_y, gamma_y_fr);
    for (size_t i = 1; i <= n; i++) {
        blst_fr yi = y[i - 1];
        blst_fr ti = t_vec[i - 1];
        blst_fr wi = fr_from_u64((i <= big_d && w[i - 1]) ? 1 : 0);
        blst_fr val = fr_add(
            fr_mul(delta_y, fr_mul(yi, wi)),
            fr_mul(fr_sub(fr_mul(delta_eq, ti), delta_y), yi)
        );
        if (i <= big_d)
            val = fr_add(val, fr_mul(delta_theta, a_theta[i - 1]));
        poly_0[n + 1 - i] = val;
    }

    std::vector<blst_fr> poly_2(n + 1, fr_zero());
    poly_2[0] = gamma_y_fr;
    for (size_t i = 1; i <= big_d; i++) {
        blst_fr yi = y[i - 1];
        blst_fr wi = fr_from_u64(w[i - 1] ? 1 : 0);
        poly_2[n + 1 - i] = fr_mul(yi, wi);
    }

    std::vector<blst_fr> poly_3(n + 1, fr_zero());
    for (size_t i = 1; i <= n; i++)
        poly_3[i] = fr_mul(delta_eq, t_vec[i - 1]);

    auto poly = poly_sub_fr(poly_mul_fr(poly_0, poly_1),
                            poly_mul_fr(poly_2, poly_3));
    if (poly.size() > n + 1) {
        poly[n + 1] = fr_sub(poly[n + 1], fr_mul(t_theta, delta_theta));
    }
    if (profile_enabled) {
        const auto now = ProfileClock::now();
        prove_profile_parts.emplace_back("poly_build_and_mul", profile_elapsed_ms(prove_stage_start, now));
        prove_stage_start = now;
    }

    blst_p1 pi_point;
    blst_p1_from_affine(&pi_point, g1_gen);
    { blst_scalar s = fr_to_scalar(poly[0]); blst_p1_mult(&pi_point, &pi_point, s.b, 256); }
    const size_t pi_terms = std::min(poly.size() - 1, 2 * n);
    blst_p1 pi_tail = p1_msm_direct(g_list, poly.data() + 1, pi_terms);
    blst_p1_add(&pi_point, &pi_point, &pi_tail);
    blst_p1_to_affine(static_cast<blst_p1_affine*>(proof.pi_ptr), &pi_point);
    if (profile_enabled) {
        const auto now = ProfileClock::now();
        prove_profile_parts.emplace_back("pi_msm", profile_elapsed_ms(prove_stage_start, now));
        emit_profile_summary("prove", prove_profile_parts);
        (void)prove_start_total;
    }

    return proof;
}

template<typename P>
bool verify(
    const Proof& proof,
    const PublicParams<P>& params,
    const PublicCommit<P>& public_commit,
    const std::vector<uint8_t>& metadata
) {
    const bool profile_enabled = zk_profile_enabled();
    auto verify_stage_start = profile_enabled ? ProfileClock::now() : ProfileClock::time_point{};
    std::vector<std::pair<const char*, double>> verify_profile_parts;

    if (!public_params_basic_usable(params)) {
        return false;
    }
    if (!params.integrity_checked && !public_params_usable(params)) {
        return false;
    }

    const size_t d = params.d;
    const size_t k = public_commit.c2.size();

    if (!proof_points_usable(proof)) {
        return false;
    }

    if (k == 0 || k > params.k ||
        public_commit.a.size() != d || public_commit.b.size() != d ||
        public_commit.c1.size() != d || public_commit.c2.size() != k) {
        return false;
    }

    __uint128_t decoded_q   = decode_q(params.q);
    uint64_t    delta_enc   = static_cast<uint64_t>(decoded_q / params.t);
    uint64_t    effective_t = params.t >> params.msbs_zero_padding_bit_count;

    size_t step_t  = ilog2_u64(effective_t);
    size_t step_b  = 1 + ilog2_u64(params.b);
    size_t step_br = 1 + ilog2_u64(params.b_r);
    size_t m_bit_constraints = k * (step_t - 1);
    size_t big_d   = d + k * step_t + (d + k) * (step_b + step_br);
    size_t n       = params.n;
    if (big_d + 1 > n) {
        return false;
    }

    const size_t row_offset = d - k;

    blst_p2 c_hat;
    blst_p2_from_affine(&c_hat, static_cast<const blst_p2_affine*>(proof.c_hat_ptr));
    blst_p1 c_y_p;
    blst_p1_from_affine(&c_y_p, static_cast<const blst_p1_affine*>(proof.c_y_ptr));
    blst_p1 pi_p;
    blst_p1_from_affine(&pi_p, static_cast<const blst_p1_affine*>(proof.pi_ptr));

    auto c_hat_bytes = p2_to_bytes(c_hat);
    auto c_y_bytes   = p1_to_bytes(c_y_p);

    blst_p1_affine* g_list     = static_cast<blst_p1_affine*>(params.g_list_ptr);
    blst_p2_affine* g_hat_list = static_cast<blst_p2_affine*>(params.g_hat_list_ptr);
    const blst_p2_affine* g2_gen = blst_p2_affine_generator();

    std::vector<uint8_t> x_bytes = build_statement_x_bytes(params, public_commit);

    std::vector<uint8_t> y_input;
    y_input.insert(y_input.end(), params.hash_ds.begin(), params.hash_ds.end());
    y_input.insert(y_input.end(), params.sid.begin(), params.sid.end());
    y_input.insert(y_input.end(), metadata.begin(), metadata.end());
    y_input.insert(y_input.end(), x_bytes.begin(), x_bytes.end());
    y_input.insert(y_input.end(), c_hat_bytes.begin(), c_hat_bytes.end());
    auto y = hash_to_fr_vec(y_input, n);
    if (profile_enabled) {
        const auto now = ProfileClock::now();
        verify_profile_parts.emplace_back("hash_y", profile_elapsed_ms(verify_stage_start, now));
        verify_stage_start = now;
    }

    std::vector<uint8_t> theta_input;
    theta_input.insert(theta_input.end(), params.hash_lmap_ds.begin(), params.hash_lmap_ds.end());
    theta_input.insert(theta_input.end(), params.sid.begin(), params.sid.end());
    theta_input.insert(theta_input.end(), metadata.begin(), metadata.end());
    theta_input.insert(theta_input.end(), x_bytes.begin(), x_bytes.end());
    theta_input.insert(theta_input.end(), c_hat_bytes.begin(), c_hat_bytes.end());
    theta_input.insert(theta_input.end(), c_y_bytes.begin(), c_y_bytes.end());
    auto theta_all = hash_to_fr_vec(theta_input, d + k + m_bit_constraints + 1);
    std::vector<blst_fr> theta0(theta_all.begin(), theta_all.begin() + d + k);
    std::vector<blst_fr> m_bit_theta(
        theta_all.begin() + d + k,
        theta_all.begin() + d + k + m_bit_constraints
    );
    blst_fr delta_theta = theta_all[d + k + m_bit_constraints];
    if (fr_is_zero(delta_theta)) {
        return false;
    }

    std::vector<blst_fr> a_theta(big_d, fr_zero());
    compute_a_theta_fn(theta0, m_bit_theta, d, k, public_commit.a, public_commit.b,
                       a_theta, row_offset, effective_t, delta_enc,
                       params.b, params.b_r, decoded_q);

    blst_fr t_theta = fr_zero();
    for (size_t i = 0; i < d; i++)
        t_theta = fr_add(t_theta, fr_mul(theta0[i], fr_from_i64(public_commit.c1[i])));
    for (size_t i = 0; i < k; i++)
        t_theta = fr_add(t_theta, fr_mul(theta0[d + i], fr_from_i64(public_commit.c2[i])));
    if (profile_enabled) {
        const auto now = ProfileClock::now();
        verify_profile_parts.emplace_back("theta_and_a_theta", profile_elapsed_ms(verify_stage_start, now));
        verify_stage_start = now;
    }

    auto y_bytes = fr_vec_to_bytes(y);
    std::vector<uint8_t> t_input;
    t_input.insert(t_input.end(), params.hash_t_ds.begin(), params.hash_t_ds.end());
    t_input.insert(t_input.end(), params.sid.begin(), params.sid.end());
    t_input.insert(t_input.end(), metadata.begin(), metadata.end());
    t_input.insert(t_input.end(), y_bytes.begin(), y_bytes.end());
    t_input.insert(t_input.end(), x_bytes.begin(), x_bytes.end());
    t_input.insert(t_input.end(), c_hat_bytes.begin(), c_hat_bytes.end());
    t_input.insert(t_input.end(), c_y_bytes.begin(), c_y_bytes.end());
    auto t_vec = hash_to_fr_vec_128bit(t_input, n);

    std::vector<uint8_t> agg_input;
    agg_input.insert(agg_input.end(), params.hash_agg_ds.begin(), params.hash_agg_ds.end());
    agg_input.insert(agg_input.end(), params.sid.begin(), params.sid.end());
    agg_input.insert(agg_input.end(), metadata.begin(), metadata.end());
    agg_input.insert(agg_input.end(), x_bytes.begin(), x_bytes.end());
    agg_input.insert(agg_input.end(), c_hat_bytes.begin(), c_hat_bytes.end());
    agg_input.insert(agg_input.end(), c_y_bytes.begin(), c_y_bytes.end());
    auto delta_agg = hash_to_fr_vec(agg_input, 2);
    blst_fr delta_eq = delta_agg[0];
    blst_fr delta_y  = delta_agg[1];
    if (fr_is_zero(delta_eq) || fr_is_zero(delta_y)) {
        return false;
    }
    if (profile_enabled) {
        const auto now = ProfileClock::now();
        verify_profile_parts.emplace_back("hash_t_and_agg", profile_elapsed_ms(verify_stage_start, now));
        verify_stage_start = now;
    }

    blst_p1 term0_g1 = c_y_p;
    { blst_scalar s = fr_to_scalar(delta_y); blst_p1_mult(&term0_g1, &term0_g1, s.b, 256); }
    blst_p1 term0_tail = p1_msm_verify_term0(
        g_list, n, big_d, y, t_vec, a_theta, delta_eq, delta_y, delta_theta
    );
    blst_p1_add(&term0_g1, &term0_g1, &term0_tail);

    blst_p2 term1_g2 = p2_msm_t_vec(g_hat_list, n, t_vec, delta_eq);
    if (profile_enabled) {
        const auto now = ProfileClock::now();
        verify_profile_parts.emplace_back("verify_msms", profile_elapsed_ms(verify_stage_start, now));
        verify_stage_start = now;
    }

    blst_fr scale = fr_mul(t_theta, delta_theta);

    blst_p1 scaled_p; 
    blst_p1_from_affine(&scaled_p, &g_list[0]); 
    blst_scalar s_scale = fr_to_scalar(scale);
    blst_p1_mult(&scaled_p, &scaled_p, s_scale.b, 256);

    blst_p1_affine pi_aff;           blst_p1_to_affine(&pi_aff, &pi_p);
    blst_p1_affine term0_g1_aff;     blst_p1_to_affine(&term0_g1_aff, &term0_g1);
    blst_p1_affine c_y_aff;          blst_p1_to_affine(&c_y_aff, &c_y_p);
    blst_p1_affine scaled_p_aff;     blst_p1_to_affine(&scaled_p_aff, &scaled_p);

    blst_p2_affine c_hat_aff  = *static_cast<const blst_p2_affine*>(proof.c_hat_ptr);
    blst_p2_affine term1_g2_aff;     blst_p2_to_affine(&term1_g2_aff, &term1_g2);

    blst_p1 neg_term0 = term0_g1; 
    blst_p1_cneg(&neg_term0, 1);
    blst_p1_affine neg_term0_aff; 
    blst_p1_to_affine(&neg_term0_aff, &neg_term0);

    blst_fp12 ml_pi, ml_neg_t0, ml_cy, ml_scale;

    blst_miller_loop(&ml_pi,     g2_gen,         &pi_aff);
    blst_miller_loop(&ml_neg_t0, &c_hat_aff,     &neg_term0_aff);
    blst_miller_loop(&ml_cy, &term1_g2_aff, &c_y_aff);
    blst_miller_loop(&ml_scale, &g_hat_list[n - 1], &scaled_p_aff);

    auto is_one = [&](const blst_fp12& x) {
        blst_fp12 tmp = x;
        blst_fp12 fe;
        blst_final_exp(&fe, &tmp);
#if defined(AUTOHOG_ZK_BACKEND_mcl)
        return fe.isOne();
#else
        return blst_fp12_is_equal(&fe, &pairing_identity_gt());
#endif
    };

    blst_fp12 acc = ml_pi;
    blst_fp12_mul(&acc, &acc, &ml_neg_t0);
    blst_fp12_mul(&acc, &acc, &ml_cy);
    blst_fp12_mul(&acc, &acc, &ml_scale);
    const bool ok = is_one(acc);
    if (profile_enabled) {
        const auto now = ProfileClock::now();
        verify_profile_parts.emplace_back("pairing", profile_elapsed_ms(verify_stage_start, now));
        emit_profile_summary("verify", verify_profile_parts);
    }
    return ok;
}

template<typename P>
void save_public_params(const PublicParams<P>& params, std::ostream& ofs) {
    ofs.write(kCrsMagicCurrent, 4);

    const std::uint8_t backend_tag = static_cast<std::uint8_t>(autohog::zk::config::backend_id);
    const std::uint8_t curve_tag = static_cast<std::uint8_t>(autohog::zk::config::curve_id);
    const std::uint16_t reserved = 0;
    const std::uint32_t p1_storage_size = static_cast<std::uint32_t>(sizeof(blst_p1_affine));
    const std::uint32_t p2_storage_size = static_cast<std::uint32_t>(sizeof(blst_p2_affine));
    ofs.write(reinterpret_cast<const char*>(&backend_tag), sizeof(backend_tag));
    ofs.write(reinterpret_cast<const char*>(&curve_tag), sizeof(curve_tag));
    ofs.write(reinterpret_cast<const char*>(&reserved), sizeof(reserved));
    ofs.write(reinterpret_cast<const char*>(&p1_storage_size), sizeof(p1_storage_size));
    ofs.write(reinterpret_cast<const char*>(&p2_storage_size), sizeof(p2_storage_size));

    ofs.write(reinterpret_cast<const char*>(&params.n), sizeof(params.n));
    ofs.write(reinterpret_cast<const char*>(&params.d), sizeof(params.d));
    ofs.write(reinterpret_cast<const char*>(&params.k), sizeof(params.k));
    ofs.write(reinterpret_cast<const char*>(&params.b), sizeof(params.b));
    ofs.write(reinterpret_cast<const char*>(&params.b_r), sizeof(params.b_r));
    ofs.write(reinterpret_cast<const char*>(&params.q), sizeof(params.q));
    ofs.write(reinterpret_cast<const char*>(&params.t), sizeof(params.t));
    ofs.write(reinterpret_cast<const char*>(&params.msbs_zero_padding_bit_count),
              sizeof(params.msbs_zero_padding_bit_count));

    ofs.write(reinterpret_cast<const char*>(params.sid.data()), 16);
    ofs.write(reinterpret_cast<const char*>(params.hash_ds.data()), 16);
    ofs.write(reinterpret_cast<const char*>(params.hash_t_ds.data()), 16);
    ofs.write(reinterpret_cast<const char*>(params.hash_agg_ds.data()), 16);
    ofs.write(reinterpret_cast<const char*>(params.hash_lmap_ds.data()), 16);
    ofs.write(reinterpret_cast<const char*>(params.hash_z_ds.data()), 16);
    ofs.write(reinterpret_cast<const char*>(params.hash_w_ds.data()), 16);
    ofs.write(reinterpret_cast<const char*>(params.hash_gamma_ds.data()), 16);

    const blst_p1_affine* g_list = static_cast<const blst_p1_affine*>(params.g_list_ptr);
    ofs.write(reinterpret_cast<const char*>(g_list),
              (2 * params.n) * sizeof(blst_p1_affine));

    const blst_p2_affine* g_hat_list = static_cast<const blst_p2_affine*>(params.g_hat_list_ptr);
    ofs.write(reinterpret_cast<const char*>(g_hat_list),
              params.n * sizeof(blst_p2_affine));

    
    auto pos = ofs.tellp();
    std::cout << "  CRS saved " << ")\n";
}

template<typename P>
void save_public_params(const PublicParams<P>& params, const std::string& filename) {
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) throw std::runtime_error("Cannot open CRS file for writing: " + filename);
    save_public_params(params, ofs);
}

template<typename P>
void save_public_params_to_memory(const PublicParams<P>& params, std::string& data) {
    std::ostringstream ofs;
    save_public_params(params, ofs);
    data = ofs.str();
}

template<typename P>
PublicParams<P> load_public_params_impl(std::istream& ifs, bool validate_structure)
{
    char magic[4];
    ifs.read(magic, 4);
    const bool is_current = std::memcmp(magic, kCrsMagicCurrent, 4) == 0;
    const bool is_legacy = std::memcmp(magic, kCrsMagicLegacy, 4) == 0;
    if (!is_current && !is_legacy) {
        throw std::runtime_error("Invalid CRS file format (bad magic)");
    }

    PublicParams<P> params;

    std::uint32_t p1_storage_size = static_cast<std::uint32_t>(sizeof(blst_p1_affine));
    std::uint32_t p2_storage_size = static_cast<std::uint32_t>(sizeof(blst_p2_affine));

    if (is_current) {
        std::uint8_t backend_tag = 0;
        std::uint8_t curve_tag = 0;
        std::uint16_t reserved = 0;

        ifs.read(reinterpret_cast<char*>(&backend_tag), sizeof(backend_tag));
        ifs.read(reinterpret_cast<char*>(&curve_tag), sizeof(curve_tag));
        ifs.read(reinterpret_cast<char*>(&reserved), sizeof(reserved));
        ifs.read(reinterpret_cast<char*>(&p1_storage_size), sizeof(p1_storage_size));
        ifs.read(reinterpret_cast<char*>(&p2_storage_size), sizeof(p2_storage_size));

        if (backend_tag != static_cast<std::uint8_t>(autohog::zk::config::backend_id) ||
            curve_tag != static_cast<std::uint8_t>(autohog::zk::config::curve_id)) {
            throw std::runtime_error(
                "CRS backend/curve mismatch: file targets " +
                std::to_string(backend_tag) + "/" + std::to_string(curve_tag) +
                " but build expects " +
                std::to_string(static_cast<std::uint8_t>(autohog::zk::config::backend_id)) + "/" +
                std::to_string(static_cast<std::uint8_t>(autohog::zk::config::curve_id))
            );
        }
        if (p1_storage_size != sizeof(blst_p1_affine) ||
            p2_storage_size != sizeof(blst_p2_affine)) {
            throw std::runtime_error("CRS point storage layout mismatch");
        }
        (void)reserved;
    } else if (autohog::zk::config::curve_id != autohog::zk::config::CurveId::bls12_381) {
        throw std::runtime_error("Legacy ZKC2 CRS files are only supported for bls12_381 builds");
    }

    ifs.read(reinterpret_cast<char*>(&params.n), sizeof(params.n));
    ifs.read(reinterpret_cast<char*>(&params.d), sizeof(params.d));
    ifs.read(reinterpret_cast<char*>(&params.k), sizeof(params.k));
    ifs.read(reinterpret_cast<char*>(&params.b), sizeof(params.b));
    ifs.read(reinterpret_cast<char*>(&params.b_r), sizeof(params.b_r));
    ifs.read(reinterpret_cast<char*>(&params.q), sizeof(params.q));
    ifs.read(reinterpret_cast<char*>(&params.t), sizeof(params.t));
    ifs.read(reinterpret_cast<char*>(&params.msbs_zero_padding_bit_count),
             sizeof(params.msbs_zero_padding_bit_count));

    if (!ifs) {
        throw std::runtime_error("CRS file truncated or corrupted in parameter header");
    }
    if (!public_params_dimensions_usable(params)) {
        throw std::runtime_error(
            "CRS parameter dimensions are invalid or inconsistent with this build");
    }

    ifs.read(reinterpret_cast<char*>(params.sid.data()), 16);
    ifs.read(reinterpret_cast<char*>(params.hash_ds.data()), 16);
    ifs.read(reinterpret_cast<char*>(params.hash_t_ds.data()), 16);
    ifs.read(reinterpret_cast<char*>(params.hash_agg_ds.data()), 16);
    ifs.read(reinterpret_cast<char*>(params.hash_lmap_ds.data()), 16);
    ifs.read(reinterpret_cast<char*>(params.hash_z_ds.data()), 16);
    ifs.read(reinterpret_cast<char*>(params.hash_w_ds.data()), 16);
    ifs.read(reinterpret_cast<char*>(params.hash_gamma_ds.data()), 16);

    params.g_list_ptr = new blst_p1_affine[2 * params.n];
    ifs.read(reinterpret_cast<char*>(params.g_list_ptr),
             (2 * params.n) * sizeof(blst_p1_affine));

    params.g_hat_list_ptr = new blst_p2_affine[params.n];
    ifs.read(reinterpret_cast<char*>(params.g_hat_list_ptr),
             params.n * sizeof(blst_p2_affine));

    if (!ifs) {
        throw std::runtime_error("CRS file truncated or corrupted");
    }

    if (validate_structure && !public_params_usable(params)) {
        throw std::runtime_error(
            "CRS contains invalid group elements or an inconsistent power chain");
    }
    params.integrity_checked = validate_structure;

    return params;
}

template<typename P>
PublicParams<P> load_public_params(const std::string& filename)
{
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) throw std::runtime_error("Cannot open CRS file for reading: " + filename);
    return load_public_params_impl<P>(ifs, true);
}

template<typename P>
PublicParams<P> load_public_params_from_memory(const char* crs, int len)
{
    constexpr int kMaxCrsBufferBytes = 512 * 1024 * 1024;
    if (crs == nullptr || len <= 0) {
        throw std::invalid_argument("CRS memory buffer must be non-null and non-empty");
    }
    if (len > kMaxCrsBufferBytes) {
        throw std::length_error("CRS memory buffer exceeds the configured limit");
    }
    std::string buffer(crs, crs + len);
    std::istringstream iss(buffer);
    if (!iss) throw std::runtime_error("Cannot open CRS stream for reading " );
    return load_public_params_impl<P>(iss, true);
}

template PublicParams<TFHEpp::lvl2param> load_public_params_from_memory<TFHEpp::lvl2param>(const char* crs, int len);

template<typename P>
PublicParams<P> load_public_params_from_memory_by_expected_hash(
    const char* crs,
    int len,
    const std::string& expected_hash_hex
) {
    if (crs == nullptr || len <= 0) {
        throw std::invalid_argument("CRS memory buffer must be non-null and non-empty");
    }
    const std::string expected_hash = normalize_sha256_hex(
        expected_hash_hex,
        "expected CRS SHA-256 digest"
    );
    const std::string actual_hash = compute_memory_sha256_hex(
        crs,
        static_cast<size_t>(len)
    );
    if (actual_hash != expected_hash) {
        throw std::runtime_error(
            "CRS SHA-256 mismatch for memory buffer: expected " + expected_hash +
            ", got " + actual_hash
        );
    }
    return load_public_params_from_memory<P>(crs, len);
}

template PublicParams<TFHEpp::lvl2param>
load_public_params_from_memory_by_expected_hash<TFHEpp::lvl2param>(
    const char* crs,
    int len,
    const std::string& expected_hash_hex
);

template<typename P>
PublicParams<P> load_public_params_by_hash(
    const std::string& filename,
    const std::string& hash_filename
) {
    std::ifstream ifs(hash_filename, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Cannot open CRS hash file: " + hash_filename);
    }

    std::string expected_hash;
    std::getline(ifs, expected_hash);
    expected_hash = normalize_sha256_hex(expected_hash, "hash file: " + hash_filename);

    return load_public_params_by_expected_hash<P>(filename, expected_hash);
}

template<typename P>
PublicParams<P> load_public_params_by_expected_hash(
    const std::string& filename,
    const std::string& expected_hash_hex
) {
    const std::string expected_hash = normalize_sha256_hex(
        expected_hash_hex,
        "expected CRS SHA-256 digest"
    );

    constexpr std::uintmax_t kMaxCrsFileBytes = std::uintmax_t{512} * 1024 * 1024;
    std::ifstream ifs(filename, std::ios::binary | std::ios::ate);
    if (!ifs) {
        throw std::runtime_error("Cannot open CRS file for reading: " + filename);
    }
    const std::streamoff file_size = ifs.tellg();
    if (file_size <= 0 || static_cast<std::uintmax_t>(file_size) > kMaxCrsFileBytes) {
        throw std::runtime_error("CRS file size is invalid or exceeds the configured limit");
    }
    ifs.seekg(0, std::ios::beg);
    std::string buffer(static_cast<size_t>(file_size), '\0');
    ifs.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    if (!ifs || ifs.gcount() != static_cast<std::streamsize>(buffer.size())) {
        throw std::runtime_error("Failed reading CRS file: " + filename);
    }

    const std::string actual_hash = compute_memory_sha256_hex(buffer.data(), buffer.size());
    if (actual_hash != expected_hash) {
        throw std::runtime_error(
            "CRS SHA-256 mismatch for " + filename +
            ": expected " + expected_hash + ", got " + actual_hash
        );
    }

    std::istringstream iss(buffer, std::ios::in | std::ios::binary);
    return load_public_params_impl<P>(iss, true);
}

template PublicParams<TFHEpp::lvl2param> generate_public_params<TFHEpp::lvl2param>(
    uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t,
    const std::array<uint8_t, SEED_LEN>&
);

template PublicCommit<TFHEpp::lvl2param> make_public_commit<TFHEpp::lvl2param>(
    const PKE::SeededRlwePkPair<TFHEpp::lvl2param>&,
    const TFHEpp::TRLWE<TFHEpp::lvl2param>&,
    uint64_t
);

template PrivateCommit<TFHEpp::lvl2param> make_private_commit<TFHEpp::lvl2param>(
    const std::vector<int64_t>&,
    const std::vector<int64_t>&,
    const std::vector<int64_t>&,
    const std::vector<int64_t>&
);

template Proof prove<TFHEpp::lvl2param>(
    const PublicParams<TFHEpp::lvl2param>&,
    const PublicCommit<TFHEpp::lvl2param>&,
    const PrivateCommit<TFHEpp::lvl2param>&,
    const std::vector<uint8_t>&,
    ComputeLoad,
    const std::array<uint8_t, SEED_LEN>&
);

template bool verify<TFHEpp::lvl2param>(
    const Proof&,
    const PublicParams<TFHEpp::lvl2param>&,
    const PublicCommit<TFHEpp::lvl2param>&,
    const std::vector<uint8_t>&
);

template void save_public_params<TFHEpp::lvl2param>(
    const PublicParams<TFHEpp::lvl2param>&,
    const std::string&
);

template void save_public_params_to_memory<TFHEpp::lvl2param>(
    const PublicParams<TFHEpp::lvl2param>&,
    std::string&
);

template PublicParams<TFHEpp::lvl2param> load_public_params<TFHEpp::lvl2param>(
    const std::string&
);

template PublicParams<TFHEpp::lvl2param> load_public_params_by_hash<TFHEpp::lvl2param>(
    const std::string&,
    const std::string&
);

template PublicParams<TFHEpp::lvl2param> load_public_params_by_expected_hash<TFHEpp::lvl2param>(
    const std::string&,
    const std::string&
);

template bool validate_public_params<TFHEpp::lvl2param>(
    const PublicParams<TFHEpp::lvl2param>&
);

}
