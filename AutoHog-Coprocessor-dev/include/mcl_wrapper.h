#pragma once

#include <cstdint>
#include <cstring>
#include <mutex>
#include <stdexcept>
#include <vector>

#include <mcl/bn.hpp>
#include <mcl/op.hpp>

#include "zk_config.h"

using blst_fr = mcl::bn::Fr;
using blst_p1 = mcl::bn::G1;
using blst_p1_affine = mcl::bn::G1;
using blst_p2 = mcl::bn::G2;
using blst_p2_affine = mcl::bn::G2;
using blst_fp12 = mcl::bn::GT;
using limb_t = std::uint64_t;

struct blst_scalar {
    std::uint8_t b[32];
};

enum BLST_ERROR {
    BLST_SUCCESS = 0,
    BLST_BAD_ENCODING = 1,
    BLST_POINT_NOT_IN_GROUP = 2,
};

namespace autohog::zk::mcl_backend::detail {

inline const mcl::bn::CurveParam& curve_param()
{
    if constexpr (config::curve_id == config::CurveId::bn_snark1) {
        printf("USE mcl::BN_SNARK1\n");
        return mcl::BN_SNARK1;
    } else {
        printf("USE mcl::BLS12_381\n");
        return mcl::BLS12_381;
    }
}

inline constexpr size_t g1_serialized_size()
{
    return config::g1_serialized_size;
}

inline constexpr size_t g2_serialized_size()
{
    return config::g2_serialized_size;
}

inline void ensure_initialized()
{
    static std::once_flag once;
    std::call_once(once, [] {
        mcl::bn::initPairing(curve_param());
        if constexpr (config::curve_id == config::CurveId::bls12_381) {
          mcl::bn::Fp::setETHserialization(true);
          mcl::bn::Fr::setETHserialization(false);
        }

        // Disable G1 subgroup verification for BN254 since its G1 cofactor is 1.
        constexpr bool _verifyOrderG1 = (config::curve_id != config::CurveId::bn_snark1);
        mcl::bn::verifyOrderG1(_verifyOrderG1);
        mcl::bn::verifyOrderG2(true);
    });
}

inline blst_fr scalar_bytes_to_fr(const std::uint8_t* bytes, size_t byte_count)
{
    ensure_initialized();
    blst_fr out;
    out.setLittleEndianMod(bytes, byte_count);
    return out;
}

inline void fr_to_le32(const blst_fr& value, std::uint8_t out[32])
{
    ensure_initialized();
    const size_t written = value.serialize(out, 32, mcl::IoSerialize);
    if (written != 32) {
        throw std::runtime_error("mcl Fr serialize size mismatch");
    }
}

inline const blst_p1_affine& g1_generator_value()
{
    ensure_initialized();
    static const blst_p1_affine gen = [] {
        blst_p1_affine point;
        static constexpr char kDst[] = "AutoHog:mcl-backend:g1";
        mcl::bn::hashAndMapToG1(point, kDst, sizeof(kDst) - 1);
        if (point.isZero()) {
            throw std::runtime_error("mcl G1 generator hash mapped to zero");
        }
        point.normalize();
        return point;
    }();
    return gen;
}

inline const blst_p2_affine& g2_generator_value()
{
    ensure_initialized();
    static const blst_p2_affine gen = [] {
        blst_p2_affine point;
        static constexpr char kDst[] = "AutoHog:mcl-backend:g2";
        mcl::bn::hashAndMapToG2(point, kDst, sizeof(kDst) - 1);
        if (point.isZero()) {
            throw std::runtime_error("mcl G2 generator hash mapped to zero");
        }
        point.normalize();
        return point;
    }();
    return gen;
}

inline bool g1_in_subgroup(const blst_p1_affine& point)
{
    ensure_initialized();
    // BN254 G1 has cofactor = 1, so subgroup checks on G1 are unnecessary and can be safely removed.
    if constexpr (config::curve_id == config::CurveId::bn_snark1) {
      return point.isValid();
    }

    return point.isValid() && point.isValidOrder();
}

inline bool g2_in_subgroup(const blst_p2_affine& point)
{
    ensure_initialized();
    return point.isValid() && point.isValidOrder();
}

template <class Point, class AffinePoint>
inline void copy_points_from_ptrs(std::vector<Point>& out,
                                  const AffinePoint* const* points,
                                  size_t n)
{
    out.resize(n);
    for (size_t i = 0; i < n; ++i) {
        out[i] = *points[i];
    }
}

template <class FrVec>
inline void collect_scalars_from_ptrs(FrVec& out,
                                      const unsigned char* const* scalar_ptrs,
                                      size_t n,
                                      size_t nbits)
{
    const size_t scalar_bytes = (nbits + 7) / 8;
    out.resize(n);
    for (size_t i = 0; i < n; ++i) {
        out[i] = scalar_bytes_to_fr(scalar_ptrs[i], scalar_bytes);
    }
}

}  // namespace autohog::zk::mcl_backend::detail

inline void blst_scalar_from_lendian(blst_scalar* out, const std::uint8_t in[32])
{
    std::memcpy(out->b, in, 32);
}

inline void blst_scalar_from_fr(blst_scalar* out, const blst_fr* in)
{
    autohog::zk::mcl_backend::detail::fr_to_le32(*in, out->b);
}

inline void blst_fr_from_scalar(blst_fr* out, const blst_scalar* in)
{
    *out = autohog::zk::mcl_backend::detail::scalar_bytes_to_fr(in->b, sizeof(in->b));
}

inline void blst_fr_from_uint64(blst_fr* out, const std::uint64_t in[4])
{
    *out = autohog::zk::mcl_backend::detail::scalar_bytes_to_fr(
        reinterpret_cast<const std::uint8_t*>(in),
        sizeof(std::uint64_t) * 4
    );
}

inline void blst_fr_add(blst_fr* out, const blst_fr* a, const blst_fr* b)
{
    blst_fr::add(*out, *a, *b);
}

inline void blst_fr_sub(blst_fr* out, const blst_fr* a, const blst_fr* b)
{
    blst_fr::sub(*out, *a, *b);
}

inline void blst_fr_mul(blst_fr* out, const blst_fr* a, const blst_fr* b)
{
    blst_fr::mul(*out, *a, *b);
}

inline void blst_fr_cneg(blst_fr* out, const blst_fr* a, bool do_neg)
{
    if (do_neg) {
        blst_fr::neg(*out, *a);
    } else {
        *out = *a;
    }
}

inline void blst_fr_inverse(blst_fr* out, const blst_fr* a)
{
    blst_fr::inv(*out, *a);
}

inline const blst_p1_affine* blst_p1_affine_generator()
{
    return &autohog::zk::mcl_backend::detail::g1_generator_value();
}

inline const blst_p2_affine* blst_p2_affine_generator()
{
    return &autohog::zk::mcl_backend::detail::g2_generator_value();
}

inline void blst_p1_from_affine(blst_p1* out, const blst_p1_affine* in)
{
    autohog::zk::mcl_backend::detail::ensure_initialized();
    *out = *in;
}

inline void blst_p2_from_affine(blst_p2* out, const blst_p2_affine* in)
{
    autohog::zk::mcl_backend::detail::ensure_initialized();
    *out = *in;
}

inline void blst_p1_to_affine(blst_p1_affine* out, const blst_p1* in)
{
    autohog::zk::mcl_backend::detail::ensure_initialized();
    *out = *in;
    out->normalize();
}

inline void blst_p2_to_affine(blst_p2_affine* out, const blst_p2* in)
{
    autohog::zk::mcl_backend::detail::ensure_initialized();
    *out = *in;
    out->normalize();
}

inline void blst_p1_add(blst_p1* out, const blst_p1* a, const blst_p1* b)
{
    blst_p1::add(*out, *a, *b);
}

inline void blst_p2_add(blst_p2* out, const blst_p2* a, const blst_p2* b)
{
    blst_p2::add(*out, *a, *b);
}

inline void blst_p1_cneg(blst_p1* out, bool do_neg)
{
    if (do_neg) {
        *out = -*out;
    }
}

inline void blst_p2_cneg(blst_p2* out, bool do_neg)
{
    if (do_neg) {
        *out = -*out;
    }
}

inline void blst_p1_mult(blst_p1* out,
                         const blst_p1* point,
                         const std::uint8_t* scalar,
                         size_t nbits)
{
    const blst_fr fr =
        autohog::zk::mcl_backend::detail::scalar_bytes_to_fr(scalar, (nbits + 7) / 8);
    blst_p1::mul(*out, *point, fr);
}

inline void blst_p2_mult(blst_p2* out,
                         const blst_p2* point,
                         const std::uint8_t* scalar,
                         size_t nbits)
{
    const blst_fr fr =
        autohog::zk::mcl_backend::detail::scalar_bytes_to_fr(scalar, (nbits + 7) / 8);
    blst_p2::mul(*out, *point, fr);
}

inline bool blst_p1_affine_is_inf(const blst_p1_affine* point)
{
    return point != nullptr && point->isZero();
}

inline bool blst_p2_affine_is_inf(const blst_p2_affine* point)
{
    return point != nullptr && point->isZero();
}

inline bool blst_p1_affine_on_curve(const blst_p1_affine* point)
{
    return point != nullptr && point->isValid();
}

inline bool blst_p2_affine_on_curve(const blst_p2_affine* point)
{
    return point != nullptr && point->isValid();
}

inline bool blst_p1_affine_in_g1(const blst_p1_affine* point)
{
    return point != nullptr && autohog::zk::mcl_backend::detail::g1_in_subgroup(*point);
}

inline bool blst_p2_affine_in_g2(const blst_p2_affine* point)
{
    return point != nullptr && autohog::zk::mcl_backend::detail::g2_in_subgroup(*point);
}

inline void blst_p1_affine_compress(std::uint8_t *out, const blst_p1_affine* in)
{
    autohog::zk::mcl_backend::detail::ensure_initialized();
    const size_t expected = autohog::zk::mcl_backend::detail::g1_serialized_size();
    const size_t written = in->serialize(out, expected, mcl::IoSerialize);
    if (written != expected) {
        throw std::runtime_error("mcl G1 serialize size mismatch");
    }
}

inline void blst_p2_affine_compress(std::uint8_t *out, const blst_p2_affine* in)
{
    autohog::zk::mcl_backend::detail::ensure_initialized();
    const size_t expected = autohog::zk::mcl_backend::detail::g2_serialized_size();
    const size_t written = in->serialize(out, expected, mcl::IoSerialize);
    if (written != expected) {
        throw std::runtime_error("mcl G2 serialize size mismatch");
    }
}

inline BLST_ERROR blst_p1_uncompress(blst_p1_affine* out, const std::uint8_t* in)
{
    autohog::zk::mcl_backend::detail::ensure_initialized();
    blst_p1_affine point;
    const size_t expected = autohog::zk::mcl_backend::detail::g1_serialized_size();
    const size_t read = point.deserialize(in, expected, mcl::IoSerialize);
    if (read != expected) {
        return BLST_BAD_ENCODING;
    }
    point.normalize();
    if (!point.isZero() && !autohog::zk::mcl_backend::detail::g1_in_subgroup(point)) {
        return BLST_POINT_NOT_IN_GROUP;
    }
    *out = point;
    return BLST_SUCCESS;
}

inline BLST_ERROR blst_p2_uncompress(blst_p2_affine* out, const std::uint8_t* in)
{
    autohog::zk::mcl_backend::detail::ensure_initialized();
    blst_p2_affine point;
    const size_t expected = autohog::zk::mcl_backend::detail::g2_serialized_size();
    const size_t read = point.deserialize(in, expected, mcl::IoSerialize);
    if (read != expected) {
        return BLST_BAD_ENCODING;
    }
    point.normalize();
    if (!point.isZero() && !autohog::zk::mcl_backend::detail::g2_in_subgroup(point)) {
        return BLST_POINT_NOT_IN_GROUP;
    }
    *out = point;
    return BLST_SUCCESS;
}

inline size_t blst_p1s_mult_pippenger_scratch_sizeof(size_t)
{
    return sizeof(limb_t);
}

inline size_t blst_p2s_mult_pippenger_scratch_sizeof(size_t)
{
    return sizeof(limb_t);
}

inline void blst_p1s_mult_pippenger(blst_p1* out,
                                    const blst_p1_affine* const* points,
                                    size_t n,
                                    const unsigned char* const* scalar_ptrs,
                                    size_t nbits,
                                    limb_t*)
{
    autohog::zk::mcl_backend::detail::ensure_initialized();
    if (n == 0) {
        out->clear();
        return;
    }

    std::vector<blst_p1> bases;
    std::vector<blst_fr> scalars;
    autohog::zk::mcl_backend::detail::copy_points_from_ptrs(bases, points, n);
    autohog::zk::mcl_backend::detail::collect_scalars_from_ptrs(scalars, scalar_ptrs, n, nbits);
    blst_p1::mulVec(*out, bases.data(), scalars.data(), n);
}

inline void blst_p2s_mult_pippenger(blst_p2* out,
                                    const blst_p2_affine* const* points,
                                    size_t n,
                                    const unsigned char* const* scalar_ptrs,
                                    size_t nbits,
                                    limb_t*)
{
    autohog::zk::mcl_backend::detail::ensure_initialized();
    if (n == 0) {
        out->clear();
        return;
    }

    std::vector<blst_p2> bases;
    std::vector<blst_fr> scalars;
    autohog::zk::mcl_backend::detail::copy_points_from_ptrs(bases, points, n);
    autohog::zk::mcl_backend::detail::collect_scalars_from_ptrs(scalars, scalar_ptrs, n, nbits);
    blst_p2::mulVec(*out, bases.data(), scalars.data(), n);
}

inline void blst_p2s_add(blst_p2* out,
                         const blst_p2_affine* const* points,
                         size_t n)
{
    autohog::zk::mcl_backend::detail::ensure_initialized();
    out->clear();
    if (n == 0) {
        return;
    }
    *out = *points[0];
    for (size_t i = 1; i < n; ++i) {
        blst_p2::add(*out, *out, *points[i]);
    }
}

inline void blst_miller_loop(blst_fp12* out,
                             const blst_p2_affine* q,
                             const blst_p1_affine* p)
{
    autohog::zk::mcl_backend::detail::ensure_initialized();
    mcl::bn::millerLoop(*out, *p, *q);
}

inline void blst_final_exp(blst_fp12* out, const blst_fp12* in)
{
    autohog::zk::mcl_backend::detail::ensure_initialized();
    mcl::bn::finalExp(*out, *in);
}

inline void blst_fp12_mul(blst_fp12* out, const blst_fp12* a, const blst_fp12* b)
{
    blst_fp12::mul(*out, *a, *b);
}

inline bool blst_fp12_is_equal(const blst_fp12* a, const blst_fp12* b)
{
    return *a == *b;
}

inline blst_fr get_primitive_root_of_unity_bn_snark1(size_t log2_n) {
    autohog::zk::mcl_backend::detail::ensure_initialized();
    constexpr size_t kTwoAdicity = 28;
    static const mcl::Fr kRoot28 = [] {
        mcl::Fr g = 5, exp, root;
        // exp = (r-1) / 2^28
        exp.setStr("81540058820840996586704275553141814055101440848469862132140264610111", 10);
        mcl::Fr::pow(root, g, exp);
        return root;
    }();
    mcl::Fr root = kRoot28;
    for (size_t i = 0; i < kTwoAdicity - log2_n; i++) root *= root;
    return root;
}
