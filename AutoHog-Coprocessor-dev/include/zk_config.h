#pragma once

#include <cstddef>
#include <cstdint>

// =======================
// Backend selection check
// =======================

#if defined(AUTOHOG_ZK_BACKEND_mcl) && defined(AUTOHOG_ZK_BACKEND_blst)
#error "Select only one ZK backend: mcl or blst"
#endif

#if !defined(AUTOHOG_ZK_BACKEND_mcl) && !defined(AUTOHOG_ZK_BACKEND_blst)
#define AUTOHOG_ZK_BACKEND_blst 1
#endif

// =======================
// Curve selection check
// =======================

#if defined(AUTOHOG_ZK_CURVE_bls12_381) && defined(AUTOHOG_ZK_CURVE_bn_snark1)
#error "Select only one ZK curve: bls12_381 or bn_snark1"
#endif

#if !defined(AUTOHOG_ZK_CURVE_bls12_381) && !defined(AUTOHOG_ZK_CURVE_bn_snark1)
#define AUTOHOG_ZK_CURVE_bls12_381 1
#endif

// BLST restriction
#if defined(AUTOHOG_ZK_BACKEND_blst) && defined(AUTOHOG_ZK_CURVE_bn_snark1)
#error                                                                         \
    "BLST backend does not support bn_snark1 curve. Use bls12_381 or switch to mcl backend."
#endif

// =======================
// Derived config values
// =======================

namespace autohog::zk::config {

enum class BackendId : std::uint8_t {
  blst = 1,
  mcl = 2,
};

enum class CurveId : std::uint8_t {
  bls12_381 = 1,
  bn_snark1 = 2,
};

#if defined(AUTOHOG_ZK_BACKEND_mcl)
inline constexpr BackendId backend_id = BackendId::mcl;
inline constexpr const char *backend_name = "mcl";
#elif defined(AUTOHOG_ZK_BACKEND_blst)
inline constexpr BackendId backend_id = BackendId::blst;
inline constexpr const char *backend_name = "blst";
#else
#error "Unknown ZK backend selection. Define AUTOHOG_ZK_BACKEND_[mcl|blst]"
#endif

#if defined(AUTOHOG_ZK_CURVE_bn_snark1)
inline constexpr CurveId curve_id = CurveId::bn_snark1;
inline constexpr const char *curve_name = "bn_snark1";
inline constexpr std::size_t g1_serialized_size = 32;
inline constexpr std::size_t g2_serialized_size = 64;
inline constexpr std::size_t fr_two_adicity = 28;
#elif defined(AUTOHOG_ZK_CURVE_bls12_381)
inline constexpr CurveId curve_id = CurveId::bls12_381;
inline constexpr const char *curve_name = "bls12_381";
inline constexpr std::size_t g1_serialized_size = 48;
inline constexpr std::size_t g2_serialized_size = 96;
inline constexpr std::size_t fr_two_adicity = 32;
#else
#error                                                                         \
    "Unknown ZK curve selection. Define AUTOHOG_ZK_CURVE_[bls12_381|bn_snark1]"
#endif

inline constexpr std::size_t proof_serialized_size =
    g2_serialized_size + 2 * g1_serialized_size;

} // namespace autohog::zk::config