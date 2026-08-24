#pragma once
#include <cstddef>
#include <cmath>
#include <limits>
#include "tfhe++.hpp"

// ============================================================================


// ============================================================================

struct DriftMsTestParamLvl0 {
    
    using P     = TFHEpp::lvl0param;
    using Torus = typename P::T;

    
    static constexpr std::size_t lwe_dim =
        static_cast<std::size_t>(P::k) * P::n;

    
    static constexpr int log_mod = 12;

    
    static constexpr std::size_t zeros_count = 2000;

    
    static constexpr double r_sigma = 13.179852282053789;
    

    
    //
    // Var(edrift) ≈ ((n+2) q^2) / (96 N^2) + (n-4)/48
    
    //   q = 2^BITS          （Torus32 -> BITS=32）
    //   N = 2^(log_mod-1) = 2048
    static double estimated_drift_variance()
    {
        const double n   = static_cast<double>(P::n);
        const double q   = std::ldexp(1.0,
                                      std::numeric_limits<Torus>::digits); // 2^BITS
        const double N   = std::ldexp(1.0, log_mod - 1);                    // 2^11 = 2048
        const double term1 = (n + 2.0) * q * q / (96.0 * N * N);
        const double term2 = (n - 4.0) / 48.0;
        return term1 + term2;
    }

    
    //   sigma ≈ sqrt(Var)，T ≈ 1.3 * r_sigma * sigma
    static double bound()
    {
        const double var   = estimated_drift_variance();
        const double sigma = std::sqrt(var);
        return 1.3 * r_sigma * sigma; 
    }

    
    
    //   input_var ≈ (α * 2^BITS)^2
    static double input_variance()
    {
        constexpr int BITS = std::numeric_limits<Torus>::digits; // 32
        const double sigma_int = P::α * std::ldexp(1.0, BITS);
        return sigma_int * sigma_int;
    }

};

// ============================================================================

// ============================================================================
struct DriftMsTestParamLvl2 {
    using P     = TFHEpp::lvl2param;
    using Torus = typename P::T; // usually uint64_t

    static constexpr std::size_t lwe_dim =
        static_cast<std::size_t>(P::k) * static_cast<std::size_t>(P::n);

    // MUST match the MS/BR indexing modulus you care about.
    // If your BR uses 2N=4096 => log_mod=12
    static constexpr int log_mod = 11;

    // lvl2 zeros are expensive: start smaller than 2000
    static constexpr std::size_t zeros_count = 3000;

    // Tail radius (kept consistent with your lvl0 experiments)
    static constexpr double r_sigma = 13.179852282053789;

    // bound is disabled (always-min policy)
    static double bound()
    {
        return std::numeric_limits<double>::infinity();
    }

    // Convert normalized alpha to torus integer variance scale
    static double input_variance()
    {
        constexpr int BITS = std::numeric_limits<Torus>::digits; // 64
        const double sigma_int = P::α * std::ldexp(1.0, BITS);
        return sigma_int * sigma_int;
    }
};


// ============================================================================



// ============================================================================

// struct DriftMsTestParamLvl0 {

//     using P     = lvl0param;
//     using Torus = typename P::T;


//     static constexpr std::size_t lwe_dim =
//         static_cast<std::size_t>(P::k) * P::n;


//     static constexpr int log_mod = 12;


//     static constexpr std::size_t zeros_count = 1449;



//     static constexpr double r_sigma = 13.179852282053789;



//     static constexpr double bound_value = 288230376151711744.0;


//     //

//     //   let modulus = ciphertext_modulus.raw_modulus_float();
//     //   let input_variance = input_variance.get_modular_variance(modulus);
//     //



//     static constexpr double input_variance_modular = 2.63039184094559e-7;


//     static double estimated_drift_variance()
//     {

//         return 0.0;
//     }


//     static double bound()
//     {
//         return bound_value;
//     }


//     static double input_variance()
//     {
//         return input_variance_modular;
//     }
// };
