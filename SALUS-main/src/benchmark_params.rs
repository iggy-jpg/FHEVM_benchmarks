//! SALUS-compatible cryptographic parameters shared by all benchmark binaries.
//!
//! These are the original parameter sets used by SALUS for its scalar (CMUX0)
//! and grouped (CMUX1) circuit-bootstrap routes. Keeping both parameter sets
//! unchanged is required for correct selector bootstrapping in hybrid mode.

use refined_tfhe_lhe::int_lhe_instance::{SALUS_CMUX0, SALUS_CMUX1};
use refined_tfhe_lhe::int_lhe_params::IntLheParam;

/// Parameters for scalar/single-output circuit bootstrapping.
pub fn benchmark_cmux0_param() -> IntLheParam<u64> {
    *SALUS_CMUX0
}

/// Parameters for grouped/multi-output circuit bootstrapping.
pub fn benchmark_cmux1_param() -> IntLheParam<u64> {
    *SALUS_CMUX1
}

/// Print both effective parameter sets with every benchmark result.
pub fn print_benchmark_parameters() {
    let cmux0 = benchmark_cmux0_param();
    let cmux1 = benchmark_cmux1_param();

    println!("parameters:");
    println!("profile: salus-compatible");
    print_route("cmux0", &cmux0);
    print_route("cmux1", &cmux1);
    println!("ciphertext_modulus: 2^64 (native)");
    println!();
}

fn print_route(name: &str, params: &IntLheParam<u64>) {
    println!("{name}_lwe_dimension: {}", params.lwe_dimension().0);
    println!(
        "{name}_lwe_noise_standard_dev: {:.17e}",
        params.lwe_modular_std_dev().0
    );
    println!("{name}_glwe_dimension: {}", params.glwe_dimension().0);
    println!(
        "{name}_glwe_noise_standard_dev: {:.17e}",
        params.glwe_modular_std_dev().0
    );
    println!("{name}_polynomial_size: {}", params.polynomial_size().0);
    println!("{name}_pbs_base_log: {}", params.pbs_base_log().0);
    println!("{name}_pbs_level: {}", params.pbs_level().0);
    println!("{name}_ks_base_log: {}", params.ks_base_log().0);
    println!("{name}_ks_level: {}", params.ks_level().0);
    println!("{name}_auto_base_log: {}", params.auto_base_log().0);
    println!("{name}_auto_level: {}", params.auto_level().0);
    println!("{name}_ss_base_log: {}", params.ss_base_log().0);
    println!("{name}_ss_level: {}", params.ss_level().0);
    println!("{name}_cbs_base_log: {}", params.cbs_base_log().0);
    println!("{name}_cbs_level: {}", params.cbs_level().0);
    println!("{name}_log_lut_count: {}", params.log_lut_count().0);
}

#[cfg(test)]
mod tests {
    use super::*;
    use tfhe::core_crypto::prelude::{
        DecompositionBaseLog, DecompositionLevelCount, GlweDimension, LweDimension,
        PolynomialSize,
    };

    #[test]
    fn cmux0_matches_original_salus_parameters() {
        let params = benchmark_cmux0_param();
        assert_eq!(params.lwe_dimension(), LweDimension(571));
        assert_eq!(params.glwe_dimension(), GlweDimension(1));
        assert_eq!(params.polynomial_size(), PolynomialSize(2048));
        assert_eq!(params.pbs_base_log(), DecompositionBaseLog(15));
        assert_eq!(params.pbs_level(), DecompositionLevelCount(2));
        assert_eq!(params.ks_base_log(), DecompositionBaseLog(2));
        assert_eq!(params.ks_level(), DecompositionLevelCount(7));
    }

    #[test]
    fn cmux1_matches_original_salus_parameters() {
        let params = benchmark_cmux1_param();
        assert_eq!(params.lwe_dimension(), LweDimension(643));
        assert_eq!(params.glwe_dimension(), GlweDimension(1));
        assert_eq!(params.polynomial_size(), PolynomialSize(2048));
        assert_eq!(params.pbs_base_log(), DecompositionBaseLog(12));
        assert_eq!(params.pbs_level(), DecompositionLevelCount(3));
        assert_eq!(params.ks_base_log(), DecompositionBaseLog(1));
        assert_eq!(params.ks_level(), DecompositionLevelCount(14));
    }

    #[test]
    fn hybrid_routes_can_share_the_glwe_secret() {
        let cmux0 = benchmark_cmux0_param();
        let cmux1 = benchmark_cmux1_param();
        assert_eq!(cmux0.glwe_dimension(), cmux1.glwe_dimension());
        assert_eq!(cmux0.polynomial_size(), cmux1.polynomial_size());
    }
}
