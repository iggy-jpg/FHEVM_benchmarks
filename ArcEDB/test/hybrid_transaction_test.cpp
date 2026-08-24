#include "ARCEDB/comparison/rgsw_ciphertext.h"
#include "ARCEDB/comparison/comparable.h"
#include "ARCEDB/conversion/packlwes.h"
#include "ARCEDB/comparison/batch_bootstrap.h"
#include "ARCEDB/utils/serialize.h"
#include "ARCEDB/conversion/repack.h"
#include <iomanip>
#include <vector>
#include <chrono>

using namespace arcedb;
using namespace seal;

void Execute_Batch_Transaction() {
    std::cout << ">>>> Starting Batch Transaction (3-way Parallel) <<<<" << std::endl;

    
    
    std::vector<uint64_t> b_vals = {1000, 500, 100};
    std::vector<uint64_t> a_vals = {200, 400, 200};
    std::vector<uint64_t> f_vals = {50, 100, 10};
    std::vector<uint64_t> thresholds = {250, 500, 210};

    
    using P = Lvl1;
    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplacebkfft<Lvl02>(sk);

    EncryptionParameters parms(scheme_type::ckks);
    size_t n = 65536;
    parms.set_poly_modulus_degree(n);
    parms.set_coeff_modulus(CoeffModulus::Create(n, {59, 42, 42, 42, 42, 42, 42, 42, 42, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 59}));
    
    SEALContext context(parms, true, seal::sec_level_type::none);
    KeyGenerator keygen(context);
    SecretKey seal_sk = keygen.secret_key();
    PublicKey public_key; keygen.create_public_key(public_key);
    RelinKeys rlk; keygen.create_relin_keys(rlk);
    GaloisKeys gk; keygen.create_galois_keys(gk);
    
    Encryptor encryptor(context, public_key);
    encryptor.set_secret_key(seal_sk);
    Evaluator evaluator(context);
    Decryptor decryptor(context, seal_sk);
    CKKSEncoder encoder(context);

    auto start_time = std::chrono::high_resolution_clock::now();

    
    
    std::vector<TLWELvl1> lt_results(3);
    for(int i=0; i<3; i++) {
        ComparableLvl1 c_b;
        exponent_encrypt<P>(b_vals[i], 64, c_b, sk);
        std::vector<TRGSWLvl1> c_t(2);
        exponent_encrypt_rgsw<P>(thresholds[i], 64, c_t, sk, true);
        less_than_tfhepp(c_b, c_t, c_b.size(), lt_results[i], ek, sk);
    }

    
    LTPreKey pre_key;
    LWEsToRLWEKeyGen(pre_key, std::pow(2., 45), seal_sk, sk, P::n, encoder, encryptor, context);
    
    Ciphertext lt_sign_batch;
    
    LWEsToRLWE(lt_sign_batch, lt_results, pre_key, std::pow(2., 29), std::pow(2., 32), std::pow(2., 13), encoder, gk, rlk, evaluator, context);
    HomRound(lt_sign_batch, lt_sign_batch.scale(), encoder, rlk, evaluator, decryptor, context);

    
    auto current_parms_id = lt_sign_batch.parms_id();
    double current_scale = lt_sign_batch.scale();

    
    Plaintext p_one;
    encoder.encode(1.0, current_parms_id, current_scale, p_one); 

    Ciphertext sign_batch;
    evaluator.sub_plain(lt_sign_batch, p_one, sign_batch);
    evaluator.negate_inplace(sign_batch); // 1 if balance >= threshold else 0

    
    Ciphertext inv_sign_batch;
    evaluator.sub_plain(sign_batch, p_one, inv_sign_batch);
    evaluator.negate_inplace(inv_sign_batch); // 1 if balance < threshold else 0

    
    std::vector<double> res_if_enough(encoder.slot_count(), 0.0);
    std::vector<double> res_if_not_enough(encoder.slot_count(), 0.0);
    for(int i=0; i<3; i++) {
        res_if_enough[i] = (double)(b_vals[i] - thresholds[i]);
        res_if_not_enough[i] = (double)b_vals[i];
    }

    Plaintext p_res_enough, p_res_not_enough;
    encoder.encode(res_if_enough, sign_batch.parms_id(), current_scale, p_res_enough);
    encoder.encode(res_if_not_enough, sign_batch.parms_id(), current_scale, p_res_not_enough);

    // term1 = sign * (b - threshold)
    Ciphertext term1;
    evaluator.multiply_plain(sign_batch, p_res_enough, term1);
    evaluator.rescale_to_next_inplace(term1);

    // term2 = (1 - sign) * b
    Ciphertext term2;
    evaluator.multiply_plain(inv_sign_batch, p_res_not_enough, term2);
    evaluator.rescale_to_next_inplace(term2);

    
    
    double final_sync_scale = term1.scale();
    term2.scale() = final_sync_scale;

    
    if (term1.parms_id() != term2.parms_id()) {
        evaluator.mod_switch_to_inplace(term1, term2.parms_id());
    }

    
    Ciphertext res_batch;
    evaluator.add(term1, term2, res_batch);

    auto end_time = std::chrono::high_resolution_clock::now();

    
    Plaintext p_final;
    std::vector<double> v_final;
    decryptor.decrypt(res_batch, p_final);
    encoder.decode(p_final, v_final);

    double elapsed = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    std::cout << "Batch Latency: " << elapsed << " ms" << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    for(int i=0; i<3; i++) {
        uint64_t expected = (b_vals[i] >= thresholds[i]) ? (b_vals[i] - thresholds[i]) : b_vals[i];
        std::cout << "Tx #" << i+1 << " | Result: " << std::fixed << std::setprecision(0) << v_final[i] 
                  << " | Expected: " << expected << (std::round(v_final[i]) == expected ? " [PASS]" : " [FAIL]") << std::endl;
    }
}



int main() {
    Execute_Batch_Transaction();
    return 0;
}