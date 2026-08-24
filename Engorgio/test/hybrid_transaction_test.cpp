#include "ENGORGIO/comp.h"
#include "ENGORGIO/utils.h"
#include <iostream>
#include <chrono>
#include <vector>

using namespace openfhe;
using namespace lbcrypto;

 






double Execute_Parallel_Transaction_64bit(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keyPair, int length) {
    uint32_t slots = cc->GetRingDimension() / 2; 
    int Bg = 7;
    double precision = (1 << (8 - 1)) - 1; 
    uint32_t polyDegree = 119;

    
    std::vector<double> A_balance_vec(slots, 0.0);
    std::vector<double> B_balance_vec(slots, 0.0);
    std::vector<double> amount_vec(slots, 0.0);

    
    
    A_balance_vec[0] = 5000.0; B_balance_vec[0] = 2000.0; amount_vec[0] = 1500.0;
    
    A_balance_vec[1] = 1000.0; B_balance_vec[1] = 2000.0; amount_vec[1] = 1500.0;
    
    A_balance_vec[2] = 3000.0; B_balance_vec[2] = 1000.0; amount_vec[2] = 500.0;

    
    auto A_bal_quant = quantization(A_balance_vec, 16, Bg); 
    auto amount_quant = quantization(amount_vec, 16, Bg);

    
    Ciphertext<DCRTPoly> c_A_balance = cc->Encrypt(keyPair.publicKey, cc->MakeCKKSPackedPlaintext(A_balance_vec));
    Ciphertext<DCRTPoly> c_B_balance = cc->Encrypt(keyPair.publicKey, cc->MakeCKKSPackedPlaintext(B_balance_vec));
    Ciphertext<DCRTPoly> c_amount = cc->Encrypt(keyPair.publicKey, cc->MakeCKKSPackedPlaintext(amount_vec));

    std::vector<Ciphertext<DCRTPoly>> c_A_quant_vec, c_amt_quant_vec;
    for (auto& q : A_bal_quant) c_A_quant_vec.push_back(cc->Encrypt(keyPair.publicKey, cc->MakeCKKSPackedPlaintext(q)));
    for (auto& q : amount_quant) c_amt_quant_vec.push_back(cc->Encrypt(keyPair.publicKey, cc->MakeCKKSPackedPlaintext(q)));

    
    auto start = std::chrono::high_resolution_clock::now();

    
    Ciphertext<DCRTPoly> sign;
    comp_greater_than_modular(c_A_quant_vec, c_amt_quant_vec, precision, polyDegree, sign, keyPair.secretKey);

    
    auto amount_true = cc->EvalMult(c_amount, sign);

    
    auto A_result = cc->EvalSub(c_A_balance, amount_true); 
    auto B_result = cc->EvalAdd(c_B_balance, amount_true); 

    auto end = std::chrono::high_resolution_clock::now();
    

    double duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    
    Plaintext resA, resB;
    cc->Decrypt(keyPair.secretKey, A_result, &resA);
    cc->Decrypt(keyPair.secretKey, B_result, &resB);
    
    std::cout << "\n[Parallel Execution Result for length=" << length << "]" << std::endl;
    for(int i = 0; i < length; i++) {
        std::cout << "Transaction " << i << ": A_Final=" << std::round(resA->GetCKKSPackedValue()[i].real()) 
                  << ", B_Final=" << std::round(resB->GetCKKSPackedValue()[i].real()) << std::endl;
    }

    return duration;
}

double Execute_Parallel_Transaction_Original_Logic(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keyPair, int length) {
    uint32_t slots = cc->GetRingDimension() / 2; 
    int Bg = 7;
    double precision = (1 << (8 - 1)) - 1; 
    uint32_t polyDegree = 119;

    
    std::vector<double> A_balance_vec(slots, 0.0);
    std::vector<double> amount_vec(slots, 0.0);

    
    A_balance_vec[0] = 5000.0; amount_vec[0] = 1500.0; 
    A_balance_vec[1] = 1000.0; amount_vec[1] = 1500.0; 
    A_balance_vec[2] = 3000.0; amount_vec[2] = 500.0;  

    
    auto A_bal_quant = quantization(A_balance_vec, 16, Bg); 
    auto amount_quant = quantization(amount_vec, 16, Bg);

    
    Ciphertext<DCRTPoly> c_A_balance = cc->Encrypt(keyPair.publicKey, cc->MakeCKKSPackedPlaintext(A_balance_vec));
    Ciphertext<DCRTPoly> c_amount = cc->Encrypt(keyPair.publicKey, cc->MakeCKKSPackedPlaintext(amount_vec));

    std::vector<Ciphertext<DCRTPoly>> c_A_quant_vec, c_amt_quant_vec;
    for (auto& q : A_bal_quant) c_A_quant_vec.push_back(cc->Encrypt(keyPair.publicKey, cc->MakeCKKSPackedPlaintext(q)));
    for (auto& q : amount_quant) c_amt_quant_vec.push_back(cc->Encrypt(keyPair.publicKey, cc->MakeCKKSPackedPlaintext(q)));

    
    auto start = std::chrono::high_resolution_clock::now();

    
    Ciphertext<DCRTPoly> sign;
    comp_greater_than_modular(c_A_quant_vec, c_amt_quant_vec, precision, polyDegree, sign, keyPair.secretKey);

    
    auto a = cc->EvalSub(c_A_balance, c_amount); 
    
    
    auto b = c_A_balance;

    
    
    auto term1 = cc->EvalMult(sign, a); 
    
    
    auto one_minus_sign = cc->EvalSub(1.0, sign);
    auto term2 = cc->EvalMult(one_minus_sign, b);
    
    
    auto final_res = cc->EvalAdd(term1, term2);

    auto end = std::chrono::high_resolution_clock::now();
    

    double duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    
    Plaintext res;
    cc->Decrypt(keyPair.secretKey, final_res, &res);
    std::cout << "\n[Logic: sign*a + (1-sign)*b | length=" << length << "]" << std::endl;
    for(int i = 0; i < length; i++) {
        std::cout << "Transaction " << i << " Result: " << std::round(res->GetCKKSPackedValue()[i].real()) << std::endl;
    }

    return duration;
}

int main() {
    
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecurityLevel(HEStd_128_classic);
    parameters.SetScalingModSize(50);
    parameters.SetFirstModSize(60);
    parameters.SetMultiplicativeDepth(30); 
    
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE); cc->Enable(ADVANCEDSHE);

    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    
    uint32_t slots = cc->GetRingDimension() / 2;
    std::vector<int> rotstep;
    for (int i = 1; i < slots; i *= 2) { rotstep.push_back(i); rotstep.push_back(-i); }
    cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);

    
    int length = 3;
    std::cout << "Starting Parallel Transaction Computation (length = " << length << ")..." << std::endl;
    
    // double compute_time = Execute_Parallel_Transaction_64bit(cc, keyPair, length);
    double compute_time = Execute_Parallel_Transaction_Original_Logic(cc, keyPair, length);

    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Total Parallel Compute Time: " << compute_time << " ms" << std::endl;
    std::cout << "Amortized Time per Transaction: " << compute_time / length << " ms" << std::endl;

    return 0;
}