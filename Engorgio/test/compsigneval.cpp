#include "openfhe.h"
#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/ckksrns/ckksrns-ser.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <functional>
#include <complex>

using namespace lbcrypto;

std::vector<double> coeff1 = {2.38502034e+01, -6.43644764e+02, 6.39413732e+03, -2.92814270e+04, 6.99651923e+04, -9.02320081e+04, 5.95218847e+04, -1.57478915e+04};
std::vector<double> coeff2 = {1.03676385e+01, -1.56686831e+02, 1.20109138e+03, -4.87084195e+03, 1.17135632e+04, -1.80263166e+04, 1.86247266e+04, -1.32876049e+04, 6.62756768e+03, -2.30508937e+03, 5.47705722e+02, -8.47355934e+01, 7.69135049e+00, -3.10721125e-01};
std::vector<double> coeff3 = {2.79543701e+00, -4.57690533e+00, 4.12196169e+00, 0.00000000e+00, -2.36649616e+00, 0.00000000e+00, 1.66909946e+00, 0.00000000e+00, -9.26141831e-01, 0.00000000e+00, 3.67665069e-01, 0.00000000e+00, -1.16917842e-01, 3.22977367e-02};

Ciphertext<DCRTPoly> MultByInteger(const ConstCiphertext<DCRTPoly>& ciphertext, const int64_t constant) {
    auto result = ciphertext->Clone();
    for (auto& c : result->GetElements())
        c *= static_cast<typename DCRTPoly::Integer>(constant);
    return result;
}

template <typename T>
inline uint64_t CeilLog2(T x) {
    return static_cast<uint64_t>(std::ceil(std::log2(x)));
}

void EvalPower(std::vector<double> coefficients, std::vector<Ciphertext<lbcrypto::DCRTPoly>>& power_basis,
               Ciphertext<lbcrypto::DCRTPoly>& result) {
    auto cc = power_basis[0]->GetCryptoContext();
    if (coefficients.size() == 1) {
        if (std::fabs(std::round(coefficients[0] * pow(2, 50))) > 1.) {
            result = power_basis[0];
            result = cc->EvalMult(result, coefficients[0]);
            return;
        }
        else {
            return;
        }
    }
    std::vector<double> quotient, remainder;
    uint64_t degree = 2 * coefficients.size() - 1;
    uint64_t m      = CeilLog2(degree + 1);
    remainder.resize((1 << (m - 1)) / 2);
    quotient.resize(coefficients.size() - remainder.size());

    for (size_t i = 0; i < remainder.size(); i++) {
        remainder[i] = coefficients[i];
    }
    for (size_t i = 0; i < quotient.size(); i++) {
        quotient[i] = coefficients[i + remainder.size()];
    }
    Ciphertext<lbcrypto::DCRTPoly> cipher_quotient, cipher_remainder;
    EvalPower(quotient, power_basis, cipher_quotient);
    EvalPower(remainder, power_basis, cipher_remainder);
    result = cc->EvalMult(cipher_quotient, power_basis[m - 1]);
    result = cc->EvalAdd(result, cipher_remainder);
}

void poly_evaluate_power(Ciphertext<DCRTPoly>& result, Ciphertext<DCRTPoly>& x,
                         std::vector<double>& coefficients) {
    auto cc         = x->GetCryptoContext();
    uint64_t degree = coefficients.size() * 2 - 1;
    uint64_t m      = CeilLog2(degree + 1);

    std::vector<Ciphertext<DCRTPoly>> power_basis(m);
    power_basis[0] = x;
    for (size_t i = 1; i < m; i++) {
        power_basis[i] = cc->EvalMult(power_basis[i - 1], power_basis[i - 1]);
    }

    EvalPower(coefficients, power_basis, result);
}

Ciphertext<DCRTPoly> compwithzero(Ciphertext<DCRTPoly>& x){
    auto cc = x->GetCryptoContext();
    Ciphertext<DCRTPoly> result, x1, x2;
    poly_evaluate_power(x1, x, coeff1);
    poly_evaluate_power(x2, x1, coeff2);
    for (auto &x : coeff3) {
    x /= 2.0;
    }
    poly_evaluate_power(result, x2, coeff3);
    result = cc->EvalAdd(result, 0.5);
    return result;
}

Ciphertext<DCRTPoly> eqzero(Ciphertext<DCRTPoly>& x){
    auto cc = x->GetCryptoContext();
    Ciphertext<DCRTPoly> result, x1, x2;
    
    
    poly_evaluate_power(x1, x, coeff1);
    poly_evaluate_power(x2, x1, coeff2);
    for (auto &c : coeff3) {
        c /= 2.0;
    }
    poly_evaluate_power(result, x2, coeff3);
    result = cc->EvalAdd(result, 0.5);
    
    
    Ciphertext<DCRTPoly> one_minus_sign = cc->EvalSub(1.0, result);
    Ciphertext<DCRTPoly> sign_times_one_minus_sign = cc->EvalMult(result, one_minus_sign);
    result = cc->EvalMult(sign_times_one_minus_sign, 4.0);
    
    return result;
}

int main() {

    CCParams<CryptoContextCKKSRNS> parameters;

    size_t BatchSize = 16;

    parameters.SetMultiplicativeDepth(16);
    parameters.SetScalingModSize(59);
    parameters.SetSecurityLevel(HEStd_128_classic);
    parameters.SetRingDim(65536);
    parameters.SetScalingTechnique(FLEXIBLEAUTO);
    parameters.SetFirstModSize(60);       
    parameters.SetBatchSize(BatchSize);
    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);
    auto keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeyGen(keyPair.secretKey);

    std::vector<double> data(BatchSize, 0);
    for(size_t i = 0; i < BatchSize; ++i){
        data[i] = -10+0.625+1.25*i;
        if(i == 0)
            data[i] = 0;
    }
    std::cout << "data: "<<data<<std::endl;
    Plaintext ptxt_data = cryptoContext->MakeCKKSPackedPlaintext(data);
    Ciphertext<DCRTPoly> ciphertext_data = cryptoContext->Encrypt(keyPair.publicKey, ptxt_data);
    Ciphertext<DCRTPoly> result;
    

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::high_resolution_clock::now();

    // result = compwithzero(ciphertext_data);
    result = eqzero(ciphertext_data);

    end = std::chrono::high_resolution_clock::now();
    double filtering_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "filtering time: " << filtering_time << " ms" << std::endl;
    std::cout << "amortized: " << filtering_time / 65536 << " ms" << std::endl;

    Plaintext result_ptxt;
    cryptoContext->Decrypt(keyPair.secretKey, result, &result_ptxt);
    result_ptxt->SetLength(BatchSize);
    auto decvec = result_ptxt->GetRealPackedValue();
    for(size_t j = 0; j < BatchSize; ++j){
        std::cout << decvec[j] << " ";
    } 
    std::cout << std::endl;

    return 0;
}