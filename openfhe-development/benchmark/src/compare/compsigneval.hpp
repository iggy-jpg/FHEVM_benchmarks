#pragma once

// Integrated from testfiles/Compare/compsigneval.cpp.  The original composite
// minimax sign evaluator is exposed as a reusable OpenFHE helper.
#include "openfhe.h"

#include <cmath>
#include <cstdint>
#include <vector>

namespace compare_minimax {

using lbcrypto::Ciphertext;
using lbcrypto::DCRTPoly;

inline std::uint64_t CeilLog2(std::uint64_t value)
{
    return static_cast<std::uint64_t>(std::ceil(std::log2(value)));
}

inline void EvalOddPower(
    const std::vector<double>& coefficients,
    std::vector<Ciphertext<DCRTPoly>>& powerBasis,
    Ciphertext<DCRTPoly>& result)
{
    const auto cc = powerBasis[0]->GetCryptoContext();
    if (coefficients.size() == 1) {
        if (std::fabs(std::round(coefficients[0] * std::pow(2.0, 50))) > 1.0)
            result = cc->EvalMult(powerBasis[0], coefficients[0]);
        return;
    }
    const std::uint64_t degree = 2 * coefficients.size() - 1;
    const std::uint64_t m = CeilLog2(degree + 1);
    const std::size_t remainderSize = (std::uint64_t{1} << (m - 1)) / 2;
    std::vector<double> remainder(coefficients.begin(), coefficients.begin() + remainderSize);
    std::vector<double> quotient(coefficients.begin() + remainderSize, coefficients.end());
    Ciphertext<DCRTPoly> quotientCt, remainderCt;
    EvalOddPower(quotient, powerBasis, quotientCt);
    EvalOddPower(remainder, powerBasis, remainderCt);
    result = cc->EvalMult(quotientCt, powerBasis[m - 1]);
    result = cc->EvalAdd(result, remainderCt);
}

inline Ciphertext<DCRTPoly> EvalOddPolynomial(
    const Ciphertext<DCRTPoly>& input,
    const std::vector<double>& coefficients)
{
    const auto cc = input->GetCryptoContext();
    const std::uint64_t degree = 2 * coefficients.size() - 1;
    const std::uint64_t m = CeilLog2(degree + 1);
    std::vector<Ciphertext<DCRTPoly>> powerBasis(m);
    powerBasis[0] = input;
    for (std::size_t i = 1; i < m; ++i)
        powerBasis[i] = cc->EvalMult(powerBasis[i - 1], powerBasis[i - 1]);
    Ciphertext<DCRTPoly> result;
    EvalOddPower(coefficients, powerBasis, result);
    return result;
}

inline Ciphertext<DCRTPoly> Step(const Ciphertext<DCRTPoly>& input)
{
    static const std::vector<double> coeff1 = {
        1.39100429e+00, -1.32778222e-01, 4.38520070e-03,
        -5.61608992e-05, 2.43818133e-07};
    static const std::vector<double> coeff2 = {
        6.40080258, -33.59524601, 83.69219412, -102.97210232,
        68.00962186, -24.59588521, 4.58765824, -0.34499159};
    static const std::vector<double> coeff3 = {
        3.20819335 / 2.0, -7.73696017 / 2.0, 14.22957007 / 2.0,
        -17.11288053 / 2.0, 13.29608953 / 2.0, -6.44914042 / 2.0,
        1.77809638 / 2.0, -0.21296897 / 2.0};
    const auto cc = input->GetCryptoContext();
    auto first = EvalOddPolynomial(input, coeff1);
    auto second = EvalOddPolynomial(first, coeff2);
    auto result = EvalOddPolynomial(second, coeff3);
    return cc->EvalAdd(result, 0.5);
}

// Inputs are integer-valued CKKS ciphertexts.  Dividing by 128 maps the
// half-integer decision boundary into the approximation interval used by the
// coefficients from Compare.
inline Ciphertext<DCRTPoly> GreaterThan(
    const Ciphertext<DCRTPoly>& lhs, const Ciphertext<DCRTPoly>& rhs)
{
    const auto cc = lhs->GetCryptoContext();
    auto shifted = cc->EvalSub(cc->EvalSub(lhs, rhs), 0.5);
    return Step(cc->EvalMult(shifted, 1.0 / 128.0));
}

inline Ciphertext<DCRTPoly> Equal(
    const Ciphertext<DCRTPoly>& lhs, const Ciphertext<DCRTPoly>& rhs)
{
    const auto cc = lhs->GetCryptoContext();
    const auto delta = cc->EvalSub(lhs, rhs);
    const auto squared = cc->EvalMult(delta, delta);
    // For integer inputs, 0.5 - (lhs-rhs)^2 is positive exactly when the
    // values are equal and is at most -0.5 otherwise.
    return Step(cc->EvalAdd(cc->EvalMult(squared, -1.0), 0.5));
}

} // namespace compare_minimax
