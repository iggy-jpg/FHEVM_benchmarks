//==================================================================================
// BSD 2-Clause License
//
// This file has been modified from the original version.
// Changes made by Jules Dumezy at CEA-List in 2025.
//
// Copyright (c) 2025, CEA-List
//
// Author TPOC: jules.dumezy@cea.fr
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*
  Control for encryption operations
 */

#include "cryptocontext.h"
#include "omp.h"
#include "key/privatekey.h"
#include "key/publickey.h"
#include "math/chebyshev.h"
#include "math/hermite.h"
#include "schemerns/rns-scheme.h"
#include "scheme/ckksrns/ckksrns-cryptoparameters.h"
#include "utils/exception.h"

namespace lbcrypto {

template <typename Element>
std::map<std::string, std::vector<EvalKey<Element>>> CryptoContextImpl<Element>::s_evalMultKeyMap{};
template <typename Element>
std::map<std::string, std::shared_ptr<std::map<usint, EvalKey<Element>>>>
    CryptoContextImpl<Element>::s_evalAutomorphismKeyMap{};

template <typename Element>
void CryptoContextImpl<Element>::SetKSTechniqueInScheme() {
    // check if the scheme is an RNS scheme
    auto schemeRNSPtr = dynamic_cast<SchemeRNS*>(&(*scheme));
    if (schemeRNSPtr != nullptr) {
        // check if the parameter object is RNS-based
        auto elPtr = dynamic_cast<const CryptoParametersRNS*>(&(*params));
        if (elPtr != nullptr) {
            schemeRNSPtr->SetKeySwitchingTechnique(elPtr->GetKeySwitchTechnique());
            return;
        }
        OPENFHE_THROW("Cannot set KeySwitchingTechnique as the parameter object is not RNS-based");
    }
}

/////////////////////////////////////////
// SHE MULTIPLICATION
/////////////////////////////////////////
template <typename Element>
void CryptoContextImpl<Element>::EvalMultKeyGen(const PrivateKey<Element> key) {
    ValidateKey(key);

    if (CryptoContextImpl<Element>::s_evalMultKeyMap.find(key->GetKeyTag()) ==
        CryptoContextImpl<Element>::s_evalMultKeyMap.end()) {
        // the key is not found in the map, so the key has to be generated
        EvalKey<Element> k                                           = GetScheme()->EvalMultKeyGen(key);
        CryptoContextImpl<Element>::s_evalMultKeyMap[k->GetKeyTag()] = {k};
    }
}

template <typename Element>
void CryptoContextImpl<Element>::EvalMultKeysGen(const PrivateKey<Element> key) {
    ValidateKey(key);

    if (CryptoContextImpl<Element>::s_evalMultKeyMap.find(key->GetKeyTag()) ==
        CryptoContextImpl<Element>::s_evalMultKeyMap.end()) {
        // the key is not found in the map, so the key has to be generated
        const std::vector<EvalKey<Element>>& evalKeys                  = GetScheme()->EvalMultKeysGen(key);
        CryptoContextImpl<Element>::s_evalMultKeyMap[key->GetKeyTag()] = evalKeys;
    }
}

template <typename Element>
void CryptoContextImpl<Element>::ClearEvalMultKeys() {
    CryptoContextImpl<Element>::s_evalMultKeyMap.clear();
}

template <typename Element>
void CryptoContextImpl<Element>::ClearEvalMultKeys(const std::string& id) {
    auto kd = CryptoContextImpl<Element>::s_evalMultKeyMap.find(id);
    if (kd != CryptoContextImpl<Element>::s_evalMultKeyMap.end())
        CryptoContextImpl<Element>::s_evalMultKeyMap.erase(kd);
}

template <typename Element>
void CryptoContextImpl<Element>::ClearEvalMultKeys(const CryptoContext<Element> cc) {
    for (auto it = CryptoContextImpl<Element>::s_evalMultKeyMap.begin();
         it != CryptoContextImpl<Element>::s_evalMultKeyMap.end();) {
        if (it->second[0]->GetCryptoContext() == cc) {
            it = CryptoContextImpl<Element>::s_evalMultKeyMap.erase(it);
        }
        else {
            ++it;
        }
    }
}

template <typename Element>
void CryptoContextImpl<Element>::InsertEvalMultKey(const std::vector<EvalKey<Element>>& vectorToInsert,
                                                   const std::string& keyTag) {
    const std::string tag = (keyTag.empty()) ? vectorToInsert[0]->GetKeyTag() : keyTag;
    if (CryptoContextImpl<Element>::s_evalMultKeyMap.find(tag) != CryptoContextImpl<Element>::s_evalMultKeyMap.end()) {
        // we do not allow to override the existing key vector if its keyTag is identical to the keyTag of the new keys
        OPENFHE_THROW("Can not save a EvalMultKeys vector as there is a key vector for the given keyTag");
    }

    CryptoContextImpl<Element>::s_evalMultKeyMap[tag] = vectorToInsert;
}

/////////////////////////////////////////
// ADVANCED SHE
/////////////////////////////////////////

template <typename Element>
void CryptoContextImpl<Element>::EvalSumKeyGen(const PrivateKey<Element> privateKey,
                                               const PublicKey<Element> publicKey) {
    ValidateKey(privateKey);
    if (publicKey != nullptr && privateKey->GetKeyTag() != publicKey->GetKeyTag()) {
        OPENFHE_THROW("Public key passed to EvalSumKeyGen does not match private key");
    }

    auto evalKeys = GetScheme()->EvalSumKeyGen(privateKey, publicKey);
    CryptoContextImpl<Element>::InsertEvalAutomorphismKey(evalKeys, privateKey->GetKeyTag());
}

template <typename Element>
std::shared_ptr<std::map<usint, EvalKey<Element>>> CryptoContextImpl<Element>::EvalSumRowsKeyGen(
    const PrivateKey<Element> privateKey, const PublicKey<Element> publicKey, usint rowSize, usint subringDim) {
    ValidateKey(privateKey);
    if (publicKey != nullptr && privateKey->GetKeyTag() != publicKey->GetKeyTag()) {
        OPENFHE_THROW("Public key passed to EvalSumKeyGen does not match private key");
    }

    std::vector<usint> indices;
    auto evalKeys = GetScheme()->EvalSumRowsKeyGen(privateKey, rowSize, subringDim, indices);
    CryptoContextImpl<Element>::InsertEvalAutomorphismKey(evalKeys, privateKey->GetKeyTag());

    return CryptoContextImpl<Element>::GetPartialEvalAutomorphismKeyMapPtr(privateKey->GetKeyTag(), indices);
}

template <typename Element>
std::shared_ptr<std::map<usint, EvalKey<Element>>> CryptoContextImpl<Element>::EvalSumColsKeyGen(
    const PrivateKey<Element> privateKey, const PublicKey<Element> publicKey) {
    ValidateKey(privateKey);
    if (publicKey != nullptr && privateKey->GetKeyTag() != publicKey->GetKeyTag()) {
        OPENFHE_THROW("Public key passed to EvalSumKeyGen does not match private key");
    }

    std::vector<usint> indices;
    auto evalKeys = GetScheme()->EvalSumColsKeyGen(privateKey, indices);
    CryptoContextImpl<Element>::InsertEvalAutomorphismKey(evalKeys, privateKey->GetKeyTag());

    return CryptoContextImpl<Element>::GetPartialEvalAutomorphismKeyMapPtr(privateKey->GetKeyTag(), indices);
}

template <typename Element>
const std::map<usint, EvalKey<Element>>& CryptoContextImpl<Element>::GetEvalSumKeyMap(const std::string& keyID) {
    return CryptoContextImpl<Element>::GetEvalAutomorphismKeyMap(keyID);
}

template <typename Element>
std::map<std::string, std::vector<EvalKey<Element>>>& CryptoContextImpl<Element>::GetAllEvalMultKeys() {
    return CryptoContextImpl<Element>::s_evalMultKeyMap;
}

template <typename Element>
const std::vector<EvalKey<Element>>& CryptoContextImpl<Element>::GetEvalMultKeyVector(const std::string& keyID) {
    auto ekv = CryptoContextImpl<Element>::s_evalMultKeyMap.find(keyID);
    if (ekv == CryptoContextImpl<Element>::s_evalMultKeyMap.end()) {
        std::string errMsg(std::string("Call EvalMultKeyGen() to have EvalMultKey available for ID [") + keyID + "].");
        OPENFHE_THROW(errMsg);
    }
    return ekv->second;
}

template <typename Element>
std::map<std::string, std::shared_ptr<std::map<usint, EvalKey<Element>>>>&
CryptoContextImpl<Element>::GetAllEvalAutomorphismKeys() {
    return CryptoContextImpl<Element>::s_evalAutomorphismKeyMap;
}

template <typename Element>
std::shared_ptr<std::map<usint, EvalKey<Element>>> CryptoContextImpl<Element>::GetEvalAutomorphismKeyMapPtr(
    const std::string& keyID) {
    auto ekv = CryptoContextImpl<Element>::s_evalAutomorphismKeyMap.find(keyID);
    if (ekv == CryptoContextImpl<Element>::s_evalAutomorphismKeyMap.end()) {
        OPENFHE_THROW("EvalAutomorphismKeys are not generated for ID [" + keyID + "].");
    }
    return ekv->second;
}

template <typename Element>
std::shared_ptr<std::map<usint, EvalKey<Element>>> CryptoContextImpl<Element>::GetPartialEvalAutomorphismKeyMapPtr(
    const std::string& keyID, const std::vector<uint32_t>& indexList) {
    if (!indexList.size())
        OPENFHE_THROW("indexList is empty");

    std::shared_ptr<std::map<usint, EvalKey<Element>>> keyMap =
        CryptoContextImpl<Element>::GetEvalAutomorphismKeyMapPtr(keyID);

    // create a return map if specific indices are provided
    std::map<usint, EvalKey<Element>> retMap;
    for (uint32_t indx : indexList) {
        const auto it = keyMap->find(indx);
        if (it == keyMap->end()) {
            OPENFHE_THROW("Key is not generated for index [" + std::to_string(indx) + "] and keyID [" + keyID + "]");
        }
        retMap.emplace(indx, it->second);
    }
    return std::make_shared<std::map<usint, EvalKey<Element>>>(retMap);
}

template <typename Element>
std::map<std::string, std::shared_ptr<std::map<usint, EvalKey<Element>>>>&
CryptoContextImpl<Element>::GetAllEvalSumKeys() {
    return CryptoContextImpl<Element>::GetAllEvalAutomorphismKeys();
}

template <typename Element>
void CryptoContextImpl<Element>::ClearEvalSumKeys() {
    CryptoContextImpl<Element>::ClearEvalAutomorphismKeys();
}

/**
 * ClearEvalMultKeys - flush EvalMultKey cache for a given id
 * @param id
 */
template <typename Element>
void CryptoContextImpl<Element>::ClearEvalSumKeys(const std::string& id) {
    CryptoContextImpl<Element>::ClearEvalAutomorphismKeys(id);
}

/**
 * ClearEvalMultKeys - flush EvalMultKey cache for a given context
 * @param cc
 */
template <typename Element>
void CryptoContextImpl<Element>::ClearEvalSumKeys(const CryptoContext<Element> cc) {
    CryptoContextImpl<Element>::ClearEvalAutomorphismKeys(cc);
}

/////////////////////////////////////////
// SHE AUTOMORPHISM
/////////////////////////////////////////

template <typename Element>
void CryptoContextImpl<Element>::EvalAtIndexKeyGen(const PrivateKey<Element> privateKey,
                                                   const std::vector<int32_t>& indexList,
                                                   const PublicKey<Element> publicKey) {
    ValidateKey(privateKey);
    if (publicKey != nullptr && privateKey->GetKeyTag() != publicKey->GetKeyTag()) {
        OPENFHE_THROW("Public key passed to EvalAtIndexKeyGen does not match private key");
    }

    auto evalKeys = GetScheme()->EvalAtIndexKeyGen(publicKey, privateKey, indexList);
    CryptoContextImpl<Element>::InsertEvalAutomorphismKey(evalKeys, privateKey->GetKeyTag());
}

template <typename Element>
void CryptoContextImpl<Element>::ClearEvalAutomorphismKeys() {
    CryptoContextImpl<Element>::s_evalAutomorphismKeyMap.clear();
}

/**
 * ClearEvalAutomorphismKeys - flush EvalAutomorphismKey cache for a given id
 * @param id
 */
template <typename Element>
void CryptoContextImpl<Element>::ClearEvalAutomorphismKeys(const std::string& id) {
    auto kd = CryptoContextImpl<Element>::s_evalAutomorphismKeyMap.find(id);
    if (kd != CryptoContextImpl<Element>::s_evalAutomorphismKeyMap.end())
        CryptoContextImpl<Element>::s_evalAutomorphismKeyMap.erase(kd);
}

/**
 * ClearEvalAutomorphismKeys - flush EvalAutomorphismKey cache for a given
 * context
 * @param cc
 */
template <typename Element>
void CryptoContextImpl<Element>::ClearEvalAutomorphismKeys(const CryptoContext<Element> cc) {
    for (auto it = CryptoContextImpl<Element>::s_evalAutomorphismKeyMap.begin();
         it != CryptoContextImpl<Element>::s_evalAutomorphismKeyMap.end();) {
        if (it->second->begin()->second->GetCryptoContext() == cc) {
            it = CryptoContextImpl<Element>::s_evalAutomorphismKeyMap.erase(it);
        }
        else {
            ++it;
        }
    }
}

template <typename Element>
std::set<uint32_t> CryptoContextImpl<Element>::GetExistingEvalAutomorphismKeyIndices(const std::string& keyTag) {
    auto keyMapIt = CryptoContextImpl<Element>::s_evalAutomorphismKeyMap.find(keyTag);
    if (keyMapIt == CryptoContextImpl<Element>::s_evalAutomorphismKeyMap.end())
        // there is no keys for the given id, return empty vector
        return std::set<uint32_t>();

    // get all inidices from the existing automorphism key map
    auto& keyMap = *(keyMapIt->second);
    std::set<uint32_t> indices;
    for (const auto& [key, _] : keyMap) {
        indices.insert(key);
    }

    return indices;
}

template <typename Element>
std::set<uint32_t> CryptoContextImpl<Element>::GetUniqueValues(const std::set<uint32_t>& oldValues,
                                                               const std::set<uint32_t>& newValues) {
    std::set<uint32_t> newUniqueValues;
    std::set_difference(newValues.begin(), newValues.end(), oldValues.begin(), oldValues.end(),
                        std::inserter(newUniqueValues, newUniqueValues.begin()));

    return newUniqueValues;
}

template <typename Element>
void CryptoContextImpl<Element>::InsertEvalAutomorphismKey(
    const std::shared_ptr<std::map<uint32_t, EvalKey<Element>>> mapToInsert, const std::string& keyTag) {
    // check if the map is empty
    if (mapToInsert->empty()) {
        return;
    }

    auto mapToInsertIt   = mapToInsert->begin();
    const std::string id = (keyTag.empty()) ? mapToInsertIt->second->GetKeyTag() : keyTag;
    std::set<uint32_t> existingIndices{CryptoContextImpl<Element>::GetExistingEvalAutomorphismKeyIndices(id)};
    if (existingIndices.empty()) {
        // there is no keys for the given id, so we insert full mapToInsert
        CryptoContextImpl<Element>::s_evalAutomorphismKeyMap[id] = mapToInsert;
    }
    else {
        // get all indices from mapToInsert
        std::set<uint32_t> newIndices;
        for (const auto& [key, _] : *mapToInsert) {
            newIndices.insert(key);
        }

        // find all indices in mapToInsert that are not in the exising map and
        // insert those new indices and their corresponding keys to the existing map
        std::set<uint32_t> indicesToInsert{CryptoContextImpl<Element>::GetUniqueValues(existingIndices, newIndices)};
        auto keyMapIt = CryptoContextImpl<Element>::s_evalAutomorphismKeyMap.find(id);
        auto& keyMap  = *(keyMapIt->second);
        for (uint32_t indx : indicesToInsert) {
            keyMap[indx] = (*mapToInsert)[indx];
        }
    }
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalSum(ConstCiphertext<Element> ciphertext, usint batchSize) const {
    ValidateCiphertext(ciphertext);

    auto evalSumKeys = CryptoContextImpl<Element>::GetEvalAutomorphismKeyMap(ciphertext->GetKeyTag());
    return GetScheme()->EvalSum(ciphertext, batchSize, evalSumKeys);
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalSumRows(ConstCiphertext<Element> ciphertext, usint numRows,
                                                            const std::map<usint, EvalKey<Element>>& evalSumKeys,
                                                            usint subringDim) const {
    ValidateCiphertext(ciphertext);

    return GetScheme()->EvalSumRows(ciphertext, numRows, evalSumKeys, subringDim);
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalSumCols(
    ConstCiphertext<Element> ciphertext, usint numCols,
    const std::map<usint, EvalKey<Element>>& evalSumKeysRight) const {
    ValidateCiphertext(ciphertext);

    auto evalSumKeys = CryptoContextImpl<Element>::GetEvalAutomorphismKeyMap(ciphertext->GetKeyTag());
    return GetScheme()->EvalSumCols(ciphertext, numCols, evalSumKeys, evalSumKeysRight);
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalAtIndex(ConstCiphertext<Element> ciphertext, int32_t index) const {
    ValidateCiphertext(ciphertext);

    // If the index is zero, no rotation is needed, copy the ciphertext and return
    // This is done after the keyMap so that it is protected if there's not a valid key.
    if (0 == index) {
        return ciphertext->Clone();
    }

    auto evalAutomorphismKeys = CryptoContextImpl<Element>::GetEvalAutomorphismKeyMap(ciphertext->GetKeyTag());
    return GetScheme()->EvalAtIndex(ciphertext, index, evalAutomorphismKeys);
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalMerge(
    const std::vector<Ciphertext<Element>>& ciphertextVector) const {
    ValidateCiphertext(ciphertextVector[0]);

    auto evalAutomorphismKeys = CryptoContextImpl<Element>::GetEvalAutomorphismKeyMap(ciphertextVector[0]->GetKeyTag());
    return GetScheme()->EvalMerge(ciphertextVector, evalAutomorphismKeys);
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalInnerProduct(ConstCiphertext<Element> ct1,
                                                                 ConstCiphertext<Element> ct2, usint batchSize) const {
    ValidateCiphertext(ct1);
    if (ct2 == nullptr || ct1->GetKeyTag() != ct2->GetKeyTag())
        OPENFHE_THROW("Information was not generated with this crypto context");

    auto evalSumKeys = CryptoContextImpl<Element>::GetEvalAutomorphismKeyMap(ct1->GetKeyTag());
    auto ek          = CryptoContextImpl<Element>::GetEvalMultKeyVector(ct1->GetKeyTag());
    return GetScheme()->EvalInnerProduct(ct1, ct2, batchSize, evalSumKeys, ek[0]);
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalInnerProduct(ConstCiphertext<Element> ct1, ConstPlaintext ct2,
                                                                 usint batchSize) const {
    ValidateCiphertext(ct1);
    if (ct2 == nullptr)
        OPENFHE_THROW("Information was not generated with this crypto context");

    auto evalSumKeys = CryptoContextImpl<Element>::GetEvalAutomorphismKeyMap(ct1->GetKeyTag());
    return GetScheme()->EvalInnerProduct(ct1, ct2, batchSize, evalSumKeys);
}

template <typename Element>
Plaintext CryptoContextImpl<Element>::GetPlaintextForDecrypt(PlaintextEncodings pte, std::shared_ptr<ParmType> evp,
                                                             EncodingParams ep) {
    auto vp = std::make_shared<typename NativePoly::Params>(evp->GetCyclotomicOrder(), ep->GetPlaintextModulus(), 1);

    if (pte == CKKS_PACKED_ENCODING)
        return PlaintextFactory::MakePlaintext(pte, evp, ep);

    return PlaintextFactory::MakePlaintext(pte, vp, ep);
}

template <typename Element>
DecryptResult CryptoContextImpl<Element>::Decrypt(ConstCiphertext<Element> ciphertext,
                                                  const PrivateKey<Element> privateKey, Plaintext* plaintext) {
    if (ciphertext == nullptr)
        OPENFHE_THROW("ciphertext is empty");
    if (plaintext == nullptr)
        OPENFHE_THROW("plaintext is empty");
    ValidateKey(privateKey);

    // determine which type of plaintext that you need to decrypt into
    // Plaintext decrypted =
    // CryptoContextImpl<Element>::GetPlaintextForDecrypt(ciphertext->GetEncodingType(),
    // this->GetElementParams(), this->GetEncodingParams());
    Plaintext decrypted = CryptoContextImpl<Element>::GetPlaintextForDecrypt(
        ciphertext->GetEncodingType(), ciphertext->GetElements()[0].GetParams(), this->GetEncodingParams());

    DecryptResult result;

    if ((ciphertext->GetEncodingType() == CKKS_PACKED_ENCODING) && (typeid(Element) != typeid(NativePoly))) {
        result = GetScheme()->Decrypt(ciphertext, privateKey, &decrypted->GetElement<Poly>());
    }
    else {
        result = GetScheme()->Decrypt(ciphertext, privateKey, &decrypted->GetElement<NativePoly>());
    }

    if (result.isValid == false)  // TODO (dsuponit): why don't we throw an exception here?
        return result;

    decrypted->SetScalingFactorInt(result.scalingFactorInt);

    if (ciphertext->GetEncodingType() == CKKS_PACKED_ENCODING) {
        auto decryptedCKKS = std::dynamic_pointer_cast<CKKSPackedEncoding>(decrypted);
        decryptedCKKS->SetNoiseScaleDeg(ciphertext->GetNoiseScaleDeg());
        decryptedCKKS->SetLevel(ciphertext->GetLevel());
        decryptedCKKS->SetScalingFactor(ciphertext->GetScalingFactor());
        decryptedCKKS->SetSlots(ciphertext->GetSlots());

        const auto cryptoParamsCKKS = std::dynamic_pointer_cast<CryptoParametersRNS>(this->GetCryptoParameters());

        decryptedCKKS->Decode(ciphertext->GetNoiseScaleDeg(), ciphertext->GetScalingFactor(),
                              cryptoParamsCKKS->GetScalingTechnique(), cryptoParamsCKKS->GetExecutionMode(),
                              this->GetFbts());
    }
    else {
        decrypted->Decode();
    }

    *plaintext = std::move(decrypted);
    return result;
}

//------------------------------------------------------------------------------
// Advanced SHE CHEBYSHEV SERIES EXAMPLES
//------------------------------------------------------------------------------

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalChebyshevFunction(std::function<double(double)> func,
                                                                      ConstCiphertext<Element> ciphertext, double a,
                                                                      double b, uint32_t degree) const {
    std::vector<double> coefficients = EvalChebyshevCoefficients(func, a, b, degree);
    return EvalChebyshevSeries(ciphertext, coefficients, a, b);
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalChebyshevFunction(std::function<std::complex<double>(double)> func,
                                                                      ConstCiphertext<Element> ciphertext, double a,
                                                                      double b, uint32_t degree) const {
    std::vector<std::complex<double>> coefficients = EvalChebyshevCoefficients(func, a, b, degree);
    return EvalChebyshevSeries(ciphertext, coefficients, a, b);
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalSin(ConstCiphertext<Element> ciphertext, double a, double b,
                                                        uint32_t degree) const {
    std::function<double(double)> func = [](double x) -> double {
        return std::sin(x);
    };
    return EvalChebyshevFunction(func, ciphertext, a, b, degree);
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalCos(ConstCiphertext<Element> ciphertext, double a, double b,
                                                        uint32_t degree) const {
    std::function<double(double)> func = [](double x) -> double {
        return std::cos(x);
    };
    return EvalChebyshevFunction(func, ciphertext, a, b, degree);
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalLogistic(ConstCiphertext<Element> ciphertext, double a, double b,
                                                             uint32_t degree) const {
    std::function<double(double)> func = [](double x) -> double {
        return std::exp(x*16);
    };
    return EvalChebyshevFunction(func, ciphertext, a, b, degree);
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalDivide(ConstCiphertext<Element> ciphertext, double a, double b,
                                                           uint32_t degree) const {
    std::function<double(double)> func = [](double x) -> double {
        return 1 / x;
    };
    return EvalChebyshevFunction(func, ciphertext, a, b, degree);
}

//------------------------------------------------------------------------------
// Advanced SHE HERMITE
//------------------------------------------------------------------------------

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalHermiteFunction(std::function<double(double)> func,
                                                                    ConstCiphertext<Element> ciphertext, int p,
                                                                    int order) const {
    std::vector<std::complex<double>> coefficients;

    switch (order) {
        case 1:
            coefficients = hermiteInterpOrder1Coeff(func, p);
            break;
        case 2:
            coefficients = hermiteInterpOrder2Coeff(func, p);
            break;
        case 3:
            coefficients = hermiteInterpOrder3Coeff(func, p);
            break;
        default:
            OPENFHE_THROW("Order of Hermite interpolation must be between 1 and 3");
            break;
    }

    return EvalPoly(ciphertext, coefficients);
}

// template <typename Element>
// Ciphertext<Element> CryptoContextImpl<Element>::EvalHermiteFunction_kx(std::function<double(double)> func,
//                                                                        ConstCiphertext<Element> ciphertext, int p,
//                                                                        int order) const {
//     std::vector<std::complex<double>> coefficients;
//     std::chrono::system_clock::time_point start, end;
//     auto cc = ciphertext->GetCryptoContext();
//     start   = std::chrono::system_clock::now();
//     switch (order) {
//         case 1:
//             coefficients = hermiteInterpOrder1Coeff(func, p);
//             break;
//         case 2:
//             coefficients = hermiteInterpOrder2Coeff(func, p);
//             break;
//         case 3:
//             coefficients = hermiteInterpOrder3Coeff(func, p);
//             break;
//         default:
//             OPENFHE_THROW("Order of Hermite interpolation must be between 1 and 3");
//             break;
//     }
//     end             = std::chrono::system_clock::now();
//     auto time_total = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

//     std::cout << "hermit coe time: " << time_total << " ms" << std::endl;
//     start = std::chrono::system_clock::now();
//     std::cout << "hermite deg: " << coefficients.size() << std::endl;
//     // for (size_t i = 0; i < (coefficients.size()); i++) {
//     //     std::cout << coefficients[i] << ", ";
//     // }
//     std::cout << std::endl;
//     Plaintext plain_coe = cc->MakeCKKSPackedPlaintext(coefficients);
//     auto res            = cc->EvalMult(ciphertext, plain_coe);
//     auto poly_degree    = coefficients.size();
//     for (int i = 0; i < std::log2(poly_degree); i++) {
//         cc->EvalAddInPlace(res, cc->EvalRotate(res, pow(2, i)));
//     }
//     end        = std::chrono::system_clock::now();
//     time_total = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

//     std::cout << "hermite_kx time: " << time_total << " ms" << std::endl;
//     return res;
// }

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalHermiteFunction_kx(
    ConstCiphertext<Element> ciphertext, std::vector<std::complex<double>> hermite_coefficients) const {
    auto cc = ciphertext->GetCryptoContext();

    auto start = std::chrono::system_clock::now();
    std::cout << "hermite deg: " << hermite_coefficients.size() << std::endl;
    std::cout << std::endl;

    Plaintext plain_coe = cc->MakeCKKSPackedPlaintext(hermite_coefficients);
    auto res            = cc->EvalMult(ciphertext, plain_coe);
    cc->ModReduceInPlace(res);
    auto poly_degree = hermite_coefficients.size();
    for (int i = 0; i < std::log2(poly_degree); i++) {
        cc->EvalAddInPlace(res, cc->EvalRotate(res, pow(2, i)));
    }

    auto end        = std::chrono::system_clock::now();
    auto time_total = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

    std::cout << "hermite_kx time: " << time_total << " ms" << std::endl;
    return res;
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalHermiteFunction_kx_parallel(
    std::vector<Ciphertext<Element>> smallPowers, std::vector<Ciphertext<Element>> bigPowers,
    std::vector<std::complex<double>> hermite_coefficients) const {
    auto cc = smallPowers[0]->GetCryptoContext();

    const size_t m = smallPowers.size();
    const size_t K = hermite_coefficients.size();
    const size_t Q = bigPowers.size();

    Ciphertext<Element> result = cc->EvalMult(smallPowers[0], 0.0);

    std::vector<Ciphertext<Element>> partial(Q);

#pragma omp parallel for schedule(dynamic)
    for (long long q = 0; q < (long long)Q; q++) {
        // inner = 0
        Ciphertext<Element> inner = cc->EvalMult(smallPowers[0], 0.0);

        size_t base = (size_t)q * m;
        for (size_t r = 0; r < m; r++) {
            size_t k = base + r;
            if (k >= K)
                break;

            const auto& coeff = hermite_coefficients[k];
            if (coeff == std::complex<double>(0.0, 0.0))
                continue;

            Ciphertext<Element> term = cc->EvalMult(smallPowers[r], coeff);
            inner                    = cc->EvalAdd(inner, term);
        }

        Ciphertext<Element> outer = cc->EvalMult(bigPowers[q], inner);
        cc->ModReduceInPlace(outer);

        partial[q] = outer;
    }

    for (size_t q = 0; q < Q; q++) {
        result = cc->EvalAdd(result, partial[q]);
    }

    return result;
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalHermiteFunction_kx_parallel_opt(
    std::vector<std::vector<Ciphertext<Element>>> bases, std::vector<std::complex<double>> hermite_coefficients) const {
    const size_t m_val = bases[0].size();
    const size_t depth = bases.size() - 1;
    auto cc            = bases[0][0]->GetCryptoContext();

    std::vector<Ciphertext<Element>> current_ciphers;
    std::vector<char> is_block_zero;
    size_t num_blocks_0 = (hermite_coefficients.size() + m_val - 1) / m_val;
    is_block_zero.resize(num_blocks_0, 0);
    current_ciphers.resize(num_blocks_0);
    int zero = 0;
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    double time_total;
#pragma omp parallel for schedule(dynamic)
    for (size_t b = 0; b < num_blocks_0; ++b) {
        std::vector<Ciphertext<Element>> terms;
        terms.reserve(m_val);

        size_t start_idx                = b * m_val;
        Ciphertext<Element> zero_cipher = cc->EvalMult(bases[0][0], 0.0);
        bool block_is_zero              = true;
        for (size_t i = 0; i < m_val; ++i) {
            size_t global_idx = start_idx + i;
            if (global_idx >= hermite_coefficients.size())
                break;

            
            if (std::abs(hermite_coefficients[global_idx].real()) < 1e-6 &&
                std::abs(hermite_coefficients[global_idx].imag()) < 1e-6) {
                zero++;
                continue;
            }
            block_is_zero = false;

            Ciphertext<Element> term;
            if (i == 0) {
                term = zero_cipher->Clone();
                cc->EvalAddInPlace(term, hermite_coefficients[global_idx]);
            }
            else {
                term = cc->EvalMult(bases[0][i], hermite_coefficients[global_idx]);
            }
            cc->ModReduceInPlace(term);
            terms.push_back(term);
        }

        if (block_is_zero || terms.empty()) {
            current_ciphers[b] = zero_cipher;  
            is_block_zero[b]   = 1;            
        }
        else {
            current_ciphers[b] = cc->EvalAddMany(terms);
            is_block_zero[b]   = 0;
        }
    }
    end        = std::chrono::system_clock::now();
    time_total = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "Phase 1 time: " << time_total << " ms" << std::endl;
    std::cout << "zero coe: " << zero << std::endl;
    // ====================================================
    
    // ====================================================

    start = std::chrono::system_clock::now();
    for (size_t level = 1; level <= depth; ++level) {
        size_t prev_count = current_ciphers.size();
        size_t next_count = (prev_count + m_val - 1) / m_val;
        std::vector<Ciphertext<Element>> next_ciphers(next_count);
        std::vector<char> next_is_zero(next_count, 0);
#pragma omp parallel for schedule(dynamic)
        for (size_t b = 0; b < next_count; ++b) {
            std::vector<Ciphertext<Element>> terms;
            terms.reserve(m_val);
            bool block_all_zeros = true;
            size_t start_idx     = b * m_val;
            Ciphertext<Element> term;
            for (size_t i = 0; i < m_val; ++i) {
                size_t input_idx = start_idx + i;
                if (input_idx >= prev_count)
                    break;
                if (is_block_zero[input_idx])
                    continue;
                block_all_zeros = false;
                if (i == 0)
                    term = cc->EvalMult(current_ciphers[input_idx], 1);
                else
                    term = cc->EvalMult(current_ciphers[input_idx], bases[level][i]);
                terms.push_back(term);
            }

            Ciphertext<Element> block_sum;
            if (block_all_zeros || terms.empty()) {
                block_sum       = cc->EvalMult(bases[level][1], 0.0);
                next_is_zero[b] = 1;
            }
            else {
                block_sum       = cc->EvalAddMany(terms);
                next_is_zero[b] = 0;
            }
            cc->ModReduceInPlace(block_sum);
            next_ciphers[b] = block_sum;
        }
        current_ciphers = std::move(next_ciphers);
        is_block_zero   = std::move(next_is_zero);
    }
    end        = std::chrono::system_clock::now();
    time_total = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "Phase 2 time: " << time_total << " ms" << std::endl;
    // std::cout << "zeros: " << zero << std::endl;
    return current_ciphers[0];
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalHermiteFunction_kx_sparse(
    Ciphertext<Element> basis1, const std::vector<std::complex<double>>& hermite_coefficients, size_t step) const {
    std::complex<double> c0 = (hermite_coefficients.size() > 0) ? hermite_coefficients[0] : 0.0;

    size_t max_k = 0;
    std::vector<size_t> active_k;  
    const double threshold = 1e-9;

    
    for (size_t i = 1; i < hermite_coefficients.size(); ++i) {
        if (std::abs(hermite_coefficients[i].real()) > threshold ||
            std::abs(hermite_coefficients[i].imag()) > threshold) {
            
            if (i % step != 0)
                continue;

            size_t k = i / step;
            active_k.push_back(k);
            if (k > max_k)
                max_k = k;
        }
    }

    // ----------------------------------------------------------------------
    
    // ----------------------------------------------------------------------
    
    
    std::vector<Ciphertext<Element>> powers(max_k + 1);

    
    if (max_k >= 1)
        powers[1] = basis1;

    for (size_t k = 2; k <= max_k; ++k) {
        if ((k & (k - 1)) == 0) {
            
            powers[k] = EvalMult(powers[k / 2], powers[k / 2]);
            RescaleInPlace(powers[k]);
        }
        else {
            
            size_t high_bit = 1;
            while ((high_bit << 1) <= k) {
                high_bit <<= 1;
            }
            size_t remainder = k - high_bit;
            powers[k]        = EvalMult(powers[high_bit], powers[remainder]);
            RescaleInPlace(powers[k]);
        }
    }

    // ----------------------------------------------------------------------
    
    // ----------------------------------------------------------------------
    std::vector<Ciphertext<Element>> terms;
    terms.reserve(active_k.size());

    size_t deepest_level = 0;

    for (size_t k : active_k) {
        size_t coeff_idx       = k * step;
        std::complex<double> c = hermite_coefficients[coeff_idx];

        
        Ciphertext<Element> term = EvalMult(powers[k], c);

        
        
        RescaleInPlace(term);

        if (term->GetLevel() > deepest_level) {
            deepest_level = term->GetLevel();
        }

        terms.push_back(term);
    }

    // ----------------------------------------------------------------------
    
    // ----------------------------------------------------------------------
    for (auto& term : terms) {
        if (term->GetLevel() < deepest_level) {
            
            size_t levels_to_drop = deepest_level - term->GetLevel();
            LevelReduceInPlace(term, nullptr, levels_to_drop);
        }
    }

    // ----------------------------------------------------------------------
    
    // ----------------------------------------------------------------------

    
    Ciphertext<Element> result = EvalAddMany(terms);

    if (std::abs(c0.real()) > 1e-9 || std::abs(c0.imag()) > 1e-9) {
        EvalAddInPlace(result, c0);
    }

    return result;
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalHermiteFunction_kx_sparse_PS(
    Ciphertext<Element> basis1, const std::vector<std::complex<double>>& hermite_coefficients, size_t step) const {
    // ------------------------------------------------------------------
    
    
    
    // ------------------------------------------------------------------
    std::vector<std::complex<double>> dense_coeffs;
    size_t sparse_size = hermite_coefficients.size();

    
    dense_coeffs.reserve(sparse_size / step + 1);

    for (size_t i = 0; i < sparse_size; i += step) {
        dense_coeffs.push_back(hermite_coefficients[i]);
    }

    while (dense_coeffs.size() > 1 && std::abs(dense_coeffs.back().real()) < 1e-9 &&
           std::abs(dense_coeffs.back().imag()) < 1e-9) {
        dense_coeffs.pop_back();
    }
    if (dense_coeffs.size() < 5)
        return EvalPolyLinear(basis1, dense_coeffs);
    else
        return EvalPolyPS(basis1, dense_coeffs);
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalHermitePrecomLinearFunction(std::function<double(double)> func,
                                                                                std::vector<Ciphertext<Element>> powers,
                                                                                int p) const {
    auto coefficients = hermiteInterpOrder1Coeff(func, p);

    return EvalPolyLinear(powers, coefficients);
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::EvalHermitePrecomPSFunction(std::function<double(double)> func,
                                                                            std::vector<Ciphertext<Element>> powers,
                                                                            std::vector<Ciphertext<Element>> powers2,
                                                                            Ciphertext<Element> power2km1,
                                                                            int p) const {
    auto coefficients = hermiteInterpOrder1Coeff(func, p);

    return EvalPolyPS(powers, powers2, power2km1, coefficients);
}

}  // namespace lbcrypto

// the code below is from cryptocontext-impl.cpp
namespace lbcrypto {

template <>
Plaintext CryptoContextImpl<DCRTPoly>::GetPlaintextForDecrypt(PlaintextEncodings pte, std::shared_ptr<ParmType> evp,
                                                              EncodingParams ep) {
    if ((pte == CKKS_PACKED_ENCODING) && (evp->GetParams().size() > 1)) {
        auto vp = std::make_shared<typename Poly::Params>(evp->GetCyclotomicOrder(), ep->GetPlaintextModulus(), 1);
        return PlaintextFactory::MakePlaintext(pte, vp, ep);
    }
    else {
        auto vp =
            std::make_shared<typename NativePoly::Params>(evp->GetCyclotomicOrder(), ep->GetPlaintextModulus(), 1);
        return PlaintextFactory::MakePlaintext(pte, vp, ep);
    }
}

template <>
DecryptResult CryptoContextImpl<DCRTPoly>::Decrypt(ConstCiphertext<DCRTPoly> ciphertext,
                                                   const PrivateKey<DCRTPoly> privateKey, Plaintext* plaintext) {
    if (ciphertext == nullptr)
        OPENFHE_THROW("ciphertext is empty");
    if (plaintext == nullptr)
        OPENFHE_THROW("plaintext is empty");
    if (privateKey == nullptr || Mismatched(privateKey->GetCryptoContext()))
        OPENFHE_THROW("Information was not generated with this crypto context");

    // determine which type of plaintext that you need to decrypt into
    // Plaintext decrypted =
    // CryptoContextImpl<Element>::GetPlaintextForDecrypt(ciphertext->GetEncodingType(),
    // this->GetElementParams(), this->GetEncodingParams());
    Plaintext decrypted = CryptoContextImpl<DCRTPoly>::GetPlaintextForDecrypt(
        ciphertext->GetEncodingType(), ciphertext->GetElements()[0].GetParams(), this->GetEncodingParams());

    DecryptResult result;

    if ((ciphertext->GetEncodingType() == CKKS_PACKED_ENCODING) &&
        (ciphertext->GetElements()[0].GetParams()->GetParams().size() > 1))  // more than one tower in DCRTPoly
        result = GetScheme()->Decrypt(ciphertext, privateKey, &decrypted->GetElement<Poly>());
    else
        result = GetScheme()->Decrypt(ciphertext, privateKey, &decrypted->GetElement<NativePoly>());

    if (result.isValid == false)
        return result;

    decrypted->SetScalingFactorInt(result.scalingFactorInt);

    if (ciphertext->GetEncodingType() == CKKS_PACKED_ENCODING) {
        auto decryptedCKKS = std::dynamic_pointer_cast<CKKSPackedEncoding>(decrypted);
        decryptedCKKS->SetNoiseScaleDeg(ciphertext->GetNoiseScaleDeg());
        decryptedCKKS->SetLevel(ciphertext->GetLevel());
        decryptedCKKS->SetScalingFactor(ciphertext->GetScalingFactor());
        decryptedCKKS->SetSlots(ciphertext->GetSlots());

        const auto cryptoParamsCKKS = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(this->GetCryptoParameters());

        decryptedCKKS->Decode(ciphertext->GetNoiseScaleDeg(), ciphertext->GetScalingFactor(),
                              cryptoParamsCKKS->GetScalingTechnique(), cryptoParamsCKKS->GetExecutionMode(),
                              this->GetFbts());
    }
    else {
        decrypted->Decode();
    }

    *plaintext = std::move(decrypted);
    return result;
}

template <>
DecryptResult CryptoContextImpl<DCRTPoly>::MultipartyDecryptFusion(
    const std::vector<Ciphertext<DCRTPoly>>& partialCiphertextVec, Plaintext* plaintext) const {
    DecryptResult result;

    // Make sure we're processing ciphertexts.
    size_t last_ciphertext = partialCiphertextVec.size();
    if (last_ciphertext < 1)
        return result;

    for (size_t i = 0; i < last_ciphertext; i++) {
        ValidateCiphertext(partialCiphertextVec[i]);
        if (partialCiphertextVec[i]->GetEncodingType() != partialCiphertextVec[0]->GetEncodingType())
            OPENFHE_THROW("Ciphertexts have mismatched encoding types");
    }

    // determine which type of plaintext that you need to decrypt into
    Plaintext decrypted = CryptoContextImpl<DCRTPoly>::GetPlaintextForDecrypt(
        partialCiphertextVec[0]->GetEncodingType(), partialCiphertextVec[0]->GetElements()[0].GetParams(),
        this->GetEncodingParams());

    if ((partialCiphertextVec[0]->GetEncodingType() == CKKS_PACKED_ENCODING) &&
        (partialCiphertextVec[0]->GetElements()[0].GetParams()->GetParams().size() > 1))
        result = GetScheme()->MultipartyDecryptFusion(partialCiphertextVec, &decrypted->GetElement<Poly>());
    else
        result = GetScheme()->MultipartyDecryptFusion(partialCiphertextVec, &decrypted->GetElement<NativePoly>());

    if (result.isValid == false)
        return result;

    decrypted->SetScalingFactorInt(result.scalingFactorInt);

    if (partialCiphertextVec[0]->GetEncodingType() == CKKS_PACKED_ENCODING) {
        auto decryptedCKKS = std::dynamic_pointer_cast<CKKSPackedEncoding>(decrypted);
        decryptedCKKS->SetSlots(partialCiphertextVec[0]->GetSlots());
        const auto cryptoParamsCKKS = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(this->GetCryptoParameters());
        decryptedCKKS->Decode(partialCiphertextVec[0]->GetNoiseScaleDeg(), partialCiphertextVec[0]->GetScalingFactor(),
                              cryptoParamsCKKS->GetScalingTechnique(), cryptoParamsCKKS->GetExecutionMode(),
                              this->GetFbts());
    }
    else {
        decrypted->Decode();
    }

    *plaintext = std::move(decrypted);

    return result;
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::IntMPBootAdjustScale(ConstCiphertext<Element> ciphertext) const {
    return GetScheme()->IntMPBootAdjustScale(ciphertext);
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::IntMPBootRandomElementGen(const PublicKey<Element> publicKey) const {
    const auto cryptoParamsCKKS = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(this->GetCryptoParameters());
    return GetScheme()->IntMPBootRandomElementGen(cryptoParamsCKKS, publicKey);
}

template <typename Element>
std::vector<Ciphertext<Element>> CryptoContextImpl<Element>::IntMPBootDecrypt(const PrivateKey<Element> privateKey,
                                                                              ConstCiphertext<Element> ciphertext,
                                                                              ConstCiphertext<Element> a) const {
    return GetScheme()->IntMPBootDecrypt(privateKey, ciphertext, a);
}

template <typename Element>
std::vector<Ciphertext<Element>> CryptoContextImpl<Element>::IntMPBootAdd(
    std::vector<std::vector<Ciphertext<Element>>>& sharesPairVec) const {
    return GetScheme()->IntMPBootAdd(sharesPairVec);
}

template <typename Element>
Ciphertext<Element> CryptoContextImpl<Element>::IntMPBootEncrypt(const PublicKey<Element> publicKey,
                                                                 const std::vector<Ciphertext<Element>>& sharesPair,
                                                                 ConstCiphertext<Element> a,
                                                                 ConstCiphertext<Element> ciphertext) const {
    return GetScheme()->IntMPBootEncrypt(publicKey, sharesPair, a, ciphertext);
}

// Function for sharing and recovery of secret for Threshold FHE with aborts
template <>
std::unordered_map<uint32_t, DCRTPoly> CryptoContextImpl<DCRTPoly>::ShareKeys(const PrivateKey<DCRTPoly>& sk, usint N,
                                                                              usint threshold, usint index,
                                                                              const std::string& shareType) const {
    // conditions on N and threshold for security with aborts
    if (N < 2)
        OPENFHE_THROW("Number of parties needs to be at least 3 for aborts");

    if (threshold <= N / 2)
        OPENFHE_THROW("Threshold required to be majority (more than N/2)");

    const auto cryptoParams = sk->GetCryptoContext()->GetCryptoParameters();
    auto elementParams      = cryptoParams->GetElementParams();
    auto vecSize            = elementParams->GetParams().size();
    auto ring_dimension     = elementParams->GetRingDimension();

    // condition for inverse in lagrange coeff to exist.
    for (size_t i = 0; i < vecSize; ++i) {
        auto modq_k = elementParams->GetParams()[i]->GetModulus();
        if (N >= modq_k)
            OPENFHE_THROW("Number of parties N needs to be less than DCRTPoly moduli");
    }

    // secret sharing
    std::unordered_map<uint32_t, DCRTPoly> SecretShares;

    if (shareType == "additive") {
        // generate a random share of N-2 elements and create the last share as sk - (sk_1 + ... + sk_N-2)
        typename DCRTPoly::DugType dug;
        DCRTPoly rsum(dug, elementParams, Format::EVALUATION);

        const uint32_t num_of_shares = N - 1;
        std::vector<DCRTPoly> SecretSharesVec;
        SecretSharesVec.reserve(num_of_shares);
        SecretSharesVec.push_back(rsum);
        for (size_t i = 1; i < num_of_shares - 1; ++i) {
            DCRTPoly r(dug, elementParams, Format::EVALUATION);  // should re-generate uniform r for each share
            rsum += r;
            SecretSharesVec.push_back(std::move(r));
        }
        SecretSharesVec.push_back(sk->GetPrivateElement() - rsum);

        for (size_t i = 1, ctr = 0; i <= N; ++i) {
            if (i != index) {
                SecretShares[i] = SecretSharesVec[ctr++];
            }
        }
    }
    else if (shareType == "shamir") {
        // vector to store columnwise randomly generated coefficients for polynomial f from Z_q for every secret key entry
        // set constant term of polynomial f_i to s_i
        std::vector<DCRTPoly> fs{sk->GetPrivateElement()};
        fs.back().SetFormat(Format::COEFFICIENT);

        // generate random coefficients
        fs.reserve(threshold);
        typename DCRTPoly::DugType dug;
        for (size_t i = 1; i < threshold; ++i) {
            fs.emplace_back(dug, elementParams, Format::COEFFICIENT);
        }

        // evaluate the polynomial at the index of the parties 1 to N
        for (size_t i = 1; i <= N; ++i) {
            if (i != index) {
                DCRTPoly feval(elementParams, Format::COEFFICIENT, true);
                for (size_t k = 0; k < vecSize; k++) {
                    auto modq_k = elementParams->GetParams()[k]->GetModulus();

                    NativePoly powtemppoly(elementParams->GetParams()[k], Format::COEFFICIENT);
                    NativePoly fevalpoly(elementParams->GetParams()[k], Format::COEFFICIENT, true);

                    NativeInteger powtemp(1);
                    for (size_t t = 1; t < threshold; t++) {
                        NativeVector powtempvec(ring_dimension, modq_k, (powtemp = powtemp.ModMul(i, modq_k)));

                        powtemppoly.SetValues(std::move(powtempvec), Format::COEFFICIENT);

                        auto& fst = fs[t].GetElementAtIndex(k);

                        for (size_t i = 0; i < ring_dimension; ++i) {
                            fevalpoly[i] += powtemppoly[i].ModMul(fst[i], modq_k);
                        }
                    }
                    fevalpoly += fs[0].GetElementAtIndex(k);

                    fevalpoly.SetFormat(Format::COEFFICIENT);
                    feval.SetElementAtIndex(k, std::move(fevalpoly));
                }
                // assign fi
                SecretShares.emplace(i, std::move(feval));
            }
        }
    }
    return SecretShares;
}

template <>
void CryptoContextImpl<DCRTPoly>::RecoverSharedKey(PrivateKey<DCRTPoly>& sk,
                                                   std::unordered_map<uint32_t, DCRTPoly>& sk_shares, usint N,
                                                   usint threshold, const std::string& shareType) const {
    if (sk_shares.size() < threshold)
        OPENFHE_THROW("Number of shares available less than threshold of the sharing scheme");

    // conditions on N and threshold for security with aborts
    if (N < 2)
        OPENFHE_THROW("Number of parties needs to be at least 3 for aborts");

    if (threshold <= N / 2)
        OPENFHE_THROW("Threshold required to be majority (more than N/2)");

    const auto& cryptoParams  = sk->GetCryptoContext()->GetCryptoParameters();
    const auto& elementParams = cryptoParams->GetElementParams();
    size_t ring_dimension     = elementParams->GetRingDimension();
    size_t vecSize            = elementParams->GetParams().size();

    // condition for inverse in lagrange coeff to exist.
    for (size_t k = 0; k < vecSize; k++) {
        auto modq_k = elementParams->GetParams()[k]->GetModulus();
        if (N >= modq_k)
            OPENFHE_THROW("Number of parties N needs to be less than DCRTPoly moduli");
    }

    // vector of indexes of the clients
    std::vector<uint32_t> client_indexes;
    client_indexes.reserve(N);
    for (uint32_t i = 1; i <= N; ++i) {
        if (sk_shares.find(i) != sk_shares.end())
            client_indexes.push_back(i);
    }
    const uint32_t client_indexes_size = client_indexes.size();

    if (client_indexes_size < threshold)
        OPENFHE_THROW("Not enough shares to recover the secret");

    if (shareType == "additive") {
        DCRTPoly sum_of_elems(elementParams, Format::EVALUATION, true);
        for (uint32_t i = 0; i < threshold; ++i) {
            sum_of_elems += sk_shares[client_indexes[i]];
        }
        sk->SetPrivateElement(std::move(sum_of_elems));
    }
    else if (shareType == "shamir") {
        // use lagrange interpolation to recover the secret
        // vector of lagrange coefficients L_j = Pdt_i ne j (i (i-j)^-1)
        std::vector<DCRTPoly> Lagrange_coeffs(client_indexes_size, DCRTPoly(elementParams, Format::EVALUATION));

        // recovery of the secret with lagrange coefficients and the secret shares
        for (uint32_t j = 0; j < client_indexes_size; j++) {
            auto cj = client_indexes[j];
            for (size_t k = 0; k < vecSize; k++) {
                auto modq_k = elementParams->GetParams()[k]->GetModulus();
                NativePoly multpoly(elementParams->GetParams()[k], Format::COEFFICIENT, true);
                multpoly.AddILElementOne();
                for (uint32_t i = 0; i < client_indexes_size; i++) {
                    auto ci = client_indexes[i];
                    if (ci != cj) {
                        auto&& denominator = (cj < ci) ? NativeInteger(ci - cj) : modq_k - NativeInteger(cj - ci);
                        auto denom_inv{denominator.ModInverse(modq_k)};
                        for (size_t d = 0; d < ring_dimension; ++d)
                            multpoly[d].ModMulFastEq(NativeInteger(ci).ModMul(denom_inv, modq_k), modq_k);
                    }
                }
                multpoly.SetFormat(Format::EVALUATION);
                Lagrange_coeffs[j].SetElementAtIndex(k, std::move(multpoly));
            }
            Lagrange_coeffs[j].SetFormat(Format::COEFFICIENT);
        }

        DCRTPoly lagrange_sum_of_elems(elementParams, Format::COEFFICIENT, true);
        for (size_t k = 0; k < vecSize; ++k) {
            NativePoly lagrange_sum_of_elems_poly(elementParams->GetParams()[k], Format::COEFFICIENT, true);
            for (uint32_t i = 0; i < client_indexes_size; ++i) {
                const auto& coeff = Lagrange_coeffs[i].GetAllElements()[k];
                const auto& share = sk_shares[client_indexes[i]].GetAllElements()[k];
                lagrange_sum_of_elems_poly += coeff.TimesNoCheck(share);
            }
            lagrange_sum_of_elems.SetElementAtIndex(k, std::move(lagrange_sum_of_elems_poly));
        }
        lagrange_sum_of_elems.SetFormat(Format::EVALUATION);
        sk->SetPrivateElement(std::move(lagrange_sum_of_elems));
    }
}

template class CryptoContextImpl<DCRTPoly>;

}  // namespace lbcrypto
