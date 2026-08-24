#pragma once

#include "ARCEDB/utils/types.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <cstdint>
#include <vector>

namespace arcedb_integer {

using Bit = arcedb::TLWELvl0;
using Word = std::vector<Bit>; // little-endian
using boost::multiprecision::cpp_int;

inline Bit Constant(bool value)
{
    Bit out;
    if (value) TFHEpp::HomCONSTANTONE(out);
    else TFHEpp::HomCONSTANTZERO(out);
    return out;
}

inline Bit And(const Bit& a, const Bit& b, const arcedb::TFHEEvalKey& ek)
{
    Bit out; TFHEpp::HomAND(out, a, b, ek); return out;
}

inline Bit Or(const Bit& a, const Bit& b, const arcedb::TFHEEvalKey& ek)
{
    Bit out; TFHEpp::HomOR(out, a, b, ek); return out;
}

inline Bit Xor(const Bit& a, const Bit& b, const arcedb::TFHEEvalKey& ek)
{
    Bit out; TFHEpp::HomXOR(out, a, b, ek); return out;
}

inline Bit Not(const Bit& a)
{
    Bit out; TFHEpp::HomNOT(out, a); return out;
}

inline Bit Select(const Bit& condition, const Bit& ifTrue, const Bit& ifFalse,
                  const arcedb::TFHEEvalKey& ek)
{
    Bit out; TFHEpp::HomMUX(out, condition, ifTrue, ifFalse, ek); return out;
}

inline Word Encrypt(cpp_int value, int width, const arcedb::TFHESecretKey& sk)
{
    const cpp_int mask = (cpp_int(1) << width) - 1;
    value &= mask;
    Word out(width);
    for (int i = 0; i < width; ++i) {
        const std::uint8_t bit = static_cast<std::uint8_t>(((value >> i) & 1) != 0);
        out[i] = TFHEpp::bootsSymEncrypt<arcedb::Lvl0>({bit}, sk)[0];
    }
    return out;
}

inline cpp_int Decrypt(const Word& value, const arcedb::TFHESecretKey& sk)
{
    cpp_int out = 0;
    for (std::size_t i = 0; i < value.size(); ++i)
        if (TFHEpp::bootsSymDecrypt<arcedb::Lvl0>({value[i]}, sk)[0] & 1)
            out |= cpp_int(1) << i;
    return out;
}

inline bool Decrypt(const Bit& value, const arcedb::TFHESecretKey& sk)
{
    return (TFHEpp::bootsSymDecrypt<arcedb::Lvl0>({value}, sk)[0] & 1) != 0;
}

inline Word Add(const Word& a, const Word& b, const arcedb::TFHEEvalKey& ek)
{
    const int width = static_cast<int>(a.size());
    Word out(width);
    Bit carry = Constant(false);
    for (int i = 0; i < width; ++i) {
        const Bit axb = Xor(a[i], b[i], ek);
        out[i] = Xor(axb, carry, ek);
        carry = Or(And(a[i], b[i], ek), And(carry, axb, ek), ek);
    }
    return out;
}

inline Word Sub(const Word& a, const Word& b, const arcedb::TFHEEvalKey& ek)
{
    const int width = static_cast<int>(a.size());
    Word out(width);
    Bit borrow = Constant(false);
    for (int i = 0; i < width; ++i) {
        const Bit axb = Xor(a[i], b[i], ek);
        out[i] = Xor(axb, borrow, ek);
        const Bit generate = And(Not(a[i]), b[i], ek);
        const Bit propagate = And(Not(axb), borrow, ek);
        borrow = Or(generate, propagate, ek);
    }
    return out;
}

inline Word Multiply(const Word& a, const Word& b, const arcedb::TFHEEvalKey& ek)
{
    const int width = static_cast<int>(a.size());
    Word result(width, Constant(false));
    for (int row = 0; row < width; ++row) {
        Word partial(width, Constant(false));
        for (int bit = row; bit < width; ++bit)
            partial[bit] = And(a[bit - row], b[row], ek);
        result = Add(result, partial, ek);
    }
    return result;
}

struct Comparison { Bit gt, lt, eq, ge, le; };

inline Comparison Compare(const Word& a, const Word& b,
                          const arcedb::TFHEEvalKey& ek)
{
    Bit gt = Constant(false), lt = Constant(false), equal = Constant(true);
    for (int i = static_cast<int>(a.size()) - 1; i >= 0; --i) {
        const Bit aGtB = And(a[i], Not(b[i]), ek);
        const Bit aLtB = And(Not(a[i]), b[i], ek);
        gt = Or(gt, And(equal, aGtB, ek), ek);
        lt = Or(lt, And(equal, aLtB, ek), ek);
        equal = And(equal, Not(Xor(a[i], b[i], ek)), ek);
    }
    return {gt, lt, equal, Or(gt, equal, ek), Or(lt, equal, ek)};
}

inline Word Select(const Bit& condition, const Word& ifTrue, const Word& ifFalse,
                   const arcedb::TFHEEvalKey& ek)
{
    Word out(ifTrue.size());
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = Select(condition, ifTrue[i], ifFalse[i], ek);
    return out;
}

inline Word ZeroExtend(const Word& value, int width)
{
    Word out(width, Constant(false));
    for (std::size_t i = 0; i < value.size() && i < out.size(); ++i) out[i] = value[i];
    return out;
}

} // namespace arcedb_integer
