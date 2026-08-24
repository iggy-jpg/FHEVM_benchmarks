#include "HomOP.h"
#include "HomUtil.h"
#include "HomPKE.h"
#include "HomDecompress.h"
#include <mutex>
#include <vector>            

namespace TFHEpp {

namespace {
constexpr TLWE<lvl2param> trivial_constant(bool value)
{
    TLWE<lvl2param> result{};
    result[lvl2param::n] = value ? lvl2param::μ : -lvl2param::μ;
    return result;
}

std::mutex constants_mutex;
}

TLWE<lvl2param> CT_ZERO = trivial_constant(false);
TLWE<lvl2param> CT_ONE = trivial_constant(true);

void init_constants(const SecretKey &sk)
{
    
    // CT_ZERO = bootsSymEncrypt<lvl2param>(0, sk);
    // CT_ONE  = bootsSymEncrypt<lvl2param>(1, sk);

    
    std::vector<uint8_t> v0{0};                                  // [CHANGED]
    std::vector<uint8_t> v1{1};                                  // [CHANGED]
    auto zero = bootsSymEncrypt<lvl2param>(v0, sk).at(0);
    auto one = bootsSymEncrypt<lvl2param>(v1, sk).at(0);
    std::lock_guard<std::mutex> lock(constants_mutex);
    CT_ZERO = std::move(zero);
    CT_ONE = std::move(one);
}

void init_constants(const PKE::SeededRlwePkPair<lvl2param>&) {
    auto zero = trivial_constant(false);
    auto one = trivial_constant(true);
    std::lock_guard<std::mutex> lock(constants_mutex);
    CT_ZERO = std::move(zero);
    CT_ONE = std::move(one);
}
} // namespace TFHEpp
