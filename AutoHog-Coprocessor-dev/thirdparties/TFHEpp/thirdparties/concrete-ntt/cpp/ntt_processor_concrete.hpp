#pragma once

#include <cxx.hpp>
#include <concrete-ntt.hpp>

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>

class NTT_Processor_Concrete {
private:
    std::unordered_map<uint64_t, size_t> planmap;
}