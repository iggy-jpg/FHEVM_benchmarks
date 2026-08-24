#pragma once
#include <vector>
#include <random>
#include <type_traits>
#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include "tfhe++.hpp" 
#include <thread> 
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

inline boost::multiprecision::cpp_int bits_to_int_signed(const std::vector<uint8_t>& v, size_t W) {
    using boost::multiprecision::cpp_int;

    cpp_int val = 0;
    for (size_t i = 0; i < W; ++i) {
        val = (val << 1) | (v[i] & 1);
    }
    if (v[0] & 1) {
        val -= (cpp_int(1) << W);
    }

    return val;
}

inline void print_bits(const std::string& tag, const std::vector<uint8_t>& bits) {
    std::cout << tag;
    for (uint8_t b : bits) std::cout << int(b) << ' ';
    std::cout << '\n';
}

template<class P>
static std::array<bool, P::n> to_bit_array(const std::vector<uint8_t>& v) {
    std::array<bool, P::n> out{};
    const size_t m = std::min(v.size(), static_cast<size_t>(P::n));
    for (size_t i = 0; i < m; ++i) out[i] = (v[i] & 1u) != 0;
    return out;
}

inline int pick_thread_count_hint(int hint) {
#ifdef _OPENMP
    const int runtime_max = std::max(1, std::min(omp_get_max_threads(),
                                                 omp_get_thread_limit()));
    if (hint > 0) return std::min(hint, runtime_max);
    return runtime_max;
#else
    if (hint > 0) return hint;
    unsigned hc = std::thread::hardware_concurrency();
    if (hc == 0) return 1;
    return std::max(1u, hc / 2);
#endif
}



//







//


#include <cstdlib>
inline bool wide_parallel_enabled(int threads_hint, int min_threads = 4) {
    static const int env_threshold = [] {
        const char* env = std::getenv("AUTOHOG_WIDE_THRESHOLD");
        return env ? std::max(1, atoi(env)) : 0;
    }();
    const int threshold = env_threshold > 0 ? env_threshold : min_threads;
    return pick_thread_count_hint(threads_hint) >= threshold;
}
