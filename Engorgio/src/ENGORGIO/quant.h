#include "openfhe.h"
#include <chrono>
#include "math/chebyshev.h"
#include <stdexcept>

int extractBits(double num, int start, int length);

std::vector<double> splitNumber(double x, int w, int Bg);

std::vector<double> extractRealParts(const std::vector<std::complex<double>> &enc);

std::vector<std::vector<double>> quantization(std::vector<double> vec, uint32_t precision, uint32_t Bg);
// num_slots: ct_slots
std::vector<std::vector<double>> split_slots(std::vector<double> vec,
                                             uint32_t num_slots);

int combineNumber(const std::vector<double> &parts, int w, int Bg);
// precision: original precision of ptxt
// Bg: quant precision
std::vector<double> merge_quant(std::vector<std::vector<double>> vec,
                                uint32_t precision, uint32_t Bg);

int test_split_merge();

int test_quant(uint32_t num_input, uint32_t num_slots, uint32_t plain_precision, uint32_t plain_bits);
