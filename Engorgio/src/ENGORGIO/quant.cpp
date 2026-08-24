#include "utils.h"
#include "sort.h"
#include <chrono>
#include "math/chebyshev.h"
#include <stdexcept>

int extractBits(double num, int start, int length)
{
    return (int(num) >> start) & ((1 << length) - 1);
}

std::vector<double> splitNumber(double x, int w, int Bg)
{
    std::vector<double> result;
    for (int i = 0; i < w; ++i)
    {
        int start_bit = i * Bg;
        int extracted = extractBits(x, start_bit, Bg);
        result.push_back(extracted);
    }
    return result;
}

std::vector<double> extractRealParts(const std::vector<std::complex<double>> &enc)
{
    std::vector<double> realParts;
    for (const auto &elem : enc)
    {
        realParts.push_back(elem.real());
    }
    return realParts;
}

std::vector<std::vector<double>> quantization(std::vector<double> vec, uint32_t precision, uint32_t Bg)
{
    int w = precision / Bg;
    int length = vec.size();
    std::vector<std::vector<double>>
        quant_plain(w, std::vector<double>(length));
    for (int i = 0; i < length; i++)
    {
        auto split_vi = splitNumber(vec[i], w, Bg);
        // quant_plain[0][i] = split_vi[0];
        // quant_plain[1][i] = split_vi[1];
        for (int j = 0; j < w; j++) {
            quant_plain[j][i] = split_vi[j];
        }

    }
    return quant_plain;
}

std::vector<std::vector<double>> split_slots(std::vector<double> vec,
                                             uint32_t num_slots)
{
    int length = vec.size();
    int ct_block_size = length / num_slots;
    std::vector<std::vector<double>>
        splited_plain(ct_block_size, std::vector<double>(num_slots));
    for (int i = 0; i < ct_block_size; i++)
    {
        // std::cout << "vec: " << vec[i * num_slots] << std::endl;
        splited_plain[i].assign(vec.begin() + i * num_slots, vec.begin() + (i + 1) * num_slots);
        // std::cout << "splited_plain: " << splited_plain[i][0] << std::endl;
    }
    return splited_plain;
}

int combineNumber(const std::vector<double> &parts, int w, int Bg)
{
    int result = 0;
    for (int i = 0; i < w; ++i)
    {
        int shift = (w - 1 - i) * Bg;
        result |= (int(parts[i]) << shift);
    }
    return result;
}
// precision: original precision of ptxt
// Bg: quant precision
std::vector<double> merge_quant(std::vector<std::vector<double>> vec,
                                uint32_t precision, uint32_t Bg)
{
    int w = precision / Bg;
    int length = vec[0].size();
    std::vector<double>
        merged_plain(length);
    for (int i = 0; i < length; i++)
    {
        std::vector<double> temp_v = {vec[0][i], vec[1][i]};
        merged_plain[i] = combineNumber(temp_v, w, Bg);
    }
    return merged_plain;
}

int test_split_merge()
{
    int x = 65535;
    int w = 2;
    int Bg = 8;

    std::vector<double> parts = splitNumber(x, w, Bg);
    std::cout << "splitNumber result: ";
    for (double part : parts)
    {
        std::cout << part << " ";
    }
    std::cout << std::endl;

    int combined = combineNumber(parts, w, Bg);
    std::cout << "combineNumber result: " << combined << std::endl;

    return 0;
}

int test_quant(uint32_t num_input, uint32_t num_slots, uint32_t plain_precision, uint32_t plain_bits)
{
    std::cout << "test_quant" << std::endl;

    // We need to enable Advanced SHE to use the Chebyshev approximation.

    double precision = (1 << (plain_bits - 1)) - 1;
    std::cout << "precision: " << precision << std::endl;

    double lowerBound = -precision - 10;
    double upperBound = precision + 10;
    // double bound      = 3;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, precision);
    std::vector<double> input(num_input);
    for (int i = 0; i < num_input; i++)
    {
        input[i] = 65535 - i;
    }

    std::vector<std::vector<double>> split_input = split_slots(input, num_slots);
    int num_ct = split_input.size();
    std::cout
        << "org size: " << num_input << ", num_solts: " << num_slots << ", split_input size: " << num_ct << std::endl;
    std::vector<std::vector<std::vector<double>>> quant_input(num_ct);
    for (int i = 0; i < num_ct; i++)
    {
        quant_input[i] = quantization(split_input[i], plain_precision, plain_bits - 1);
    }

    std::vector<std::vector<double>> merge_input(num_ct);
    for (int i = 0; i < quant_input.size(); i++)
    {
        merge_input[i] = merge_quant(quant_input[i], plain_precision, plain_bits - 1);
    }
    for (int i = 0; i < input.size(); i++)
    {
        std::cout << "input:" << input[i] << ", quant_input: (" << quant_input[0][0][i] << "," << quant_input[0][1][i] << "), merge_input: " << merge_input[0][i] << std::endl;
    }

    return 0;
}
