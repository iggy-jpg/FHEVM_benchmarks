#include "ENGORGIO/comp.h"
#include "ENGORGIO/sort.h"
#include "ENGORGIO/utils.h"
#include <chrono>
#include <omp.h>
#include "math/chebyshev.h"
using namespace openfhe;
using namespace lbcrypto;

void Eval_hybrid_query_1()
{
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits = 78;
    usint firstMod = 89;
#else
    usint dcrtBits = 59;
    usint firstMod = 60;
#endif
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetFirstModSize(firstMod);
    std::vector<std::uint32_t> levelBudget = {4, 4};
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    std::uint32_t polyDegree = 495;
    std::uint32_t numIterations = 2;
    std::uint32_t levelsAvailableAfterBootstrap = 25;
    usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    // std::cout << "RingDimension " << cc->GetRingDimension() << std::endl;
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(FHE);
    cc->Enable(ADVANCEDSHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.
    int plain_bits = 8;
    int block = 2;
    double precision = (1 << (plain_bits - 1)) - 1;

    double lowerBound = -precision - 10;
    double upperBound = precision + 10;
    // double bound      = 3;
    auto keyPair = cc->KeyGen();
    usint ringDim = cc->GetRingDimension();
    int length = ringDim / 2;

    std::ifstream csvFile("../../engorgio_hybrid_hq1_2.csv", std::ios::app);
    std::vector<std::vector<std::string>> csvData;
    std::string line;
    while (std::getline(csvFile, line))
    {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ','))
        {
            row.push_back(cell);
        }
        csvData.push_back(row);
    }
    csvFile.close();
    auto sort_time_table = bitonic_sort_modular_query(8, 2, 4096);
    for (int num = 64; num < 4097; num *= 4)
    {
        std::cout << "--------------- Eval HQ1 with " << num << "records ---------------" << std::endl;
        // filter
        double filter_time = 0.;
        // 8*num comp
        std::cout << "Eval filtering..." << std::endl;
        for (int i = 0; i < 4; i++)
            filter_time += Eval_modular_Strictly_greater_test(plain_bits, block);
        for (int i = 0; i < 4; i++)
            filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
        filter_time *= num;
        std::cout << num << "filtering time: " << filter_time << "ms" << std::endl;
        double sort_time = sort_time_table[std::log2(num) - 1];
        std::cout
            << "Eval data ordering..." << std::endl;
        std::cout << num << "data ordering time: " << sort_time << "ms" << std::endl;

        std::cout << "Eval agg..." << std::endl;
        double agg_time = Eval_Agg_1(plain_bits, num);
        std::cout << num << "agg time: " << agg_time << "ms" << std::endl;

        std::cout << "Eval sync..." << std::endl;
        double sync_time = sync_test_small(plain_bits, num);
        std::cout << num << "sync time: " << sync_time << "ms" << std::endl;
        double total_time = filter_time + agg_time + sort_time + sync_time;
        std::cout << num << "HQ1 total_time time: " << total_time << "ms" << std::endl
                  << std::endl
                  << std::endl;
        for (size_t i = 1; i < csvData.size(); i++)
        {
            if (std::stoi(csvData[i][0]) == num)
            {
                csvData[i][4] = std::to_string(total_time);
                break;
            }
        }
    }
    std::ofstream outFile("../../engorgio_hybrid_hq1_2.csv");
    for (const auto &row : csvData)
    {
        for (size_t i = 0; i < row.size(); i++)
        {
            outFile << row[i];
            if (i < row.size() - 1)
            {
                outFile << ",";
            }
        }
        outFile << std::endl;
    }
    outFile.close();
}
void Eval_hybrid_query_2()
{
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits = 78;
    usint firstMod = 89;
#else
    usint dcrtBits = 59;
    usint firstMod = 60;
#endif
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetFirstModSize(firstMod);
    std::vector<std::uint32_t> levelBudget = {4, 4};
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    std::uint32_t polyDegree = 495;
    std::uint32_t numIterations = 2;
    std::uint32_t levelsAvailableAfterBootstrap = 25;
    usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    // std::cout << "RingDimension " << cc->GetRingDimension() << std::endl;
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(FHE);
    cc->Enable(ADVANCEDSHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.
    int plain_bits = 8;
    int block = 2;
    double precision = (1 << (plain_bits - 1)) - 1;

    double lowerBound = -precision - 10;
    double upperBound = precision + 10;
    // double bound      = 3;
    auto keyPair = cc->KeyGen();
    usint ringDim = cc->GetRingDimension();
    int length = ringDim / 2;
    int k = 4;

    std::ifstream csvFile("../../engorgio_hybrid_hq2_2.csv", std::ios::app);
    std::vector<std::vector<std::string>> csvData;
    std::string line;
    while (std::getline(csvFile, line))
    {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ','))
        {
            row.push_back(cell);
        }
        csvData.push_back(row);
    }
    csvFile.close();

    for (int num = 64; num < 4097; num *= 4)
    { // filter
        std::cout << "--------------- Eval HQ2 with " << num << "records ---------------" << std::endl;
        double filter_time = 0.;
        // 8*num comp
        std::cout << "Eval filtering..." << std::endl;
        for (int i = 0; i < 4; i++)
            filter_time += Eval_modular_Strictly_greater_test(plain_bits, block);
        for (int i = 0; i < 4; i++)
            filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
        filter_time *= num;
        std::cout << num << "filtering time: " << filter_time << "ms" << std::endl;

        std::cout << "Eval data ordering..." << std::endl;
        std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> sync_matrix;
        auto topk_time = topk_sort_test(plain_bits, num, k);

        std::cout << num << " data ordering time: " << topk_time << "ms" << std::endl;

        std::cout << "Eval agg..." << std::endl;
        double agg_time = Eval_Agg_1(plain_bits, num);
        std::cout << num << "agg_time time: " << agg_time << "ms" << std::endl;
        std::cout << "Eval sync..." << std::endl;
        double sync_time = sync_test_small(plain_bits, num);
        std::cout << num << "sync time: " << sync_time << "ms" << std::endl;
        double total_time = filter_time + agg_time + topk_time + sync_time;
        std::cout << num << "HQ2 total_time time: " << total_time << "ms" << std::endl
                  << std::endl
                  << std::endl;
        for (size_t i = 1; i < csvData.size(); i++)
        {
            if (std::stoi(csvData[i][0]) == num)
            {
                csvData[i][4] = std::to_string(total_time);
                break;
            }
        }
    }
    std::ofstream outFile("../../engorgio_hybrid_hq2_2.csv");
    for (const auto &row : csvData)
    {
        for (size_t i = 0; i < row.size(); i++)
        {
            outFile << row[i];
            if (i < row.size() - 1)
            {
                outFile << ",";
            }
        }
        outFile << std::endl;
    }
    outFile.close();
}

int main(int argc, char *argv[])
{
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    std::cout << "--------------------- Eval HQ1 --------------------- " << std::endl;
    Eval_hybrid_query_1();
    end = std::chrono::system_clock::now();
    double total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
    std::cout << "HQ1 compute time:" << total_time << std::endl;

    return 0;
}