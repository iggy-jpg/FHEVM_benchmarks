#include "ENGORGIO/comp.h"
#include "ENGORGIO/sort.h"
#include "ENGORGIO/utils.h"
#include <chrono>
#include <omp.h>
#include "math/chebyshev.h"
using namespace openfhe;
using namespace lbcrypto;
void Eval_VQ1()
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

    std::ifstream vq_file("../../engorgio_top4_query.csv", std::ios::app);
    std::vector<std::vector<std::string>> csvData_vq;
    std::string line;
    while (std::getline(vq_file, line))
    {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ','))
        {
            row.push_back(cell);
        }
        csvData_vq.push_back(row);
    }
    vq_file.close();
    std::vector<double> vq_time;
    for (int num_slots = 128; num_slots < 1025; num_slots *= 2)
    {
        std::cout << "---------- Eval VQ1 with " << num_slots << " records ----------" << std::endl;

        // dist
        std::cout << "Eval Homdis ..." << std::endl;
        double dis_time = 0.;
        dis_time = Euclid_distance(plain_bits) * block * num_slots / 65536;
        std::cout << num_slots << " dist time: " << dis_time << " ms" << std::endl;
        // knn
        std::cout << "Eval Homknn ..." << std::endl;
        std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> sync_matrix;
        auto topk_time = topk_sort_test(plain_bits, num_slots, k);
        std::cout << num_slots << " knn time: " << topk_time << " ms" << std::endl;
        int sync_col_num = 9;
        double sync_time = 0.;
        std::cout << "Eval Homsync ..." << std::endl;
        sync_time = sync_test_small(plain_bits, num_slots) * sync_col_num;
        std::cout << num_slots << " sync time: " << sync_time << " ms" << std::endl;
        double total_time = sync_time + topk_time;
        std::cout << "VQ1 with " << num_slots << " records total time: " << total_time << "ms" << std::endl
                  << std::endl
                  << std::endl;
        vq_time.push_back(total_time);
    }
    for (size_t i = 1; i < csvData_vq.size(); i++)
    {
        csvData_vq[i][3] = std::to_string(vq_time[i - 1]);
    }
    std::ofstream outFile_vq("../../engorgio_top4_query.csv");
    for (const auto &row : csvData_vq)
    {
        for (size_t i = 0; i < row.size(); i++)
        {
            outFile_vq << row[i];
            if (i < row.size() - 1)
            {
                outFile_vq << ",";
            }
        }
        outFile_vq << std::endl;
    }
    outFile_vq.close();
}

void Eval_VQ2()
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
    int k = 16;

    std::ifstream vq_file("../../engorgio_top16_query.csv", std::ios::app);
    std::vector<std::vector<std::string>> csvData_vq;
    std::string line;
    while (std::getline(vq_file, line))
    {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ','))
        {
            row.push_back(cell);
        }
        csvData_vq.push_back(row);
    }
    vq_file.close();
    std::vector<double> vq_time;
    for (int num_slots = 128; num_slots < 1025; num_slots *= 2)
    {
        std::cout << "---------- Eval VQ2 with " << num_slots << " records ----------" << std::endl;

        // dist
        std::cout << "Eval Homdis " << std::endl;
        double dis_time = 0.;
        dis_time = Euclid_distance(plain_bits) * block * num_slots / 65536;
        std::cout << num_slots << " distance time: " << dis_time << " ms" << std::endl;
        // knn
        std::cout << "Eval Homknn " << std::endl;
        std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> sync_matrix;
        auto topk_time = topk_sort_test(plain_bits, num_slots, k);
        std::cout << num_slots << " knn time: " << topk_time << " ms" << std::endl;
        int sync_col_num = 5;
        double sync_time = 0.;
        std::cout << "Eval HomSync " << std::endl;
        sync_time += sync_test_small(plain_bits, num_slots) * sync_col_num;
        std::cout << num_slots << " sync time: " << sync_time << " ms" << std::endl;
        double total_time = sync_time + topk_time;
        std::cout << num_slots << "VQ2 with" << num_slots << " records total time: " << total_time << "ms" << std::endl;
        vq_time.push_back(total_time);
    }
    for (size_t i = 1; i < csvData_vq.size(); i++)
    {
        csvData_vq[i][3] = std::to_string(vq_time[i - 1]);
    }
    std::ofstream outFile_vq("../../engorgio_top16_query.csv");
    for (const auto &row : csvData_vq)
    {
        for (size_t i = 0; i < row.size(); i++)
        {
            outFile_vq << row[i];
            if (i < row.size() - 1)
            {
                outFile_vq << ",";
            }
        }
        outFile_vq << std::endl;
    }
    outFile_vq.close();
}

int main(int argc, char *argv[])
{
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    Eval_VQ1();
    // Eval_VQ2();
    end = std::chrono::system_clock::now();
    double total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "VQ1 test compute time:" << total_time << std::endl;
    return 0;
}