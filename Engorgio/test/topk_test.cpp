#include "ENGORGIO/comp.h"
#include "ENGORGIO/sort.h"
#include "ENGORGIO/quant.h"
#include "ENGORGIO/utils.h"
#include <stdexcept>
using namespace openfhe;
using namespace lbcrypto;
void Eval_Topk_256(std::uint32_t plain_bits, std::uint32_t num)
{

    std::ifstream topkfile("../../engorgio_topk.csv", std::ios::app);
    std::vector<std::vector<std::string>> csvData;
    std::string line;
    while (std::getline(topkfile, line))
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
    topkfile.close();

    std::vector<double> topk_time;
    std::cout << "Eval Top-1/256" << std::endl;
    topk_time.push_back(topk_sort_test(plain_bits, num, 1));
    std::cout << "---------------------------------------------" << std::endl
              << std::endl;
    std::cout << "Eval Top-4/256" << std::endl;
    topk_time.push_back(topk_sort_test(plain_bits, num, 4));
    std::cout << "---------------------------------------------" << std::endl
              << std::endl;
    std::cout << "Eval Top-8/256" << std::endl;
    topk_time.push_back(topk_sort_test(plain_bits, num, 8));
    std::cout << "---------------------------------------------" << std::endl
              << std::endl;
    std::cout << "Eval Top-16/256" << std::endl;
    topk_time.push_back(topk_sort_test(plain_bits, num, 16));

    for (size_t i = 1; i < csvData.size(); i++)
    {
        csvData[i][3] = std::to_string(topk_time[i - 1]);
    }
    std::ofstream outFile_topk("../../engorgio_topk.csv");
    for (const auto &row : csvData)
    {
        for (size_t i = 0; i < row.size(); i++)
        {
            outFile_topk << row[i];
            if (i < row.size() - 1)
            {
                outFile_topk << ",";
            }
        }
        outFile_topk << std::endl;
    }
    outFile_topk.close();
}
int main(int argc, char *argv[])
{
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    Eval_Topk_256(8, 256);
    end = std::chrono::system_clock::now();
    double total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
    std::cout << "Topk test compute time:" << total_time << std::endl;
    return 0;
}