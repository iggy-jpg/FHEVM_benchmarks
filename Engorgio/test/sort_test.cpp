#include "ENGORGIO/comp.h"
#include "ENGORGIO/sort.h"
#include "ENGORGIO/quant.h"
#include "ENGORGIO/utils.h"
#include <stdexcept>
using namespace openfhe;
using namespace lbcrypto;

int main(int argc, char *argv[])
{
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    // std::ifstream Sortfile("../../engorgio_sort.csv", std::ios::app);
    // std::vector<std::vector<std::string>> csvData;
    // std::string line;
    // while (std::getline(Sortfile, line))
    // {
    //     std::vector<std::string> row;
    //     std::stringstream ss(line);
    //     std::string cell;
    //     while (std::getline(ss, cell, ','))
    //     {
    //         row.push_back(cell);
    //     }
    //     csvData.push_back(row);
    // }
    // Sortfile.close();
    // std::vector<double> sort_time;
    // sort_time = bitonic_sort_modular(8, 2, 32768);
    bitonic_sort(8, 4);

    // sort_time.push_back(bitonic_sort_full_table_scan(8, 65536 * 2) + sort_time[15] * 2);

    // for (size_t i = 1; i < csvData.size(); i++)
    // {
    //     csvData[i][5] = std::to_string(sort_time[i + 1]);
    // }
    // std::ofstream outFile_sort("../../engorgio_sort.csv");
    // for (const auto &row : csvData)
    // {
    //     for (size_t i = 0; i < row.size(); i++)
    //     {
    //         outFile_sort << row[i];
    //         if (i < row.size() - 1)
    //         {
    //             outFile_sort << ",";
    //         }
    //     }
    //     outFile_sort << std::endl;
    // }
    // outFile_sort.close();
    end = std::chrono::system_clock::now();
    double total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
    std::cout << "Sort test compute time:" << total_time << std::endl;
    return 0;
}