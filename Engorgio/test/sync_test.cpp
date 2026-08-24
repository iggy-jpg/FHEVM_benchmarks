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
    std::ifstream Syncfile("../../engorgio_sync.csv", std::ios::app);
    std::vector<std::vector<std::string>> csvData;
    std::string line;
    while (std::getline(Syncfile, line))
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
    Syncfile.close();
    std::vector<double> sync_time;
    for (int i = 8; i < 8193; i *= 2)
    {
        sync_time.push_back(sync_test_small(8, i));
    }
    for (size_t i = 1; i < csvData.size(); i++)
    {
        csvData[i][3] = std::to_string(sync_time[i - 1]);
    }
    std::ofstream outFile_sync("../../engorgio_sync.csv");
    for (const auto &row : csvData)
    {
        for (size_t i = 0; i < row.size(); i++)
        {
            outFile_sync << row[i];
            if (i < row.size() - 1)
            {
                outFile_sync << ",";
            }
        }
        outFile_sync << std::endl;
    }
    outFile_sync.close();
    end = std::chrono::system_clock::now();
    double total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
    std::cout << "Sync test compute time:" << total_time << std::endl;
    return 0;
}