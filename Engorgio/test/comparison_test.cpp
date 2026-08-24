#include "ENGORGIO/comp.h"
#include "ENGORGIO/sort.h"
#include "ENGORGIO/quant.h"
#include <chrono>
using namespace openfhe;
using namespace lbcrypto;

void Eval_modular_Strictly_greater_test(std::uint32_t Bg)
{
    std::ifstream csvFile("../../engorgio_ge.csv", std::ios::app);
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
    std::cout << " --------------------- Greater (relational) comparison 8-64 bit --------------------- " << std::endl;

    double comp_time8 = Eval_modular_Strictly_greater_than(Bg, 1);
    std::cout << std::endl;
    double comp_time16 = Eval_modular_Strictly_greater_than(Bg, 2);
    std::cout << std::endl;
    double comp_time20 = Eval_modular_Strictly_greater_than(Bg, 3);
    std::cout << std::endl;
    double comp_time32 = Eval_modular_Strictly_greater_than(Bg, 4);
    std::cout << std::endl;
    double comp_time64 = Eval_modular_Strictly_greater_than(Bg, 8);
    std::cout << std::endl;
    for (size_t i = 1; i < csvData.size(); i++)
    {
        if (std::stoi(csvData[i][0]) == 8)
        {
            csvData[i][8] = std::to_string(comp_time8);
        }
        else if (std::stoi(csvData[i][0]) == 16)
        {
            csvData[i][8] = std::to_string(comp_time16);
        }
        else if (std::stoi(csvData[i][0]) == 20)
        {
            csvData[i][8] = std::to_string(comp_time20);
        }
        else if (std::stoi(csvData[i][0]) == 32)
        {
            csvData[i][8] = std::to_string(comp_time32);
        }
        else if (std::stoi(csvData[i][0]) == 64)
        {
            csvData[i][8] = std::to_string(comp_time64);
        }
    }
    std::ofstream outFile("../../engorgio_ge.csv");
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

void Eval_modular_Strictly_equal_test(std::uint32_t Bg)
{
    std::ifstream csvFile("../../engorgio_eq.csv", std::ios::app);
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
    std::cout << " --------------------- Equality comparison 8-64 bit --------------------- " << std::endl;
    double comp_time8 = Eval_modular_Strictly_equal(Bg, 1);
    std::cout << std::endl;
    double comp_time16 = Eval_modular_Strictly_equal(Bg, 2);
    std::cout << std::endl;
    double comp_time20 = Eval_modular_Strictly_equal(Bg, 3);
    std::cout << std::endl;
    double comp_time32 = Eval_modular_Strictly_equal(Bg, 4);
    std::cout << std::endl;
    double comp_time64 = Eval_modular_Strictly_equal(Bg, 8);
    std::cout << std::endl;
    for (size_t i = 1; i < csvData.size(); i++)
    {
        if (std::stoi(csvData[i][0]) == 8)
        {
            csvData[i][8] = std::to_string(comp_time8);
        }
        else if (std::stoi(csvData[i][0]) == 16)
        {
            csvData[i][8] = std::to_string(comp_time16);
        }
        else if (std::stoi(csvData[i][0]) == 20)
        {
            csvData[i][8] = std::to_string(comp_time20);
        }
        else if (std::stoi(csvData[i][0]) == 32)
        {
            csvData[i][8] = std::to_string(comp_time32);
        }
        else if (std::stoi(csvData[i][0]) == 64)
        {
            csvData[i][8] = std::to_string(comp_time64);
        }
    }
    std::ofstream outFile("../../engorgio_eq.csv");
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
    Eval_modular_Strictly_greater_test(8);
    Eval_modular_Strictly_equal_test(8);
    end = std::chrono::system_clock::now();
    double total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
    std::cout << "Comparison test compute time:" << total_time << std::endl;
    return 0;
}
