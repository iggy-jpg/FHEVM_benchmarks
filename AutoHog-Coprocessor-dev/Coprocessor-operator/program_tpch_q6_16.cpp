// TPC-H Query 6 over encrypted 16-bit columns.
// Usage: program_tpch_q6_16 [row-count]
#include "HomOP.h"
#include "workload_common.hpp"

#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace TFHEpp;
using namespace autohog_workload;

struct PlainRow {
    uint16_t extended_price; // fixed-point monetary units
    uint16_t discount;       // basis points; 600 means 6.00%
    uint16_t ship_date;      // days since an application-defined epoch
    uint16_t quantity;
};

static TLWE<lvl2param> and_bit(
    const TLWE<lvl2param> &a, const TLWE<lvl2param> &b,
    const EvalKey &ek, int threads)
{
    return HomAND({a}, {b}, 1, ek, threads)[0];
}

int main(int argc, char **argv)
{
    try {
        constexpr int field_width = 16;
        constexpr int revenue_width = 32;
        const int row_count = argc > 1 ? std::stoi(argv[1]) : 2;
        if (row_count <= 0) throw std::runtime_error("row-count must be positive");
        const int threads = thread_count();

        // Q6 parameters: one-year date window, discount +/- 1 percentage point.
        const uint16_t date = 10000;
        const uint16_t date_next_year = date + 365;
        const uint16_t target_discount = 600;
        const uint16_t discount_low = target_discount - 100;
        const uint16_t discount_high = target_discount + 100;
        const uint16_t quantity_limit = 24;

        std::vector<PlainRow> rows;
        rows.reserve(row_count);
        for (int i = 0; i < row_count; ++i) {
            if ((i & 1) == 0)
                rows.push_back({static_cast<uint16_t>(1000 + i), 600, 10100, 10});
            else
                rows.push_back({static_cast<uint16_t>(2000 + i), 900, 10100, 30});
        }

        SecretKey sk;
        EvalKey ek;
        ek.emplacebkfft<lvl02param>(sk);
        ek.emplaceiksk<lvl20param>(sk);
        init_constants(sk);

        const auto enc_date = encrypt_uint(date, field_width, sk);
        const auto enc_next_year = encrypt_uint(date_next_year, field_width, sk);
        const auto enc_discount_low = encrypt_uint(discount_low, field_width, sk);
        const auto enc_discount_high = encrypt_uint(discount_high, field_width, sk);
        const auto enc_quantity_limit = encrypt_uint(quantity_limit, field_width, sk);

        struct EncryptedRow {
            std::vector<TLWE<lvl2param>> price, discount, ship_date, quantity;
        };
        std::vector<EncryptedRow> encrypted_rows;
        encrypted_rows.reserve(rows.size());
        for (const auto &row : rows) {
            encrypted_rows.push_back({
                encrypt_uint(row.extended_price, field_width, sk),
                encrypt_uint(row.discount, field_width, sk),
                encrypt_uint(row.ship_date, field_width, sk),
                encrypt_uint(row.quantity, field_width, sk)});
        }

        auto revenue = encrypted_zero(revenue_width);
        const auto zero_revenue = encrypted_zero(revenue_width);
        const auto begin = std::chrono::steady_clock::now();
        for (const auto &row : encrypted_rows) {
            auto selected = HomGE_unsigned_pal(
                row.ship_date, enc_date, field_width, ek, threads);
            selected = and_bit(selected, HomLT_unsigned_pal(
                row.ship_date, enc_next_year, field_width, ek, threads), ek, threads);
            selected = and_bit(selected, HomGE_unsigned_pal(
                row.discount, enc_discount_low, field_width, ek, threads), ek, threads);
            selected = and_bit(selected, HomLE_unsigned_pal(
                row.discount, enc_discount_high, field_width, ek, threads), ek, threads);
            selected = and_bit(selected, HomLT_unsigned_pal(
                row.quantity, enc_quantity_limit, field_width, ek, threads), ek, threads);

            // 16x16 produces the natural 32-bit fixed-point revenue term.
            const auto term = HomMul_pal(
                row.price, row.discount, field_width, ek, threads);
            const auto filtered_term = HomSelect(
                zero_revenue, term, selected, revenue_width, ek, threads);
            revenue = HomAdd_unsigned_pal(
                revenue, filtered_term, revenue_width, ek, false, threads).SUM;
        }
        const auto end = std::chrono::steady_clock::now();

        uint64_t expected = 0;
        for (const auto &row : rows) {
            const bool match = row.ship_date >= date && row.ship_date < date_next_year &&
                row.discount >= discount_low && row.discount <= discount_high &&
                row.quantity < quantity_limit;
            if (match) expected += uint64_t(row.extended_price) * row.discount;
        }
        expected &= 0xffffffffULL;
        const cpp_int actual = decrypt_uint(revenue, sk);
        const double elapsed_ms =
            std::chrono::duration<double, std::milli>(end - begin).count();
        const bool correct = actual == expected;
        std::cout << std::fixed << std::setprecision(3)
                  << "TPCH_Q6_16 rows=" << row_count
                  << " revenue=" << actual
                  << " expected=" << expected
                  << " ms=" << elapsed_ms
                  << " correct=" << (correct ? "yes" : "no") << '\n';
        return correct ? 0 : 2;
    } catch (const std::exception &error) {
        std::cerr << "program_tpch_q6_16: " << error.what() << '\n';
        return 1;
    }
}
