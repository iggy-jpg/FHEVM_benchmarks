#include "Autohog_Cop.hpp"
#include "HomOP.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <ctime>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace TFHEpp;

namespace {

struct Config {
    int width = 64;
    std::string op = "add";
    std::string variant = "parallel";
    std::string input_mode = "trivial";
    int op_threads = 1;
    int dag_parallel = 1;
    int warmup = 1;
    int repeat = 10;
    std::string csv_path;
    bool quiet = false;
};

struct CipherPair {
    std::vector<TLWE<lvl2param>> a;
    std::vector<TLWE<lvl2param>> b;
    TLWE<lvl2param> sel;
};

struct WaveResult {
    double ms = 0.0;
    std::uint64_t checksum = 0;
};

[[noreturn]] void usage(const char *prog) {
    std::cerr
        << "Usage: " << prog << " [options]\n"
        << "Options:\n"
        << "  --width N             Operand width, usually 64 or 256\n"
        << "  --op NAME             add, sub, add_unsigned, sub_unsigned, and, or, xor,\n"
        << "                        select, gt, ge, lt, le, eq, ne, gt_unsigned,\n"
        << "                        ge_unsigned, lt_unsigned, le_unsigned\n"
        << "  --variant NAME        serial or parallel. Default: parallel\n"
        << "  --op-threads N        OpenMP threads used inside one operator\n"
        << "  --dag-parallel N      Independent operators launched at the DAG layer\n"
        << "  --warmup N            Warmup repeats excluded from timing\n"
        << "  --repeat N            Timed repeats (default: 10)\n"
        << "  --input-mode NAME     trivial or std. std uses fresh TLWE encryption. Default: trivial\n"
        << "  --csv PATH            Append one summary row to PATH\n"
        << "  --quiet               Print only the CSV summary line\n";
    std::exit(1);
}

int parse_positive(const char *value, const char *name) {
    char *end = nullptr;
    long v = std::strtol(value, &end, 10);
    if (end == value || *end != '\0' || v <= 0 || v > 1000000) {
        throw std::invalid_argument(std::string("invalid positive integer for ") + name);
    }
    return static_cast<int>(v);
}

int parse_nonnegative(const char *value, const char *name) {
    char *end = nullptr;
    long v = std::strtol(value, &end, 10);
    if (end == value || *end != '\0' || v < 0 || v > 1000000) {
        throw std::invalid_argument(std::string("invalid nonnegative integer for ") + name);
    }
    return static_cast<int>(v);
}

Config parse_args(int argc, char **argv) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto need_value = [&](const char *name) -> const char * {
            if (i + 1 >= argc) usage(argv[0]);
            ++i;
            return argv[i];
        };

        if (arg == "--width") cfg.width = parse_positive(need_value("--width"), "--width");
        else if (arg == "--op") cfg.op = need_value("--op");
        else if (arg == "--variant") cfg.variant = need_value("--variant");
        else if (arg == "--op-threads") cfg.op_threads = parse_positive(need_value("--op-threads"), "--op-threads");
        else if (arg == "--dag-parallel") cfg.dag_parallel = parse_positive(need_value("--dag-parallel"), "--dag-parallel");
        else if (arg == "--warmup") cfg.warmup = parse_nonnegative(need_value("--warmup"), "--warmup");
        else if (arg == "--repeat" || arg == "--waves") cfg.repeat = parse_positive(need_value("--repeat"), "--repeat");
        else if (arg == "--input-mode") cfg.input_mode = need_value("--input-mode");
        else if (arg == "--csv") cfg.csv_path = need_value("--csv");
        else if (arg == "--quiet") cfg.quiet = true;
        else if (arg == "--help" || arg == "-h") usage(argv[0]);
        else usage(argv[0]);
    }

    if (cfg.width <= 0) throw std::invalid_argument("width must be positive");
    if (cfg.variant != "serial" && cfg.variant != "parallel") {
        throw std::invalid_argument("variant must be serial or parallel");
    }
    if (cfg.input_mode == "secret") cfg.input_mode = "std";
    if (cfg.input_mode != "trivial" && cfg.input_mode != "std") {
        throw std::invalid_argument("input-mode must be trivial or std");
    }
    return cfg;
}

TLWE<lvl2param> make_bit(uint8_t bit, const SecretKey &sk, const std::string &input_mode) {
    if (input_mode == "trivial") {
        TLWE<lvl2param> ct{};
        ct[lvl2param::n] = bit ? lvl2param::μ : -lvl2param::μ;
        return ct;
    }
    return tlweSymEncrypt<lvl2param>(
        bit ? lvl2param::μ : -lvl2param::μ,
        lvl2param::α,
        sk.key.get<lvl2param>());
}

std::vector<CipherPair> make_dataset(const Config &cfg, const SecretKey &sk) {
    std::vector<CipherPair> dataset;
    dataset.reserve(cfg.dag_parallel);

    std::mt19937_64 rng(0x20260508ULL);
    std::uniform_int_distribution<int> bit01(0, 1);

    for (int item = 0; item < cfg.dag_parallel; ++item) {
        CipherPair pair;
        pair.a.resize(cfg.width);
        pair.b.resize(cfg.width);
        for (int i = 0; i < cfg.width; ++i) {
            pair.a[i] = make_bit(static_cast<uint8_t>(bit01(rng)), sk, cfg.input_mode);
            pair.b[i] = make_bit(static_cast<uint8_t>(bit01(rng)), sk, cfg.input_mode);
        }
        pair.sel = make_bit(static_cast<uint8_t>(bit01(rng)), sk, cfg.input_mode);
        dataset.push_back(std::move(pair));
    }
    return dataset;
}

std::uint64_t mix(std::uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

std::uint64_t checksum_tlwe(const TLWE<lvl2param> &ct) {
    return mix(static_cast<std::uint64_t>(ct[lvl2param::n]));
}

std::uint64_t checksum_vec(const std::vector<TLWE<lvl2param>> &cts) {
    std::uint64_t acc = 0x9e3779b97f4a7c15ULL;
    if (!cts.empty()) acc ^= checksum_tlwe(cts.front());
    if (cts.size() > 1) acc ^= mix(checksum_tlwe(cts.back()) + cts.size());
    return acc;
}

std::uint64_t checksum_add(const HomAddResult &res) {
    std::uint64_t acc = checksum_vec(res.SUM);
    if (res.ovf_ct.has_value()) acc ^= checksum_tlwe(*res.ovf_ct);
    return acc;
}

std::uint64_t run_one_op(const Config &cfg, const CipherPair &in, const EvalKey &ek) {
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(std::max(1, cfg.op_threads));
#endif
    const bool serial = cfg.variant == "serial";
    const int threads = serial ? 1 : cfg.op_threads;

    if (cfg.op == "add") {
        return checksum_add(serial ? HomAdd(in.a, in.b, cfg.width, ek, false)
                                   : HomAdd_pal(in.a, in.b, cfg.width, ek, false, threads));
    }
    if (cfg.op == "sub") {
        return checksum_add(serial ? HomSub(in.a, in.b, cfg.width, ek, false)
                                   : HomSub_pal(in.a, in.b, cfg.width, ek, false, threads));
    }
    if (cfg.op == "add_unsigned") {
        return checksum_add(serial ? HomAdd_unsigned(in.a, in.b, cfg.width, ek, false)
                                   : HomAdd_unsigned_pal(in.a, in.b, cfg.width, ek, false));
    }
    if (cfg.op == "sub_unsigned") {
        return checksum_add(serial ? HomSub_unsigned(in.a, in.b, cfg.width, ek, false)
                                   : HomSub_unsigned_pal(in.a, in.b, cfg.width, ek, false));
    }
    if (cfg.op == "and") return checksum_vec(HomAND(in.a, in.b, cfg.width, ek, threads));
    if (cfg.op == "or") return checksum_vec(HomOR(in.a, in.b, cfg.width, ek, threads));
    if (cfg.op == "xor") return checksum_vec(HomXOR(in.a, in.b, cfg.width, ek, threads));
    if (cfg.op == "select") return checksum_vec(HomSelect(in.a, in.b, in.sel, cfg.width, ek, threads));

    if (cfg.op == "gt") {
        return checksum_tlwe(serial ? HomGT(in.a, in.b, cfg.width, ek)
                                    : HomGT_pal(in.a, in.b, cfg.width, ek, threads));
    }
    if (cfg.op == "ge") {
        return checksum_tlwe(serial ? HomGE(in.a, in.b, cfg.width, ek)
                                    : HomGE_pal(in.a, in.b, cfg.width, ek, threads));
    }
    if (cfg.op == "lt") {
        return checksum_tlwe(serial ? HomLT(in.a, in.b, cfg.width, ek)
                                    : HomLT_pal(in.a, in.b, cfg.width, ek, threads));
    }
    if (cfg.op == "le") {
        return checksum_tlwe(serial ? HomLE(in.a, in.b, cfg.width, ek)
                                    : HomLE_pal(in.a, in.b, cfg.width, ek, threads));
    }
    if (cfg.op == "eq") {
        return checksum_tlwe(serial ? HomEQ(in.a, in.b, cfg.width, ek)
                                    : HomEQ_pal(in.a, in.b, cfg.width, ek, threads));
    }
    if (cfg.op == "ne") {
        return checksum_tlwe(serial ? HomNE(in.a, in.b, cfg.width, ek)
                                    : HomNE_pal(in.a, in.b, cfg.width, ek, threads));
    }
    if (cfg.op == "gt_unsigned") {
        return checksum_tlwe(serial ? HomGT_unsigned(in.a, in.b, cfg.width, ek)
                                    : HomGT_unsigned_pal(in.a, in.b, cfg.width, ek, threads));
    }
    if (cfg.op == "ge_unsigned") {
        return checksum_tlwe(serial ? HomGE_unsigned(in.a, in.b, cfg.width, ek)
                                    : HomGE_unsigned_pal(in.a, in.b, cfg.width, ek, threads));
    }
    if (cfg.op == "lt_unsigned") {
        return checksum_tlwe(serial ? HomLT_unsigned(in.a, in.b, cfg.width, ek)
                                    : HomLT_unsigned_pal(in.a, in.b, cfg.width, ek, threads));
    }
    if (cfg.op == "le_unsigned") {
        return checksum_tlwe(serial ? HomLE_unsigned(in.a, in.b, cfg.width, ek)
                                    : HomLE_unsigned_pal(in.a, in.b, cfg.width, ek, threads));
    }

    throw std::invalid_argument("unknown op: " + cfg.op);
}

WaveResult run_wave(const Config &cfg, const std::vector<CipherPair> &dataset, const EvalKey &ek) {
    std::vector<std::thread> workers;
    std::vector<std::uint64_t> checks(static_cast<size_t>(cfg.dag_parallel), 0);
    std::atomic<int> ready{0};
    std::atomic<bool> go{false};
    std::atomic<bool> failed{false};
    std::mutex error_mutex;
    std::string error_msg;

    workers.reserve(cfg.dag_parallel);
    for (int i = 0; i < cfg.dag_parallel; ++i) {
        workers.emplace_back([&, i]() {
            ready.fetch_add(1, std::memory_order_release);
            while (!go.load(std::memory_order_acquire)) {
                std::this_thread::yield();
            }
            try {
                checks[static_cast<size_t>(i)] =
                    run_one_op(cfg, dataset[static_cast<size_t>(i % dataset.size())], ek);
            } catch (const std::exception &e) {
                failed.store(true, std::memory_order_release);
                std::lock_guard<std::mutex> lock(error_mutex);
                error_msg = e.what();
            }
        });
    }

    while (ready.load(std::memory_order_acquire) < cfg.dag_parallel) {
        std::this_thread::yield();
    }

    const auto t0 = std::chrono::steady_clock::now();
    go.store(true, std::memory_order_release);
    for (auto &worker : workers) worker.join();
    const auto t1 = std::chrono::steady_clock::now();

    if (failed.load(std::memory_order_acquire)) {
        throw std::runtime_error(error_msg.empty() ? "worker failed" : error_msg);
    }

    std::uint64_t checksum = 0;
    for (std::uint64_t c : checks) checksum ^= mix(c + checksum);

    WaveResult result;
    result.ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    result.checksum = checksum;
    return result;
}

std::string now_compact() {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
    return oss.str();
}

bool file_exists(const std::string &path) {
    std::ifstream f(path);
    return f.good();
}

std::string csv_header() {
    return "timestamp,width,op,variant,input_mode,op_threads,dag_parallel,warmup,repeat,"
           "total_ops,total_ms,avg_repeat_ms,avg_op_ms,throughput_ops_s,checksum,"
           "hw_threads,omp_max_threads";
}

std::string make_csv_row(const Config &cfg,
                         double total_ms,
                         double avg_repeat_ms,
                         double avg_op_ms,
                         double throughput,
                         std::uint64_t checksum,
                         int hw_threads,
                         int omp_max_threads) {
    std::ostringstream row;
    row << now_compact() << ','
        << cfg.width << ','
        << cfg.op << ','
        << cfg.variant << ','
        << cfg.input_mode << ','
        << cfg.op_threads << ','
        << cfg.dag_parallel << ','
        << cfg.warmup << ','
        << cfg.repeat << ','
        << (cfg.repeat * cfg.dag_parallel) << ','
        << std::fixed << std::setprecision(3)
        << total_ms << ','
        << avg_repeat_ms << ','
        << avg_op_ms << ','
        << throughput << ','
        << checksum << ','
        << hw_threads << ','
        << omp_max_threads;
    return row.str();
}

} // namespace

int main(int argc, char **argv) {
    try {
        Config cfg = parse_args(argc, argv);

#ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(std::max(1, cfg.op_threads));
        const int omp_max_threads = omp_get_max_threads();
#else
        const int omp_max_threads = 1;
#endif

        const int hw_threads = std::max(1u, std::thread::hardware_concurrency());

        if (!cfg.quiet) {
            std::cout << "[init] width=" << cfg.width
                      << " op=" << cfg.op
                      << " variant=" << cfg.variant
                      << " input_mode=" << cfg.input_mode
                      << " op_threads=" << cfg.op_threads
                      << " dag_parallel=" << cfg.dag_parallel
                      << " repeat=" << cfg.repeat
                      << " warmup=" << cfg.warmup
                      << " hw_threads=" << hw_threads
                      << std::endl;
        }

        SecretKey sk;
        EvalKey ek;
        ek.emplacebkfft<lvl02param>(sk);
        ek.emplaceiksk<lvl20param>(sk);
        init_constants(sk);

        auto dataset = make_dataset(cfg, sk);

        for (int i = 0; i < cfg.warmup; ++i) {
            (void)run_wave(cfg, dataset, ek);
        }

        std::vector<double> repeat_ms;
        repeat_ms.reserve(cfg.repeat);
        std::uint64_t checksum = 0;

        const auto total_t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < cfg.repeat; ++i) {
            WaveResult wave = run_wave(cfg, dataset, ek);
            repeat_ms.push_back(wave.ms);
            checksum ^= mix(wave.checksum + static_cast<std::uint64_t>(i));
            if (!cfg.quiet) {
                std::cout << "[repeat " << (i + 1) << "/" << cfg.repeat
                          << "] " << std::fixed << std::setprecision(3)
                          << wave.ms << " ms" << std::endl;
            }
        }
        const auto total_t1 = std::chrono::steady_clock::now();

        const double total_ms =
            std::chrono::duration<double, std::milli>(total_t1 - total_t0).count();
        const double avg_repeat_ms = total_ms / static_cast<double>(cfg.repeat);
        const double avg_op_ms = avg_repeat_ms / static_cast<double>(cfg.dag_parallel);
        const double total_ops = static_cast<double>(cfg.repeat * cfg.dag_parallel);
        const double throughput = total_ops * 1000.0 / total_ms;

        std::string row = make_csv_row(
            cfg, total_ms, avg_repeat_ms, avg_op_ms, throughput,
            checksum, hw_threads, omp_max_threads);

        if (cfg.quiet) {
            std::cout << row << std::endl;
        } else {
            std::cout << "\n" << csv_header() << "\n" << row << std::endl;
        }

        if (!cfg.csv_path.empty()) {
            const bool exists = file_exists(cfg.csv_path);
            std::ofstream csv(cfg.csv_path, std::ios::app);
            if (!csv) throw std::runtime_error("cannot open csv: " + cfg.csv_path);
            if (!exists) csv << csv_header() << '\n';
            csv << row << '\n';
        }

        return 0;
    } catch (const std::exception &e) {
        std::cerr << "[error] " << e.what() << std::endl;
        return 1;
    }
}
