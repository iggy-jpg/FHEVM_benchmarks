#include "Autohog_Cop.hpp"
#include "HomOP.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <unistd.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace TFHEpp;
using boost::multiprecision::cpp_int;

namespace {

struct Options {
    int width = 256;
    int repeat = 100;
    int threads = 48;
    int progress_every = 1;
    std::uint64_t seed = 0x20260726ULL;
    long stop_cgroup_kib = 15L * 1024L * 1024L;
    std::string csv_path;
};

struct MemorySnapshot {
    long vm_rss_kib = 0;
    long vm_hwm_kib = 0;
    long ru_maxrss_kib = 0;
    long cgroup_current_kib = 0;
    long cgroup_peak_kib = 0;
    int process_threads = 1;
};

struct ProcessStatus {
    long vm_rss_kib = 0;
    long vm_hwm_kib = 0;
    int threads = 1;
};

Options parse_options(int argc, char **argv)
{
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto value = [&]() -> std::string {
            if (++i >= argc) throw std::runtime_error(arg + " requires a value");
            return argv[i];
        };
        if (arg == "--width") {
            options.width = std::stoi(value());
        } else if (arg == "--repeat") {
            options.repeat = std::stoi(value());
        } else if (arg == "--threads") {
            options.threads = std::stoi(value());
        } else if (arg == "--progress-every") {
            options.progress_every = std::stoi(value());
        } else if (arg == "--seed") {
            options.seed = std::stoull(value(), nullptr, 0);
        } else if (arg == "--stop-cgroup-kib") {
            options.stop_cgroup_kib = std::stol(value());
        } else if (arg == "--csv") {
            options.csv_path = value();
        } else if (arg == "--help" || arg == "-h") {
            std::cout
                << "usage: test_DIV_MEMORY_PROFILE [--width 256] [--repeat 100] "
                   "[--threads 48] [--progress-every 1] [--seed N] "
                   "[--stop-cgroup-kib N] --csv PATH\n";
            std::exit(0);
        } else {
            throw std::runtime_error("unknown argument: " + arg);
        }
    }
    if (options.width < 4 || options.width > 256)
        throw std::runtime_error("width must be in [4, 256]");
    if (options.repeat <= 0 || options.threads <= 0 || options.progress_every < 0)
        throw std::runtime_error("repeat/threads must be positive and progress non-negative");
    if (options.stop_cgroup_kib < 0)
        throw std::runtime_error("stop-cgroup-kib must be non-negative");
    if (options.csv_path.empty()) throw std::runtime_error("--csv is required");
    return options;
}

long read_long_file(const std::string &path)
{
    std::ifstream input(path);
    long value = 0;
    return input >> value ? value : 0;
}

std::string cgroup_v2_path()
{
#ifdef __linux__
    std::ifstream input("/proc/self/cgroup");
    std::string line;
    while (std::getline(input, line)) {
        const auto marker = line.find("0::");
        if (marker == 0) return "/sys/fs/cgroup" + line.substr(3);
    }
#endif
    return {};
}

ProcessStatus read_process_status()
{
    ProcessStatus result;
#ifdef __linux__
    std::ifstream input("/proc/self/status");
    std::string key;
    while (input >> key) {
        if (key == "VmRSS:") {
            input >> result.vm_rss_kib;
        } else if (key == "VmHWM:") {
            input >> result.vm_hwm_kib;
        } else if (key == "Threads:") {
            input >> result.threads;
        }
        input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
#endif
    return result;
}

MemorySnapshot snapshot(const std::string &cgroup_path)
{
    const ProcessStatus status = read_process_status();
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    MemorySnapshot result;
    result.vm_rss_kib = status.vm_rss_kib;
    result.vm_hwm_kib = status.vm_hwm_kib;
#ifdef __APPLE__
    result.ru_maxrss_kib = usage.ru_maxrss / 1024;
#else
    result.ru_maxrss_kib = usage.ru_maxrss;
#endif
    result.process_threads = status.threads;
    if (!cgroup_path.empty()) {
        result.cgroup_current_kib = read_long_file(cgroup_path + "/memory.current") / 1024;
        result.cgroup_peak_kib = read_long_file(cgroup_path + "/memory.peak") / 1024;
    }
    return result;
}

std::string wall_clock_time()
{
    const std::time_t now = std::time(nullptr);
    std::tm local{};
#ifdef _WIN32
    localtime_s(&local, &now);
#else
    localtime_r(&now, &local);
#endif
    std::ostringstream out;
    out << std::put_time(&local, "%Y-%m-%dT%H:%M:%S");
    return out.str();
}

cpp_int mask_bits(int width) { return (cpp_int(1) << width) - 1; }

cpp_int random_value(int width, std::mt19937_64 &rng)
{
    cpp_int value = 0;
    for (int generated = 0; generated < width; generated += 64) {
        value <<= 64;
        value |= rng();
    }
    return value & mask_bits(width);
}

std::string to_hex(const cpp_int &value)
{
    std::ostringstream out;
    out << std::hex << value;
    return out.str();
}

std::vector<TLWE<lvl2param>> encrypt_value(
    const cpp_int &value, int width, const SecretKey &sk)
{
    std::vector<TLWE<lvl2param>> encrypted(width);
    for (int i = 0; i < width; ++i) {
        const auto bit = static_cast<std::uint8_t>(
            ((value >> (width - 1 - i)) & 1) != 0);
        encrypted[i] = bootsSymEncrypt<lvl2param>({bit}, sk)[0];
    }
    return encrypted;
}

cpp_int decrypt_value(
    const std::vector<TLWE<lvl2param>> &encrypted, const SecretKey &sk)
{
    cpp_int value = 0;
    for (const auto &bit : encrypted)
        value = (value << 1) | (bootsSymDecrypt<lvl2param>({bit}, sk)[0] & 1u);
    return value;
}

void write_event(std::ofstream &csv, int width, int sample, const std::string &stage,
                 int threads_hint, double div_ms, int correct,
                 const MemorySnapshot &memory, const cpp_int *a, const cpp_int *b,
                 const std::string &status)
{
    csv << wall_clock_time() << ',' << width << ',' << sample << ',' << stage << ','
        << threads_hint << ',' << memory.process_threads << ',' << std::fixed
        << std::setprecision(6) << div_ms << ',' << correct << ','
        << memory.vm_rss_kib << ',' << memory.vm_hwm_kib << ','
        << memory.ru_maxrss_kib << ',' << memory.cgroup_current_kib << ','
        << memory.cgroup_peak_kib << ',' << (a == nullptr ? "" : to_hex(*a)) << ','
        << (b == nullptr ? "" : to_hex(*b)) << ',' << status << '\n';
    csv.flush();
}

bool over_memory_stop(const MemorySnapshot &memory, const Options &options)
{
    return options.stop_cgroup_kib > 0 && memory.cgroup_current_kib > 0 &&
           memory.cgroup_current_kib >= options.stop_cgroup_kib;
}

} // namespace

int main(int argc, char **argv)
{
    try {
        const Options options = parse_options(argc, argv);
#ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(options.threads);
#endif
        const std::string cgroup_path = cgroup_v2_path();
        std::ofstream csv(options.csv_path);
        if (!csv) throw std::runtime_error("cannot open CSV: " + options.csv_path);
        csv << "timestamp,width,sample,stage,threads_hint,process_threads,div_ms,"
               "correct,vm_rss_kib,vm_hwm_kib,ru_maxrss_kib,cgroup_current_kib,"
               "cgroup_peak_kib,a_hex,b_hex,status\n";
        csv.flush();

        std::cout << "[CONFIG] width=" << options.width << " repeat=" << options.repeat
                  << " threads=" << options.threads << " seed=" << options.seed
                  << " stop_cgroup_kib=" << options.stop_cgroup_kib
                  << " cgroup=" << (cgroup_path.empty() ? "unavailable" : cgroup_path)
                  << std::endl;

        SecretKey sk;
        EvalKey ek;
        ek.emplacebkfft<lvl02param>(sk);
        ek.emplaceiksk<lvl20param>(sk);
        init_constants(sk);
        write_event(csv, options.width, 0, "after_key", options.threads, 0.0, -1,
                    snapshot(cgroup_path), nullptr, nullptr, "ok");

        std::mt19937_64 rng(options.seed);
        int correct_count = 0;
        for (int sample = 1; sample <= options.repeat; ++sample) {
            const cpp_int a = random_value(options.width, rng);
            const cpp_int b = random_value(options.width, rng);
            double div_ms = 0.0;
            bool correct = false;

            write_event(csv, options.width, sample, "before_encrypt", options.threads,
                        div_ms, -1, snapshot(cgroup_path), &a, &b, "running");
            {
                const auto encrypted_a = encrypt_value(a, options.width, sk);
                const auto encrypted_b = encrypt_value(b, options.width, sk);
                write_event(csv, options.width, sample, "before_div", options.threads,
                            div_ms, -1, snapshot(cgroup_path), &a, &b, "running");

                const auto start = std::chrono::steady_clock::now();
                HomDivResult result = HomDiv_pal(
                    encrypted_a, encrypted_b, options.width, ek, options.threads);
                const auto end = std::chrono::steady_clock::now();
                div_ms = std::chrono::duration<double, std::milli>(end - start).count();
                write_event(csv, options.width, sample, "after_div", options.threads,
                            div_ms, -1, snapshot(cgroup_path), &a, &b, "running");

                const cpp_int expected_q = b == 0 ? cpp_int(0) : a / b;
                const cpp_int expected_r = b == 0 ? a : a % b;
                correct = result.quotient.size() == static_cast<std::size_t>(options.width) &&
                          result.remainder.size() == static_cast<std::size_t>(options.width) &&
                          decrypt_value(result.quotient, sk) == expected_q &&
                          decrypt_value(result.remainder, sk) == expected_r;
                write_event(csv, options.width, sample, "after_decrypt", options.threads,
                            div_ms, correct ? 1 : 0, snapshot(cgroup_path), &a, &b,
                            correct ? "ok" : "incorrect");
            }

            const MemorySnapshot after_destroy = snapshot(cgroup_path);
            write_event(csv, options.width, sample, "after_destroy", options.threads,
                        div_ms, correct ? 1 : 0, after_destroy, &a, &b,
                        correct ? "ok" : "incorrect");
            if (!correct)
                throw std::runtime_error("DIV correctness failure at sample " +
                                         std::to_string(sample));
            ++correct_count;
            if (after_destroy.process_threads > options.threads)
                throw std::runtime_error(
                    "DIV exceeded thread budget: observed " +
                    std::to_string(after_destroy.process_threads) + ", budget " +
                    std::to_string(options.threads));

            if (options.progress_every > 0 &&
                (sample == options.repeat || sample % options.progress_every == 0)) {
                std::cout << "[PROGRESS] sample=" << sample << '/' << options.repeat
                          << " correct=" << correct_count << '/' << sample
                          << " div_ms=" << std::fixed << std::setprecision(3) << div_ms
                          << " rss_kib=" << after_destroy.vm_rss_kib
                          << " hwm_kib=" << after_destroy.vm_hwm_kib
                          << " cgroup_kib=" << after_destroy.cgroup_current_kib
                          << " cgroup_peak_kib=" << after_destroy.cgroup_peak_kib
                          << std::endl;
            }
            if (over_memory_stop(after_destroy, options)) {
                std::cerr << "[MEMORY_LIMIT] cgroup memory reached soft stop at sample "
                          << sample << std::endl;
                return 75;
            }
        }
        std::cout << "[PASS] DIV memory profile completed correct=" << correct_count << '/'
                  << options.repeat << std::endl;
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "[FAIL] " << error.what() << std::endl;
        return 1;
    }
}
