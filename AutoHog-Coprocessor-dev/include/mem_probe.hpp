// mem_probe.hpp
#pragma once
// Minimal, header-only memory snapshot helper for Linux.
// Provides: print_mem(const char*), print_mem(const std::string&),
//           struct MemScope { MemScope(const char* tag); ~MemScope(); }

#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#if !defined(_WIN32)
#include <unistd.h>
#endif

namespace _memprobe_detail
{

    struct MemInfoKB
    {
        long rss_kb = 0;      // Resident Set Size
        long vms_kb = 0;      // VmSize (virtual)
        long data_kb = 0;     // VmData (heap/anon)
        long peak_rss_kb = 0; // VmHWM
    };

    inline long _read_kb_from_status(const std::string &key)
    {
        std::ifstream ifs("/proc/self/status");
        if (!ifs.good())
            return 0;
        std::string line;
        while (std::getline(ifs, line))
        {
            if (line.rfind(key, 0) == 0)
            {
                // format: Key:\t  12345 kB
                std::istringstream iss(line.substr(key.size()));
                long val = 0;
                std::string unit;
                iss >> val >> unit;
                return val; // already KB
            }
        }
        return 0;
    }

    // VmSize is not in /proc/self/status? (it is, as "VmSize").
    // We also provide a fallback via /proc/self/statm (pages).
    inline long _read_vms_kb_statm()
    {
#if defined(_WIN32)
        return 0;
#else
        std::ifstream ifs("/proc/self/statm");
        if (!ifs.good())
            return 0;
        long size_pages = 0;
        ifs >> size_pages;
        long page_kb = long(sysconf(_SC_PAGESIZE) / 1024);
        return size_pages * page_kb;
#endif
    }

    inline MemInfoKB read_meminfo()
    {
        MemInfoKB mi;
        mi.rss_kb = _read_kb_from_status("VmRSS:");
        mi.vms_kb = _read_kb_from_status("VmSize:");
        if (mi.vms_kb <= 0)
            mi.vms_kb = _read_vms_kb_statm();
        mi.data_kb = _read_kb_from_status("VmData:");
        mi.peak_rss_kb = _read_kb_from_status("VmHWM:");
        return mi;
    }

    inline double kb_to_mb(long kb) { return double(kb) / 1024.0; }
    inline double kb_to_gb(long kb) { return double(kb) / (1024.0 * 1024.0); }

    inline void print_line(const char *tag, const MemInfoKB &mi)
    {
        // match your existing style as closely as possible
        std::fprintf(
            stderr,
            "%s RSS=%.3f MB  VMS=%.2f MB  DATA≈%.0f MB  PeakRSS=%.0f MB\n",
            tag,
            kb_to_mb(mi.rss_kb),
            kb_to_mb(mi.vms_kb),
            kb_to_mb(mi.data_kb),
            kb_to_mb(mi.peak_rss_kb));
        std::fflush(stderr);
    }

} // namespace _memprobe_detail

// ---- Public API ----

inline void print_mem(const char *tag)
{
    using namespace _memprobe_detail;
    print_line(tag, read_meminfo());
}

inline void print_mem(const std::string &tag)
{
    print_mem(tag.c_str());
}

// RAII scope logger: prints enter/exit with deltas
struct MemScope
{
    std::string tag_;
    _memprobe_detail::MemInfoKB start_;
    std::chrono::steady_clock::time_point t0_;

    explicit MemScope(const char *tag)
        : tag_(tag), start_(_memprobe_detail::read_meminfo()), t0_(std::chrono::steady_clock::now())
    {
        std::fprintf(stderr, "[MemScope enter] %s  ", tag_.c_str());
        _memprobe_detail::print_line("", start_); // line already contains numbers
    }

    explicit MemScope(const std::string &tag)
        : MemScope(tag.c_str()) {}

    ~MemScope()
    {
        using namespace _memprobe_detail;
        auto end = read_meminfo();
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0_).count();

        double d_rss_mb = kb_to_mb(end.rss_kb - start_.rss_kb);
        double d_vms_mb = kb_to_mb(end.vms_kb - start_.vms_kb);

        std::fprintf(stderr,
                     "[MemScope exit ] %s  RSS=%.3f MB  VMS=%.2f MB  DATA≈%.0f MB  PeakRSS=%.0f MB  ΔRSS=%.2f MB  ΔVMS=%.2f MB  (%.0f ms)\n",
                     tag_.c_str(),
                     kb_to_mb(end.rss_kb),
                     kb_to_mb(end.vms_kb),
                     kb_to_mb(end.data_kb),
                     kb_to_mb(end.peak_rss_kb),
                     d_rss_mb, d_vms_mb, ms);
        std::fflush(stderr);
    }
};

// Optional no-op switches if you need to temporarily silence logs
#ifdef MEM_PROBE_SILENT
#undef print_mem
inline void print_mem(const char *) {}
inline void print_mem(const std::string &) {}
struct MemScope
{
    explicit MemScope(const char *) {}
    explicit MemScope(const std::string &) {}
    ~MemScope() {}
};
#endif

// ---- add to mem_probe.hpp (anywhere after the functions are defined) ----
#ifndef MEM_SNAPSHOT
#define MEM_SNAPSHOT(tag) ::print_mem(tag)
#endif

#ifndef MEM_SCOPE

#define MEM_SCOPE(tag) \
    ::MemScope _memscope_##__LINE__ { tag }
#endif

#ifdef MEM_PROBE_SILENT
#undef MEM_SNAPSHOT
#undef MEM_SCOPE
#define MEM_SNAPSHOT(tag) \
    do                    \
    {                     \
    } while (0)
#define MEM_SCOPE(tag) \
    if (true)          \
        ;              \
    else               \
        struct _dummy  \
        {              \
        }
#endif
// ------------------------------------------------------------------------
