#pragma once

#include <algorithm>
#include <cstdint>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace TFHEpp {

struct HomTraceLayer {
    int level = 0;
    std::string phase;
    std::uint64_t scheduled_br = 0;
    std::uint64_t actual_br = 0;
    std::map<std::string, std::uint64_t> scheduled_primitives;
    std::map<std::string, std::uint64_t> actual_primitives;
};

class HomTrace {
public:
    explicit HomTrace(bool schedule_only = false) : schedule_only_(schedule_only) {}

    void reset(bool schedule_only)
    {
        std::lock_guard<std::mutex> lock(mu_);
        schedule_only_ = schedule_only;
        current_level_ = 0;
        output_depth_ = 0;
        prefix_.clear();
        current_phase_.clear();
        layers_.clear();
    }

    bool schedule_only() const
    {
        std::lock_guard<std::mutex> lock(mu_);
        return schedule_only_;
    }

    int next_layer(const std::string &phase)
    {
        std::lock_guard<std::mutex> lock(mu_);
        ++current_level_;
        current_phase_ = phase;
        if (schedule_only_) return current_level_;
        if (layers_.size() < static_cast<std::size_t>(current_level_))
            layers_.resize(static_cast<std::size_t>(current_level_));
        HomTraceLayer &layer = layers_[static_cast<std::size_t>(current_level_ - 1)];
        layer.level = current_level_;
        const std::string full_phase = prefix_ + phase;
        if (layer.phase.empty()) layer.phase = full_phase;
        else if (layer.phase.find(full_phase) == std::string::npos)
            layer.phase += "+" + full_phase;
        return current_level_;
    }

    void set_level(int level)
    {
        if (level < 0) throw std::runtime_error("HomTrace: logical level must be non-negative");
        std::lock_guard<std::mutex> lock(mu_);
        current_level_ = level;
    }

    void ensure_level(int level)
    {
        std::lock_guard<std::mutex> lock(mu_);
        current_level_ = std::max(current_level_, level);
    }

    int current_level() const
    {
        std::lock_guard<std::mutex> lock(mu_);
        return current_level_;
    }

    void record_scheduled(const std::string &primitive, std::uint64_t count)
    {
        std::lock_guard<std::mutex> lock(mu_);
        HomTraceLayer &layer = current_layer_locked();
        layer.scheduled_br += count;
        layer.scheduled_primitives[primitive] += count;
    }

    void record_scheduled_at(int level, const std::string &primitive,
                             std::uint64_t count = 1)
    {
        if (level <= 0) throw std::runtime_error("HomTrace: BR level must be positive");
        std::lock_guard<std::mutex> lock(mu_);
        if (layers_.size() < static_cast<std::size_t>(level))
            layers_.resize(static_cast<std::size_t>(level));
        HomTraceLayer &layer = layers_[static_cast<std::size_t>(level - 1)];
        layer.level = level;
        const std::string full_phase = prefix_ + current_phase_;
        if (!full_phase.empty()) {
            if (layer.phase.empty()) layer.phase = full_phase;
            else if (layer.phase.find(full_phase) == std::string::npos)
                layer.phase += "+" + full_phase;
        }
        layer.scheduled_br += count;
        layer.scheduled_primitives[primitive] += count;
    }

    void record_actual(const std::string &primitive, std::uint64_t count)
    {
        std::lock_guard<std::mutex> lock(mu_);
        HomTraceLayer &layer = current_layer_locked();
        layer.actual_br += count;
        layer.actual_primitives[primitive] += count;
    }

    std::string set_prefix(const std::string &prefix)
    {
        std::lock_guard<std::mutex> lock(mu_);
        std::string previous = prefix_;
        prefix_ = prefix;
        return previous;
    }

    int depth() const
    {
        std::lock_guard<std::mutex> lock(mu_);
        return output_depth_;
    }

    void mark_outputs_at_current_level()
    {
        std::lock_guard<std::mutex> lock(mu_);
        output_depth_ = std::max(output_depth_, current_level_);
    }

    void mark_output_level(int level)
    {
        std::lock_guard<std::mutex> lock(mu_);
        output_depth_ = std::max(output_depth_, level);
    }

    std::vector<HomTraceLayer> layers() const
    {
        std::lock_guard<std::mutex> lock(mu_);
        return layers_;
    }

    std::uint64_t scheduled_total() const
    {
        std::lock_guard<std::mutex> lock(mu_);
        std::uint64_t total = 0;
        for (const auto &layer : layers_) total += layer.scheduled_br;
        return total;
    }

    std::uint64_t actual_total() const
    {
        std::lock_guard<std::mutex> lock(mu_);
        std::uint64_t total = 0;
        for (const auto &layer : layers_) total += layer.actual_br;
        return total;
    }

private:
    HomTraceLayer &current_layer_locked()
    {
        if (layers_.empty())
            throw std::runtime_error("HomTrace: BR recorded without a logical layer");
        return layers_[static_cast<std::size_t>(current_level_ - 1)];
    }

    mutable std::mutex mu_;
    bool schedule_only_ = false;
    int current_level_ = 0;
    int output_depth_ = 0;
    std::string prefix_;
    std::string current_phase_;
    std::vector<HomTraceLayer> layers_;
};

inline HomTrace *g_active_hom_trace = nullptr;
inline std::mutex g_active_hom_trace_mu;

class ScopedHomTrace {
public:
    explicit ScopedHomTrace(HomTrace &trace)
    {
        std::lock_guard<std::mutex> lock(g_active_hom_trace_mu);
        if (g_active_hom_trace != nullptr)
            throw std::runtime_error("HomTrace: nested or concurrent trace is not supported");
        g_active_hom_trace = &trace;
    }

    ~ScopedHomTrace()
    {
        std::lock_guard<std::mutex> lock(g_active_hom_trace_mu);
        g_active_hom_trace = nullptr;
    }

    ScopedHomTrace(const ScopedHomTrace &) = delete;
    ScopedHomTrace &operator=(const ScopedHomTrace &) = delete;
};

class ScopedHomTracePrefix {
public:
    explicit ScopedHomTracePrefix(const std::string &prefix)
        : trace_(g_active_hom_trace)
    {
        if (trace_ != nullptr) previous_ = trace_->set_prefix(prefix);
    }

    ~ScopedHomTracePrefix()
    {
        if (trace_ != nullptr) trace_->set_prefix(previous_);
    }

private:
    HomTrace *trace_ = nullptr;
    std::string previous_;
};

inline bool hom_trace_active() { return g_active_hom_trace != nullptr; }
inline bool hom_trace_schedule_only()
{
    return g_active_hom_trace != nullptr && g_active_hom_trace->schedule_only();
}
inline void hom_trace_next_layer(const std::string &phase)
{
    if (g_active_hom_trace != nullptr) g_active_hom_trace->next_layer(phase);
}
inline int hom_trace_current_level()
{
    return g_active_hom_trace != nullptr ? g_active_hom_trace->current_level() : 0;
}
inline void hom_trace_set_level(int level)
{
    if (g_active_hom_trace != nullptr) g_active_hom_trace->set_level(level);
}
inline void hom_trace_ensure_level(int level)
{
    if (g_active_hom_trace != nullptr) g_active_hom_trace->ensure_level(level);
}
inline void hom_trace_scheduled(const std::string &primitive, std::uint64_t count = 1)
{
    if (g_active_hom_trace != nullptr)
        g_active_hom_trace->record_scheduled(primitive, count);
}
inline void hom_trace_scheduled_at(int level, const std::string &primitive,
                                   std::uint64_t count = 1)
{
    if (g_active_hom_trace != nullptr)
        g_active_hom_trace->record_scheduled_at(level, primitive, count);
}
inline void hom_trace_actual(const std::string &primitive, std::uint64_t count = 1)
{
    if (g_active_hom_trace != nullptr)
        g_active_hom_trace->record_actual(primitive, count);
}
inline void hom_trace_mark_outputs()
{
    if (g_active_hom_trace != nullptr)
        g_active_hom_trace->mark_outputs_at_current_level();
}
inline void hom_trace_mark_output_level(int level)
{
    if (g_active_hom_trace != nullptr)
        g_active_hom_trace->mark_output_level(level);
}

inline constexpr std::uint64_t hom_trace_level_magic = UINT64_C(0x6a17000000000000);
inline constexpr std::uint64_t hom_trace_level_mask = UINT64_C(0xffff000000000000);

inline std::uint64_t hom_trace_encode_level(int level)
{
    return level <= 0 ? 0 : hom_trace_level_magic | static_cast<std::uint64_t>(level);
}

inline int hom_trace_decode_level(std::uint64_t value)
{
    auto decode = [](std::uint64_t candidate) -> int {
        if ((candidate & hom_trace_level_mask) != hom_trace_level_magic) return 0;
        return static_cast<int>(candidate & ~hom_trace_level_mask);
    };
    const int direct = decode(value);
    return direct != 0 ? direct : decode(UINT64_C(0) - value);
}

} // namespace TFHEpp
