#pragma once

#include <chrono>
#include <string>

namespace logging
{

enum class Level : char
{
    Debug,
    Info,
    Warn,
    Error
};

void init(Level level, const std::string& file);
void set_level(Level level);
void debug(const std::string& msg) noexcept;
void info(const std::string& msg) noexcept;
void warn(const std::string& msg) noexcept;
void error(const std::string& msg) noexcept;

[[nodiscard]] Level parse_level(const std::string& name);

class ScopedTimer
{
  public:
    explicit ScopedTimer(std::string name);
    ~ScopedTimer();

    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;
    ScopedTimer(ScopedTimer&&) = delete;
    ScopedTimer& operator=(ScopedTimer&&) = delete;

  private:
    std::string name_;
    std::chrono::steady_clock::time_point start_;
};

} // namespace logging
