#pragma once

#include <chrono>
#include <string>

namespace log
{

enum class Level : char
{
    Debug,
    Info,
    Warn,
    Error
};

void Init(Level level, const std::string& file);
void SetLevel(Level level);
void Debug(const std::string& msg) noexcept;
void Info(const std::string& msg) noexcept;
void Warn(const std::string& msg) noexcept;
void Error(const std::string& msg) noexcept;

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

} // namespace log
