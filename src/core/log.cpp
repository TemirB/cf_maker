#include "core/log.h"

#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>

namespace
{

std::mutex gMutex;
logging::Level gLevel = logging::Level::Info;
std::unique_ptr<std::ofstream> gFile;

const char* level_name(logging::Level level)
{
    switch (level) {
    case logging::Level::Debug:
        return "DEBUG";
    case logging::Level::Info:
        return "INFO";
    case logging::Level::Warn:
        return "WARN";
    case logging::Level::Error:
        return "ERROR";
    }
    return "?";
}

std::string timestamp()
{
    using namespace std::chrono;

    const auto now = system_clock::now();
    const auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;
    const std::time_t t = system_clock::to_time_t(now);

    std::tm tmBuf{};
    localtime_r(&t, &tmBuf);

    char buf[48];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d %02d:%02d:%02d.%03d", tmBuf.tm_year + 1900,
                  tmBuf.tm_mon + 1, tmBuf.tm_mday, tmBuf.tm_hour, tmBuf.tm_min, tmBuf.tm_sec,
                  static_cast<int>(ms.count()));
    return buf;
}

void write(logging::Level level, const std::string& msg) noexcept
{
    try {
        std::lock_guard<std::mutex> lock(gMutex);
        if (level < gLevel) {
            return;
        }

        std::ostringstream line;
        line << timestamp() << " [" << level_name(level) << "] [t=" << std::this_thread::get_id()
             << "] " << msg << "\n";
        std::cerr << line.str();
        if (gFile && gFile->is_open()) {
            *gFile << line.str();
            gFile->flush();
        }
    } catch (...) {
        std::cerr << "log write failed\n";
    }
}

} // namespace

namespace logging
{

void init(Level level, const std::string& file)
{
    std::lock_guard<std::mutex> lock(gMutex);
    gLevel = level;
    gFile.reset();
    if (!file.empty()) {
        gFile = std::make_unique<std::ofstream>(file, std::ios::out | std::ios::app);
    }
}

void set_level(Level level)
{
    std::lock_guard<std::mutex> lock(gMutex);
    gLevel = level;
}

void debug(const std::string& msg) noexcept
{
    write(Level::Debug, msg);
}

void info(const std::string& msg) noexcept
{
    write(Level::Info, msg);
}

void warn(const std::string& msg) noexcept
{
    write(Level::Warn, msg);
}

void error(const std::string& msg) noexcept
{
    write(Level::Error, msg);
}

Level parse_level(const std::string& name)
{
    if (name == "debug") {
        return Level::Debug;
    }
    if (name == "info") {
        return Level::Info;
    }
    if (name == "warning" || name == "warn") {
        return Level::Warn;
    }
    if (name == "error") {
        return Level::Error;
    }
    throw std::runtime_error("unknown log level: " + name);
}

ScopedTimer::ScopedTimer(std::string name)
    : name_(std::move(name)), start_(std::chrono::steady_clock::now())
{
}

ScopedTimer::~ScopedTimer()
{
    try {
        const double seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
        info(name_ + ": done in " + std::to_string(seconds) + " s");
    } catch (const std::exception& e) {
        std::cerr << "ScopedTimer: " << e.what() << "\n";
    } catch (...) {
        std::cerr << "ScopedTimer: unknown error\n";
    }
}

} // namespace logging
