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
log::Level gLevel = log::Level::Info;
std::unique_ptr<std::ofstream> gFile;

const char* LevelName(log::Level level)
{
    switch (level) {
    case log::Level::Debug:
        return "DEBUG";
    case log::Level::Info:
        return "INFO";
    case log::Level::Warn:
        return "WARN";
    case log::Level::Error:
        return "ERROR";
    }
    return "?";
}

std::string Timestamp()
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

void Write(log::Level level, const std::string& msg) noexcept
{
    try {
        std::lock_guard<std::mutex> lock(gMutex);
        if (level < gLevel) {
            return;
        }

        std::ostringstream line;
        line << Timestamp() << " [" << LevelName(level) << "] [t=" << std::this_thread::get_id()
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

namespace log
{

void Init(Level level, const std::string& file)
{
    std::lock_guard<std::mutex> lock(gMutex);
    gLevel = level;
    gFile.reset();
    if (!file.empty()) {
        gFile = std::make_unique<std::ofstream>(file, std::ios::out | std::ios::app);
    }
}

void SetLevel(Level level)
{
    std::lock_guard<std::mutex> lock(gMutex);
    gLevel = level;
}

void Debug(const std::string& msg) noexcept
{
    Write(Level::Debug, msg);
}

void Info(const std::string& msg) noexcept
{
    Write(Level::Info, msg);
}

void Warn(const std::string& msg) noexcept
{
    Write(Level::Warn, msg);
}

void Error(const std::string& msg) noexcept
{
    Write(Level::Error, msg);
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
        Info(name_ + ": done in " + std::to_string(seconds) + " s");
    } catch (const std::exception& e) {
        std::cerr << "ScopedTimer: " << e.what() << "\n";
    } catch (...) {
        std::cerr << "ScopedTimer: unknown error\n";
    }
}

} // namespace log
