#pragma once

#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

template <typename F> void ParallelFor(std::size_t count, std::size_t nThreads, F&& f)
{
    if (count == 0) {
        return;
    }

    if (nThreads == 0) {
        nThreads = std::thread::hardware_concurrency();
    }
    if (nThreads > count) {
        nThreads = count;
    }

    if (nThreads <= 1) {
        for (std::size_t i = 0; i < count; ++i) {
            f(i);
        }
        return;
    }

    std::atomic<std::size_t> index{0};
    std::mutex errorMutex;
    std::exception_ptr error;

    auto worker = [&]() {
        try {
            for (;;) {
                const std::size_t i = index.fetch_add(1, std::memory_order_relaxed);
                if (i >= count) {
                    return;
                }
                f(i);
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(errorMutex);
            if (!error) {
                error = std::current_exception();
            }
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(nThreads - 1);
    for (std::size_t t = 1; t < nThreads; ++t) {
        workers.emplace_back(worker);
    }
    worker();
    for (auto& w : workers) {
        w.join();
    }

    if (error) {
        std::rethrow_exception(error);
    }
}
