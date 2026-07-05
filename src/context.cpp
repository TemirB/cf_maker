#include "context.h"

#include <limits.h>
#include <libgen.h>
#include <filesystem>
#include <string>

#include "helpers.h"

#if defined(__APPLE__)
#include <mach-o/dyld.h>
#elif defined(__linux__)
#include <unistd.h>
#endif

std::filesystem::path GetExecutablePath() {
#if defined(__APPLE__)
    uint32_t size = 0;
    _NSGetExecutablePath(nullptr, &size);

    std::string buffer(size, '\0');

    if (_NSGetExecutablePath(buffer.data(), &size) != 0) {
        throw std::runtime_error("Failed to get executable path on macOS");
    }

    // buffer may contain trailing '\0'
    buffer.resize(std::char_traits<char>::length(buffer.c_str()));

    return std::filesystem::canonical(buffer);

#elif defined(__linux__)
    char buf[4096];

    ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (len == -1) {
        throw std::runtime_error("Failed to read /proc/self/exe");
    }

    buf[len] = '\0';
    return std::filesystem::canonical(buf);

#else
    throw std::runtime_error("Unsupported platform");
#endif
}

std::string GetExeDir() {
    std::string path =  GetExecutablePath();

    auto pos = path.find_last_of('/');
    if (pos == std::string::npos) {
        return ".";
    }

    return path.substr(0, pos);
}

Context BuildContext(char** argv) {
    Context ctx;

    std::string resDir  = GetExeDir() + "/../results";

    EnsureDir(resDir);

    ctx.inputFile = argv[1];
    ctx.outDir    = resDir + "/" + argv[2];
    EnsureDir(ctx.outDir);


    ctx.bining.type = argv[3];
    if (std::strcmp(ctx.bining.type, "kt") == 0) {
        ctx.bining.count = Kt::kCount;
        ctx.bining.names = Kt::kNames;
        ctx.bining.values = Kt::kValues;
        ctx.bining.fileNames = Kt::kFileNames;
    } else if (std::strcmp(ctx.bining.type, "rapidity") == 0) {
        ctx.bining.count = Rapidity::kCount;
        ctx.bining.names = Rapidity::kNames;
        ctx.bining.values = Rapidity::kValues;
        ctx.bining.fileNames = Rapidity::kFileNames;
    }


    ctx.fitRes = FitGrid(
        Charge::kCount,
        std::vector<std::vector<FitResult>>(
            Centrality::kCount,
            std::vector<FitResult>(ctx.bining.count)
        )
    );
    return ctx;
}