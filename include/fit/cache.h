#include <string>

#include "fit/types.h"
#include "helpers.h"

enum class CacheStatus {
    Ok,
    NoFiles,
    HashMismatch,
    CorruptJson
};

CacheStatus LoadFitCache(
    const std::string& inputFile, 
    const std::string& outDir,
    FitResult fitRes[chargeSize][centralitySize][ktSize][rapiditySize],
    std::string& reason
);