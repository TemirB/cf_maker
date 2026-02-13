#include <string>

#include "fit/types.h"

enum class CacheStatus {
    Ok,
    NoFiles,
    HashMismatch,
    CorruptJson
};

CacheStatus LoadFitCache(
    const std::string& inputFile, 
    const std::string& outDir,
    FitGrid fitRes,
    std::string& reason
);