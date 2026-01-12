#include "fit/cache.h"
#include "fit/json.h"
#include "helpers.h"

#include <fstream>
#include <TSystem.h>

using json = nlohmann::json;

static bool FileExists(const std::string& f) {
    return !gSystem->AccessPathName(f.c_str());
}

CacheStatus LoadFitCache(
    const std::string& inputFile,
    const std::string& outDir,
    FitResult fitRes[chargeSize][centralitySize][ktSize][rapiditySize],
    std::string& reason
) {
    std::string cfFile   = outDir + "/cf3d.root";
    std::string jsonFile = outDir + "/fit3d.json";

    if (!FileExists(cfFile) || !FileExists(jsonFile)) {
        reason = "cache files missing";
        return CacheStatus::NoFiles;
    }

    json J;
    try {
        std::ifstream f(jsonFile);
        f >> J;
    } catch (...) {
        reason = "cannot parse fit3d.json";
        return CacheStatus::CorruptJson;
    }

    std::string want = Sha256OfFile(inputFile);
    std::string got  = J["meta"]["input"]["sha256"].get<std::string>();

    if (want != got) {
        reason = "input file hash mismatch";
        return CacheStatus::HashMismatch;
    }

    for (auto& p : J["data"]) {
        int ch = (p["charge"].get<std::string>() == "pos") ? 0 : 1;

        int cent = FindCentralityIndex(p["centrality"]["name"].get<std::string>());
        if (cent < 0 || cent >= centralitySize) continue;

        double kt0 = p["kt"]["range"][0];
        double kt1 = p["kt"]["range"][1];
        int kt = FindKtIndex(kt0, kt1);
        if (kt < 0 || kt >= ktSize) continue;

        double yLow = p["rapidity"]["range"][0];
        int y = int(std::floor((yLow + 1.0) / 0.2 + 1e-6));
        if (y < 0 || y >= rapiditySize) continue;

        fitRes[ch][cent][kt][y] = FitResultFromJson(p["fit"]);
    }

    reason = "ok";
    return CacheStatus::Ok;
}
