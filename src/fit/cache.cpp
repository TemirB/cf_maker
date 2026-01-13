#include "fit/cache.h"
#include "fit/json.h"
#include "helpers.h"

#include <fstream>
#include <TSystem.h>
#include <iostream>

using json = nlohmann::json;

static bool FileExists(const std::string& f) {
    return !gSystem->AccessPathName(f.c_str());
}

int FindRapidityIndex(double y0, double y1)
{
    // центр диапазона
    double y = 0.5 * (y0 + y1);

    // защитимся от мусора double
    if (y < -1.000001 || y > 1.000001)
        return -1;

    // переводим [-1, +1] -> [0, 10)
    int idx = int( std::floor( (y + 1.0) / 0.2 ) );

    // защита от y=+1.0
    if (idx == rapiditySize)
        idx = rapiditySize - 1;

    return (idx >= 0 && idx < rapiditySize) ? idx : -1;
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

    int ok = 0;
    int dropCharge = 0, dropCent = 0, dropKt = 0, dropY = 0;

    for (auto& p : J["data"]) {

        // charge
        std::string c = p["charge"].get<std::string>();
        int ch;
        if (c == "pos" || c == "+") ch = 0;
        else if (c == "neg" || c == "-") ch = 1;
        else { dropCharge++; continue; }

        // centrality
        int cent = FindCentralityIndex(p["centrality"]["name"].get<std::string>());
        if (cent < 0 || cent >= centralitySize) {
            dropCent++;
            continue;
        }

        // kt
        double kt0 = p["kt"]["range"][0];
        double kt1 = p["kt"]["range"][1];
        int kt = FindKtIndex(kt0, kt1);
        if (kt < 0 || kt >= ktSize) {
            dropKt++;
            continue;
        }

        // rapidity
        double y0 = p["rapidity"]["range"][0];
        double y1 = p["rapidity"]["range"][1];
        int y = FindRapidityIndex(y0, y1);
        if (y < 0 || y >= rapiditySize) {
            dropY++;
            continue;
        }

        // store
        fitRes[ch][cent][kt][y] = FitResultFromJson(p["fit"]);
        ok++;
    }

    std::cout
        << "Cache load summary:\n"
        << "  ok         = " << ok << "\n"
        << "  dropCharge = " << dropCharge << "\n"
        << "  dropCent   = " << dropCent << "\n"
        << "  dropKt     = " << dropKt << "\n"
        << "  dropY      = " << dropY << std::endl;

    if (ok == 0) {
        reason = "cache JSON loaded but no bins matched";
        return CacheStatus::CorruptJson;
    }

    reason = "ok";
    return CacheStatus::Ok;
}
