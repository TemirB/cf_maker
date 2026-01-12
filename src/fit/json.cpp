#include "fit/json.h"

#include <fstream>
#include <iomanip>
#include <ctime>
#include <cmath>

#include "helpers.h"   // NowUTC(), Sha256OfFile(), chargeNames, centralityNames, ktNames …

using json = nlohmann::json;

// ----------------------
// helpers
// ----------------------

static std::string NowUTC() {
    std::time_t t = std::time(nullptr);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%FT%TZ", std::gmtime(&t));
    return buf;
}

// ----------------------
// FitResult → JSON
// ----------------------

json FitResultToJson(const FitResult& r) {
    return {
        {"ok", r.ok},
        {"valid", r.IsValid()},
        {"R",  {r.R[0],  r.R[1],  r.R[2]}},
        {"eR", {r.eR[0], r.eR[1], r.eR[2]}},
        {"lambda",  r.lambda},
        {"elambda", r.elambda},
        {"chi2",    r.chi2},
        {"ndf",     r.ndf},
        {"chi2_ndf", r.Chi2Ndf()},
        {"pvalue",  r.pvalue}
    };
}

// ----------------------
// Metadata
// ----------------------

json BuildMeta(const std::string& inputFile) {
    json m;

    m["schema_version"] = "1.0";
    m["analysis"] = "StHbtMaker → BPLCMS 3D CF → Gaussian HBT fit";

    m["generator"] = {
        {"program", "cf3d"},
        {"version", "v1.0"}
    };

    m["date_utc"] = NowUTC();

    m["input"] = {
        {"file", inputFile},
        {"sha256", Sha256OfFile(inputFile)}
    };

    // ---- axes ----
    m["axes"]["lcms"] = {"out","side","long"};
    m["axes"]["charge"] = {"pos","neg"};

    m["axes"]["centrality"] = {
        {{"name","0-10%"},{"impact_b_fm",{0.0,4.6}}},
        {{"name","10-30%"},{"impact_b_fm",{4.6,8.0}}},
        {{"name","30-50%"},{"impact_b_fm",{8.0,10.4}}},
        {{"name","50-80%"},{"impact_b_fm",{10.4,13.1}}}
    };

    m["axes"]["kt"] = {
        {0.15,0.25},
        {0.25,0.35},
        {0.35,0.45},
        {0.45,0.60}
    };

    m["axes"]["rapidity"] = {
        {"range", {-1.0,1.0}},
        {"bins", 10},
        {"step", 0.2}
    };

    m["cf"] = {
        {"type","BPLCMS3D"},
        {"q_axes",{"q_out","q_side","q_long"}},
        {"q_range",{-0.4,0.4}},
        {"q_bins",80}
    };

    m["fit_model"] = {
        {"name","Gaussian3D"},
        {"parameters",{"R_out","R_side","R_long","lambda"}}
    };

    return m;
}

// ----------------------
// One data point
// ----------------------

json BuildDataPoint(int ch, int cent, int kt, int y, const FitResult& r) {
    json d;

    d["charge"] = (ch==0) ? "pos" : "neg";

    static const std::array<std::pair<std::string,std::array<double,2>>,4> centInfo = {{
        {"0-10%",{0.0,4.6}},
        {"10-30%",{4.6,8.0}},
        {"30-50%",{8.0,10.4}},
        {"50-80%",{10.4,13.1}}
    }};

    d["centrality"] = {
        {"name",centInfo[cent].first},
        {"impact_b_fm",centInfo[cent].second}
    };

    static const double ktLow[]  = {0.15,0.25,0.35,0.45};
    static const double ktHigh[] = {0.25,0.35,0.45,0.60};

    d["kt"] = {
        {"range",{ktLow[kt],ktHigh[kt]}},
        {"mean",0.5*(ktLow[kt]+ktHigh[kt])}
    };

    double yLow  = -1.0 + 0.2*y;
    double yHigh = yLow + 0.2;

    d["rapidity"] = {
        {"range",{yLow,yHigh}},
        {"mean",0.5*(yLow+yHigh)}
    };

    d["fit"] = FitResultToJson(r);

    return d;
}

// ----------------------
// Write full file
// ----------------------

void WriteFitJson(const std::string& inputFile,
                  const std::string& outDir,
                  FitResult fitRes[chargeSize][centralitySize][ktSize][rapiditySize])
{
    json J;
    J["meta"] = BuildMeta(inputFile);
    J["data"] = json::array();

    for (int ch=0; ch<chargeSize; ch++)
    for (int c =0; c <centralitySize; c++)
    for (int k =0; k <ktSize; k++)
    for (int y =0; y <rapiditySize; y++) {
        J["data"].push_back(BuildDataPoint(ch,c,k,y, fitRes[ch][c][k][y]));
    }

    std::ofstream f(outDir + "/fit3d.json");
    f << std::setw(2) << J << "\n";
}

FitResult FitResultFromJson(const nlohmann::json& j) {
    FitResult r;

    r.ok      = j.at("ok").get<bool>();
    r.lambda  = j.at("lambda").get<double>();
    r.elambda = j.at("elambda").get<double>();
    r.chi2    = j.at("chi2").get<double>();
    r.ndf     = j.at("ndf").get<int>();
    r.pvalue  = j.value("pvalue", 0.0);

    for (int i = 0; i < 3; i++) {
        r.R[i]  = j.at("R")[i].get<double>();
        r.eR[i] = j.at("eR")[i].get<double>();
    }

    return r;
}