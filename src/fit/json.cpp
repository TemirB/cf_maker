#include "fit/json.h"

#include <fstream>
#include <iomanip>
#include <ctime>
#include <cmath>

#include "helpers/utils.h"

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
        {"R",  {r.R[0],  r.R[1],  r.R[2],  r.R[3],  r.R[4],  r.R[5]}},
        {"eR", {r.eR[0], r.eR[1], r.eR[2], r.eR[3], r.eR[4], r.eR[5]}},
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

    // m["axes"]["rapidity"] = {
    //     {"range", {-1.0,1.0}},
    //     {"bins", 10},
    //     {"step", 0.2}
    // };

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

json BuildDataPoint(int ch, int cent, int kt, const FitResult& r) {
    json d;

    d["charge"] = (ch==0) ? "pos" : "neg";

    static const std::array<std::pair<std::string,std::array<double,2>>,4> centInfo = {{
        {"0-10",{0.0,4.6}},
        {"10-30",{4.6,8.0}},
        {"30-50",{8.0,10.4}},
        {"50-80",{10.4,13.1}}
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

    // double yLow  = -1.0 + 0.2*y;
    // double yHigh = yLow + 0.2;

    // d["rapidity"] = {
    //     {"range",{yLow,yHigh}},
    //     {"mean",0.5*(yLow+yHigh)}
    // };

    d["fit"] = FitResultToJson(r);

    return d;
}

// ----------------------
// Write full file
// ----------------------

void WriteFitJson(const std::string& inputFile,
                  const std::string& outDir,
                  FitResult fitRes[chargeSize][centralitySize][ktSize])
{
    json J;
    J["meta"] = BuildMeta(inputFile);
    J["data"] = json::array();

    for (int chIdx = 0; chIdx < chargeSize; chIdx++)
    for (int centrIdx = 0; centrIdx < centralitySize; centrIdx++)
    for (int ktIdx = 0; ktIdx < ktSize; ktIdx++) {
        J["data"].push_back(BuildDataPoint(chIdx,chIdx,ktIdx, fitRes[chIdx][centrIdx][ktIdx]));
    }

    std::ofstream f(outDir + "/fit3d.json");
    f << std::setw(2) << J << "\n";
}

FitResult FitResultFromJson(const nlohmann::json& j) {
    FitResult r;

    if (j.at("ok").is_boolean()) {
        r.ok = j.at("ok").get<bool>();
    } else {
        r.ok = false;
    }

    if (j.contains("lambda") && j.at("lambda").is_number()) {
        r.lambda = j.at("lambda").get<double>();
    } else {
        r.lambda = 0.0;
    }

    if (j.contains("elambda") && !j.at("elambda").is_null()) {
        r.elambda = j.at("elambda").get<double>();
    } else {
        r.elambda = 0.0;
    }

    if (j.contains("chi2") && j.at("chi2").is_number()) {
        r.chi2 = j.at("chi2").get<double>();
    } else {
        r.chi2 = 0.0;
    }

    if (j.contains("ndf") && j.at("ndf").is_number()) {
        r.ndf = j.at("ndf").get<int>();
    } else {
        r.ndf = 0;
    }

    r.pvalue = j.value("pvalue", 0.0);

    if (j.contains("R") && j.at("R").is_array()) {
        size_t r_size = j.at("R").size();
        for (int i = 0; i < 6; i++) {
            if (i < r_size) {
                if (!j.at("R")[i].is_null()) {
                    r.R[i] = j.at("R")[i].get<double>();
                } else {
                    r.R[i] = 0.0;
                }
            } else {
                r.R[i] = 0.0;
            }
        }
    } else {
        for (int i = 0; i < 6; i++) {
            r.R[i] = 0.0;
        }
    }

    if (j.contains("eR") && j.at("eR").is_array()) {
        size_t er_size = j.at("eR").size();
        for (int i = 0; i < 6; i++) {
            if (i < er_size) {
                if (!j.at("eR")[i].is_null()) {
                    r.eR[i] = j.at("eR")[i].get<double>();
                } else {
                    r.eR[i] = 0.0;
                }
            } else {
                r.eR[i] = 0.0;
            }
        }
    } else {
        for (int i = 0; i < 6; i++) {
            r.eR[i] = 0.0;
        }
    }

    std::cout << "Read FitResult: R[0]=" << r.R[0] << ", lambda=" << r.lambda << ", ndf=" << r.ndf << ", chi2=" << r.chi2 << std::endl;
    return r;
}