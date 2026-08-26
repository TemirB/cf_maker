#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>

#include "config/config.h"

namespace
{

#define CHECK(cond)                                                                                \
    do {                                                                                           \
        if (!(cond)) {                                                                             \
            std::cerr << "FAILED: " #cond " at " << __FILE__ << ":" << __LINE__ << "\n";           \
            std::exit(1);                                                                          \
        }                                                                                          \
    } while (0)

#define CHECK_THROWS(expr)                                                                         \
    do {                                                                                           \
        bool thrown = false;                                                                       \
        try {                                                                                      \
            (void)(expr);                                                                          \
        } catch (const std::exception&) {                                                          \
            thrown = true;                                                                         \
        }                                                                                          \
        if (!thrown) {                                                                             \
            std::cerr << "FAILED: expected exception at " << __FILE__ << ":" << __LINE__ << "\n";  \
            std::exit(1);                                                                          \
        }                                                                                          \
    } while (0)

std::filesystem::path TmpDir()
{
    const auto dir = std::filesystem::temp_directory_path() / "cf_maker_test_config";
    std::filesystem::create_directories(dir);
    return dir;
}

void WriteConfig(const std::filesystem::path& path, const std::string& content)
{
    std::ofstream out(path);
    out << content;
    out.close();
    CHECK(std::filesystem::exists(path));
}

const std::string kBase = TmpDir().string();

std::string Wrap(const std::string& extra)
{
    std::string json = "{\n"
                       "    \"machine\": {\n"
                       "        \"base_input\": \"" +
                       kBase +
                       "/input\",\n"
                       "        \"base_output\": \"" +
                       kBase +
                       "/output\"\n"
                       "    },\n"
                       "    \"vars\": {\n"
                       "        \"input\": {\n"
                       "            \"file\": \"merged.root\",\n"
                       "            \"type\": \"kt\"\n"
                       "        },\n"
                       "        \"output\": {\n"
                       "            \"dir\": \"results\"\n"
                       "        }\n"
                       "    }";
    if (!extra.empty()) {
        json += ",\n    " + extra;
    }
    json += "\n}\n";
    return json;
}

void TestDefaults()
{
    const auto path = TmpDir() / "defaults.json";
    WriteConfig(path, Wrap(""));

    const Config cfg = load(path.string());

    CHECK(cfg.fit.q_max == 0.20);
    CHECK(!cfg.fit.use_default_ip);
    CHECK(cfg.fit.options == "RQS0");
    CHECK(cfg.fit.minimizer == "Minuit2");
    CHECK(cfg.fit.retry_with_defaults);
    CHECK(!cfg.fit.use_integral);
    CHECK(!cfg.fit.minos_errors);
    CHECK(cfg.fit.radius_sq_min == 0.0);
    CHECK(cfg.fit.radius_sq_max == 100.0);
    CHECK(cfg.fit.cross_min == -30.0);
    CHECK(cfg.fit.lambda_max == 1.0);

    CHECK(!cfg.fit.freeze[0].has_value());
    CHECK(!cfg.fit.freeze[1].has_value());
    CHECK(!cfg.fit.freeze[2].has_value());
    CHECK(cfg.fit.freeze[3].has_value() && *cfg.fit.freeze[3] == 0.0);
    CHECK(cfg.fit.freeze[4].has_value() && *cfg.fit.freeze[4] == 0.0);
    CHECK(cfg.fit.freeze[5].has_value() && *cfg.fit.freeze[5] == 0.0);
    CHECK(!cfg.fit.freeze[6].has_value());

    CHECK(cfg.projections.slice_1d == 0.05);
    CHECK(cfg.projections.slice_2d == 0.2);
    CHECK(cfg.projections.crop_2d == 0.1);
    CHECK(cfg.projections.slice_ratio == 0.05);
    CHECK(cfg.projections.fit_over_cf_range == 0.08);
    CHECK(cfg.projections.axis_range_1d == 0.2);
    CHECK(cfg.projections.cf_y_min_1d == 0.9);
    CHECK(cfg.projections.cf_y_max_1d == 1.7);
    CHECK(cfg.projections.fit_over_cf_y_min_1d == 0.9);
    CHECK(cfg.projections.fit_over_cf_y_max_1d == 1.1);

    CHECK(cfg.general.images.need);
    CHECK(cfg.general.images.format == "pdf");

    CHECK(cfg.stages.cf3d && cfg.stages.dependency && cfg.stages.projections_1d &&
          cfg.stages.projections_2d && cfg.stages.ratios);

    CHECK(cfg.selection.charges == std::vector<int>({0, 1}));
    CHECK(cfg.selection.centralities == std::vector<int>({0, 1, 2, 3}));

    CHECK(cfg.logging.level == "info");
    CHECK(cfg.logging.file == "run.log");
    CHECK(cfg.threads == 0);

    CHECK(cfg.binning.count == kt::kCount);
    CHECK(cfg.binning.names.size() == static_cast<std::size_t>(kt::kCount));
    CHECK(cfg.binning.names[0] == "0.15-0.25");
    CHECK(cfg.input.file == kBase + "/input/merged.root");
    CHECK(cfg.output.dir == kBase + "/output/results");
}

void TestOverrides()
{
    const auto path = TmpDir() / "overrides.json";
    WriteConfig(path, Wrap("\"fit\": {\n"
                           "            \"q_max\": 0.30,\n"
                           "            \"use_default_ip\": true,\n"
                           "            \"options\": \"R\",\n"
                           "            \"minimizer\": \"Minuit\",\n"
                           "            \"retry_with_defaults\": false,\n"
                           "            \"use_integral\": true,\n"
                           "            \"minos_errors\": true,\n"
                           "            \"limits\": {\n"
                           "                \"radius_sq\": [1.0, 50.0],\n"
                           "                \"cross\": [-5.0, 5.0],\n"
                           "                \"lambda\": [0.1, 0.9]\n"
                           "            },\n"
                           "            \"freeze\": {\n"
                           "                \"r_out\": 5.0,\n"
                           "                \"lambda\": 0.6\n"
                           "            }\n"
                           "        },\n"
                           "        \"projections\": {\n"
                           "            \"slice_1d\": 0.1,\n"
                           "            \"fit_over_cf_range\": 0.15,\n"
                           "            \"axis_range_1d\": 0.3,\n"
                           "            \"cf_y_range_1d\": [0.8, 1.6]\n"
                           "        },\n"
                           "        \"binning\": {\n"
                           "            \"values\": [0.0, 0.5, 1.0],\n"
                           "            \"names\": [\"lo\", \"hi\"]\n"
                           "        },\n"
                           "        \"images\": {\n"
                           "            \"need\": false,\n"
                           "            \"format\": \"png\"\n"
                           "        },\n"
                           "        \"stages\": {\n"
                           "            \"cf3d\": false,\n"
                           "            \"ratios\": false\n"
                           "        },\n"
                           "        \"selection\": {\n"
                           "            \"charges\": [1],\n"
                           "            \"centralities\": [0, 2]\n"
                           "        },\n"
                           "        \"logging\": {\n"
                           "            \"level\": \"debug\",\n"
                           "            \"file\": \"\"\n"
                           "        }"));

    const Config cfg = load(path.string());

    CHECK(cfg.fit.q_max == 0.30);
    CHECK(cfg.fit.use_default_ip);
    CHECK(cfg.fit.options == "R");
    CHECK(cfg.fit.minimizer == "Minuit");
    CHECK(!cfg.fit.retry_with_defaults);
    CHECK(cfg.fit.use_integral);
    CHECK(cfg.fit.minos_errors);
    CHECK(cfg.fit.radius_sq_min == 1.0);
    CHECK(cfg.fit.radius_sq_max == 50.0);
    CHECK(cfg.fit.cross_min == -5.0);
    CHECK(cfg.fit.lambda_max == 0.9);

    CHECK(cfg.fit.freeze[0].has_value() && *cfg.fit.freeze[0] == 5.0);
    CHECK(cfg.fit.freeze[6].has_value() && *cfg.fit.freeze[6] == 0.6);
    CHECK(!cfg.fit.freeze[3].has_value());
    CHECK(!cfg.fit.freeze[4].has_value());
    CHECK(!cfg.fit.freeze[5].has_value());

    CHECK(cfg.projections.slice_1d == 0.1);
    CHECK(cfg.projections.fit_over_cf_range == 0.15);
    CHECK(cfg.projections.slice_2d == 0.2);
    CHECK(cfg.projections.axis_range_1d == 0.3);
    CHECK(cfg.projections.cf_y_min_1d == 0.8);
    CHECK(cfg.projections.cf_y_max_1d == 1.6);

    CHECK(cfg.binning.count == 2);
    CHECK(cfg.binning.values == std::vector<double>({0.0, 0.5, 1.0}));
    CHECK(cfg.binning.names[0] == "lo");
    CHECK(cfg.binning.file_names == cfg.binning.names);

    CHECK(!cfg.general.images.need);
    CHECK(cfg.general.images.format == "png");

    CHECK(!cfg.stages.cf3d);
    CHECK(cfg.stages.dependency);
    CHECK(!cfg.stages.ratios);

    CHECK(cfg.selection.charges == std::vector<int>({1}));
    CHECK(cfg.selection.centralities == std::vector<int>({0, 2}));

    CHECK(cfg.logging.level == "debug");
    CHECK(cfg.logging.file.empty());
}

void TestSnapshot()
{
    const auto path = TmpDir() / "snapshot.json";
    WriteConfig(path, Wrap(""));

    const Config cfg = load(path.string());

    const std::filesystem::path snapshotPath =
        std::filesystem::path(cfg.output.dir) / "run_config.json";
    CHECK(std::filesystem::exists(snapshotPath));

    std::ifstream in(snapshotPath);
    CHECK(in.good());
    const std::string content((std::istreambuf_iterator<char>(in)),
                              std::istreambuf_iterator<char>());
    CHECK(content.find("\"threads\"") != std::string::npos);
    CHECK(content.find("\"q_max\"") != std::string::npos);
    CHECK(content.find("\"r_os\"") != std::string::npos);
    CHECK(content.find("\"binning\"") != std::string::npos);
}

void TestErrors()
{
    const auto dir = TmpDir();

    const auto unknownKey = dir / "unknown_key.json";
    WriteConfig(unknownKey, Wrap("\"foo\": 1"));
    CHECK_THROWS(load(unknownKey.string()));

    const auto badType = dir / "bad_type.json";
    WriteConfig(badType, "{\"machine\": {\"base_input\": \"/in\", \"base_output\": \"/out\"}, "
                         "\"vars\": {\"input\": {\"file\": \"f.root\", \"type\": \"nope\"}, "
                         "\"output\": {\"dir\": \"r\"}}}");
    CHECK_THROWS(load(badType.string()));

    const auto badLimits = dir / "bad_limits.json";
    WriteConfig(badLimits, Wrap("\"fit\": {\"limits\": {\"lambda\": [0.5]}}"));
    CHECK_THROWS(load(badLimits.string()));

    const auto badFreeze = dir / "bad_freeze.json";
    WriteConfig(badFreeze, Wrap("\"fit\": {\"freeze\": {\"r_foo\": 1.0}}"));
    CHECK_THROWS(load(badFreeze.string()));

    const auto badBinning = dir / "bad_binning.json";
    WriteConfig(badBinning, Wrap("\"binning\": {\"values\": [0.0, 1.0, 2.0], \"names\": [\"a\"]}"));
    CHECK_THROWS(load(badBinning.string()));

    const auto badSelection = dir / "bad_selection.json";
    WriteConfig(badSelection, Wrap("\"selection\": {\"charges\": [5]}"));
    CHECK_THROWS(load(badSelection.string()));

    const auto badCharge = dir / "bad_charge.json";
    WriteConfig(badCharge, Wrap("\"selection\": {\"centralities\": [-1]}"));
    CHECK_THROWS(load(badCharge.string()));

    const auto badLogLevel = dir / "bad_log_level.json";
    WriteConfig(badLogLevel, Wrap("\"logging\": {\"level\": \"verbose\"}"));
    CHECK_THROWS(load(badLogLevel.string()));
}

} // namespace

int main()
{
    TestDefaults();
    TestOverrides();
    TestSnapshot();
    TestErrors();
    std::cout << "All tests passed\n";
    return 0;
}
