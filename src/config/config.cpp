#include "config/config.h"

#include <array>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "core/binning.h"
#include "core/fs.h"
#include "core/log.h"

namespace
{

void CheckObject(const nlohmann::json& j, const std::string& path)
{
    if (!j.is_object()) {
        throw std::runtime_error(path + " must be an object");
    }
}

void ForbidUnknownKeys(const nlohmann::json& j, const std::set<std::string>& allowed,
                       const std::string& path)
{
    CheckObject(j, path);

    for (const auto& item : j.items()) {
        if (allowed.count(item.key()) == 0) {
            throw std::runtime_error("unknown config field: " + path + "." + item.key());
        }
    }
}

const nlohmann::json& RequiredObject(const nlohmann::json& j, const char* key,
                                     const std::string& path)
{
    if (!j.contains(key)) {
        throw std::runtime_error("missing required field: " + path + "." + key);
    }

    if (!j.at(key).is_object()) {
        throw std::runtime_error("field must be object: " + path + "." + key);
    }

    return j.at(key);
}

bool OptionalBool(const nlohmann::json& j, const char* key, bool fallback, const std::string& path)
{
    if (!j.contains(key)) {
        return fallback;
    }

    if (!j.at(key).is_boolean()) {
        throw std::runtime_error("field must be a boolean: " + path + "." + key);
    }

    return j.at(key).get<bool>();
}
double OptionalDouble(const nlohmann::json& j, const char* key, double fallback,
                      const std::string& path)
{
    if (!j.contains(key)) {
        return fallback;
    }

    if (!j.at(key).is_number()) {
        throw std::runtime_error("field must be a number: " + path + "." + key);
    }

    return j.at(key).get<double>();
}

int OptionalInt(const nlohmann::json& j, const char* key, int fallback, const std::string& path)
{
    if (!j.contains(key)) {
        return fallback;
    }

    if (!j.at(key).is_number_integer()) {
        throw std::runtime_error("field must be an integer: " + path + "." + key);
    }

    return j.at(key).get<int>();
}
std::string RequiredString(const nlohmann::json& j, const char* key, const std::string& path)
{
    if (!j.contains(key)) {
        throw std::runtime_error("missing required field: " + path + "." + key);
    }

    if (!j.at(key).is_string()) {
        throw std::runtime_error("field must be string: " + path + "." + key);
    }

    return j.at(key).get<std::string>();
}

std::string OptionalString(const nlohmann::json& j, const char* key, const std::string& fallback,
                           const std::string& path)
{
    if (!j.contains(key)) {
        return fallback;
    }

    if (!j.at(key).is_string()) {
        throw std::runtime_error("field must be string: " + path + "." + key);
    }

    return j.at(key).get<std::string>();
}

std::vector<double> RequiredDoubleArray(const nlohmann::json& j, const char* key,
                                        const std::string& path)
{
    if (!j.contains(key)) {
        throw std::runtime_error("missing required field: " + path + "." + key);
    }

    if (!j.at(key).is_array()) {
        throw std::runtime_error("field must be an array: " + path + "." + key);
    }

    std::vector<double> out;
    for (const auto& v : j.at(key)) {
        if (!v.is_number()) {
            throw std::runtime_error("field must be an array of numbers: " + path + "." + key);
        }
        out.push_back(v.get<double>());
    }
    return out;
}

std::vector<std::string> RequiredStringArray(const nlohmann::json& j, const char* key,
                                             const std::string& path)
{
    if (!j.contains(key)) {
        throw std::runtime_error("missing required field: " + path + "." + key);
    }

    if (!j.at(key).is_array()) {
        throw std::runtime_error("field must be an array: " + path + "." + key);
    }

    std::vector<std::string> out;
    for (const auto& v : j.at(key)) {
        if (!v.is_string()) {
            throw std::runtime_error("field must be an array of strings: " + path + "." + key);
        }
        out.push_back(v.get<std::string>());
    }
    return out;
}

std::vector<std::string> OptionalStringArray(const nlohmann::json& j, const char* key,
                                             const std::string& path)
{
    if (!j.contains(key)) {
        return {};
    }

    return RequiredStringArray(j, key, path);
}

std::vector<int> OptionalIntArray(const nlohmann::json& j, const char* key,
                                  const std::vector<int>& fallback, const std::string& path)
{
    if (!j.contains(key)) {
        return fallback;
    }

    if (!j.at(key).is_array()) {
        throw std::runtime_error("field must be an array: " + path + "." + key);
    }

    std::vector<int> out;
    for (const auto& v : j.at(key)) {
        if (!v.is_number_integer()) {
            throw std::runtime_error("field must be an array of integers: " + path + "." + key);
        }
        out.push_back(v.get<int>());
    }
    return out;
}
void ParseLimits(const nlohmann::json& fit, FitConfig& out)
{
    if (!fit.contains("limits")) {
        return;
    }

    const auto& limits = RequiredObject(fit, "limits", "config.fit");
    const auto radius_sq = RequiredDoubleArray(limits, "radius_sq", "config.fit.limits");
    const auto cross = RequiredDoubleArray(limits, "cross", "config.fit.limits");
    const auto lambda = RequiredDoubleArray(limits, "lambda", "config.fit.limits");

    if (radius_sq.size() != 2 || cross.size() != 2 || lambda.size() != 2) {
        throw std::runtime_error("config.fit.limits arrays must contain exactly 2 values");
    }

    out.radius_sq_min = radius_sq[0];
    out.radius_sq_max = radius_sq[1];
    out.cross_min = cross[0];
    out.cross_max = cross[1];
    out.lambda_min = lambda[0];
    out.lambda_max = lambda[1];
}

int param_index(const std::string& name)
{
    static const std::map<std::string, int> kIndex = {{"r_out", 0}, {"r_side", 1}, {"r_long", 2},
                                                      {"r_os", 3},  {"r_ol", 4},   {"r_sl", 5},
                                                      {"lambda", 6}};
    const auto it = kIndex.find(name);
    if (it == kIndex.end()) {
        throw std::runtime_error("unknown fit parameter in freeze: " + name);
    }
    return it->second;
}

void ParseFreeze(const nlohmann::json& fit, FitConfig& out)
{
    if (!fit.contains("freeze")) {
        return;
    }

    const auto& freeze = RequiredObject(fit, "freeze", "config.fit");

    for (auto& v : out.freeze) {
        v.reset();
    }

    for (const auto& item : freeze.items()) {
        if (!item.value().is_number()) {
            throw std::runtime_error("field must be a number: config.fit.freeze." + item.key());
        }
        out.freeze[static_cast<std::size_t>(param_index(item.key()))] = item.value().get<double>();
    }
}

void ParseOptionalRange(const nlohmann::json& j, const char* key, double& min, double& max,
                        const std::string& path)
{
    if (!j.contains(key)) {
        return;
    }

    const auto range = RequiredDoubleArray(j, key, path);
    if (range.size() != 2) {
        throw std::runtime_error(path + "." + key + " must contain exactly 2 values");
    }

    min = range[0];
    max = range[1];
}

nlohmann::json ToJson(const Config& cfg)
{
    static constexpr std::array<const char*, 7> kParamNames = {"r_out", "r_side", "r_long", "r_os",
                                                               "r_ol",  "r_sl",   "lambda"};

    nlohmann::json j;
    j["machine"]["base_input"] = cfg.machine.base_input;
    j["machine"]["base_output"] = cfg.machine.base_output;
    j["vars"]["input"]["file"] = cfg.input.file;
    j["vars"]["input"]["type"] = cfg.input.type;
    j["vars"]["output"]["dir"] = cfg.output.dir;

    j["fit"]["q_max"] = cfg.fit.q_max;
    j["fit"]["use_default_ip"] = cfg.fit.use_default_ip;
    j["fit"]["options"] = cfg.fit.options;
    j["fit"]["minimizer"] = cfg.fit.minimizer;
    j["fit"]["retry_with_defaults"] = cfg.fit.retry_with_defaults;
    j["fit"]["use_integral"] = cfg.fit.use_integral;
    j["fit"]["minos_errors"] = cfg.fit.minos_errors;
    j["fit"]["limits"]["radius_sq"] = {cfg.fit.radius_sq_min, cfg.fit.radius_sq_max};
    j["fit"]["limits"]["cross"] = {cfg.fit.cross_min, cfg.fit.cross_max};
    j["fit"]["limits"]["lambda"] = {cfg.fit.lambda_min, cfg.fit.lambda_max};
    for (std::size_t i = 0; i < kParamNames.size(); ++i) {
        const auto& frozen = cfg.fit.freeze[i];
        if (frozen.has_value()) {
            j["fit"]["freeze"][kParamNames[i]] = frozen.value();
        }
    }

    j["projections"]["slice_1d"] = cfg.projections.slice_1d;
    j["projections"]["slice_2d"] = cfg.projections.slice_2d;
    j["projections"]["crop_2d"] = cfg.projections.crop_2d;
    j["projections"]["slice_ratio"] = cfg.projections.slice_ratio;
    j["projections"]["fit_over_cf_range"] = cfg.projections.fit_over_cf_range;
    j["projections"]["axis_range_1d"] = cfg.projections.axis_range_1d;
    j["projections"]["cf_y_range_1d"] = {cfg.projections.cf_y_min_1d, cfg.projections.cf_y_max_1d};
    j["projections"]["fit_over_cf_y_range_1d"] = {cfg.projections.fit_over_cf_y_min_1d,
                                                  cfg.projections.fit_over_cf_y_max_1d};

    j["binning"]["values"] = cfg.binning.values;
    j["binning"]["names"] = cfg.binning.names;
    j["binning"]["file_names"] = cfg.binning.file_names;

    j["images"]["need"] = cfg.general.images.need;
    j["images"]["format"] = cfg.general.images.format;

    j["stages"]["cf3d"] = cfg.stages.cf3d;
    j["stages"]["dependency"] = cfg.stages.dependency;
    j["stages"]["projections_1d"] = cfg.stages.projections_1d;
    j["stages"]["projections_2d"] = cfg.stages.projections_2d;
    j["stages"]["ratios"] = cfg.stages.ratios;

    j["selection"]["charges"] = cfg.selection.charges;
    j["selection"]["centralities"] = cfg.selection.centralities;

    j["logging"]["level"] = cfg.logging.level;
    j["logging"]["file"] = cfg.logging.file;

    j["threads"] = cfg.threads;

    return j;
}
void Validate(const Config& cfg)
{
    if (cfg.machine.base_input.empty()) {
        throw std::runtime_error("machine.base_input must not be empty");
    }

    if (cfg.machine.base_output.empty()) {
        throw std::runtime_error("machine.base_output must not be empty");
    }

    if (cfg.input.file.empty()) {
        throw std::runtime_error("vars.input.file must not be empty");
    }

    if (cfg.output.dir.empty()) {
        throw std::runtime_error("vars.output.dir must not be empty");
    }

    if (cfg.input.type != "kt" && cfg.input.type != "rapidity") {
        throw std::runtime_error("vars.input.type must be either \"kt\" or \"rapidity\"");
    }
    if (cfg.fit.q_max <= 0) {
        throw std::runtime_error("fit.q_max must be positive");
    }

    if (cfg.fit.minimizer.empty()) {
        throw std::runtime_error("fit.minimizer must not be empty");
    }

    if (cfg.fit.radius_sq_min >= cfg.fit.radius_sq_max) {
        throw std::runtime_error("fit.limits.radius_sq: min must be less than max");
    }

    if (cfg.fit.cross_min >= cfg.fit.cross_max) {
        throw std::runtime_error("fit.limits.cross: min must be less than max");
    }

    if (cfg.fit.lambda_min >= cfg.fit.lambda_max) {
        throw std::runtime_error("fit.limits.lambda: min must be less than max");
    }
    if (cfg.general.images.format.empty()) {
        throw std::runtime_error("images.format must not be empty");
    }

    if (cfg.projections.axis_range_1d <= 0) {
        throw std::runtime_error("projections.axis_range_1d must be positive");
    }

    if (cfg.projections.cf_y_min_1d >= cfg.projections.cf_y_max_1d) {
        throw std::runtime_error("projections.cf_y_range_1d: min must be less than max");
    }

    if (cfg.projections.fit_over_cf_y_min_1d >= cfg.projections.fit_over_cf_y_max_1d) {
        throw std::runtime_error("projections.fit_over_cf_y_range_1d: min must be less than max");
    }

    if (!cfg.binning.values.empty()) {
        if (cfg.binning.values.size() < 2) {
            throw std::runtime_error("binning.values must contain at least 2 edges");
        }

        const std::size_t nBins = cfg.binning.values.size() - 1;
        if (cfg.binning.names.size() != nBins) {
            throw std::runtime_error("binning.names must contain values.size() - 1 entries");
        }

        if (!cfg.binning.file_names.empty() && cfg.binning.file_names.size() != nBins) {
            throw std::runtime_error("binning.file_names must contain values.size() - 1 entries");
        }
    }

    if (cfg.selection.charges.empty()) {
        throw std::runtime_error("selection.charges must not be empty");
    }

    for (const int ch : cfg.selection.charges) {
        if (ch < 0 || ch >= charge::kCount) {
            throw std::runtime_error("selection.charges values must be in [0, " +
                                     std::to_string(charge::kCount) + ")");
        }
    }

    if (cfg.selection.centralities.empty()) {
        throw std::runtime_error("selection.centralities must not be empty");
    }

    for (const int centr : cfg.selection.centralities) {
        if (centr < 0 || centr >= centrality::kCount) {
            throw std::runtime_error("selection.centralities values must be in [0, " +
                                     std::to_string(centrality::kCount) + ")");
        }
    }
    if (cfg.threads < 0) {
        throw std::runtime_error("threads must be >= 0 (0 = auto)");
    }

    static_cast<void>(log::parse_level(cfg.logging.level));
}

} // namespace

void build(Config& cfg)
{
    if (cfg.binning.values.empty()) {
        if (cfg.input.type == "kt") {
            cfg.binning.count = kt::kCount;
            cfg.binning.names.assign(kt::kNames.begin(), kt::kNames.end());
            cfg.binning.values.assign(kt::kValues.begin(), kt::kValues.end());
            cfg.binning.file_names.assign(kt::kFileNames.begin(), kt::kFileNames.end());
        } else if (cfg.input.type == "rapidity") {
            cfg.binning.count = rapidity::kCount;
            cfg.binning.names.assign(rapidity::kNames.begin(), rapidity::kNames.end());
            cfg.binning.values.assign(rapidity::kValues.begin(), rapidity::kValues.end());
            cfg.binning.file_names.assign(rapidity::kFileNames.begin(), rapidity::kFileNames.end());
        }
    } else {
        cfg.binning.count = static_cast<int>(cfg.binning.values.size()) - 1;
        if (cfg.binning.file_names.empty()) {
            cfg.binning.file_names = cfg.binning.names;
        }
    }

    cfg.fit_results =
        FitGrid(charge::kCount, std::vector<std::vector<FitResult>>(
                                    centrality::kCount, std::vector<FitResult>(cfg.binning.count)));
}

void StrictConfig(const nlohmann::json& raw)
{
    const auto& machine = RequiredObject(raw, "machine", "config");
    const auto& vars = RequiredObject(raw, "vars", "config");
    const auto& input = RequiredObject(vars, "input", "config.vars");
    const auto& output = RequiredObject(vars, "output", "config.vars");

    ForbidUnknownKeys(raw,
                      {"machine", "vars", "fit", "projections", "binning", "images", "stages",
                       "selection", "threads", "logging"},
                      "config");
    ForbidUnknownKeys(machine, {"base_input", "base_output"}, "config.machine");
    ForbidUnknownKeys(vars, {"input", "output"}, "config.vars");
    ForbidUnknownKeys(input, {"file", "type"}, "config.vars.input");
    ForbidUnknownKeys(output, {"dir"}, "config.vars.output");
    if (raw.contains("fit")) {
        const auto& fit = RequiredObject(raw, "fit", "config");
        ForbidUnknownKeys(fit,
                          {"q_max", "use_default_ip", "options", "limits", "freeze", "minimizer",
                           "retry_with_defaults", "use_integral", "minos_errors"},
                          "config.fit");
        if (fit.contains("limits")) {
            const auto& limits = RequiredObject(fit, "limits", "config.fit");
            ForbidUnknownKeys(limits, {"radius_sq", "cross", "lambda"}, "config.fit.limits");
        }
        if (fit.contains("freeze")) {
            ForbidUnknownKeys(RequiredObject(fit, "freeze", "config.fit"),
                              {"r_out", "r_side", "r_long", "r_os", "r_ol", "r_sl", "lambda"},
                              "config.fit.freeze");
        }
    }
    if (raw.contains("projections")) {
        ForbidUnknownKeys(RequiredObject(raw, "projections", "config"),
                          {"slice_1d", "slice_2d", "crop_2d", "slice_ratio", "fit_over_cf_range",
                           "axis_range_1d", "cf_y_range_1d", "fit_over_cf_y_range_1d"},
                          "config.projections");
    }
    if (raw.contains("binning")) {
        ForbidUnknownKeys(RequiredObject(raw, "binning", "config"),
                          {"values", "names", "file_names"}, "config.binning");
    }

    if (raw.contains("images")) {
        ForbidUnknownKeys(RequiredObject(raw, "images", "config"), {"need", "format"},
                          "config.images");
    }

    if (raw.contains("stages")) {
        ForbidUnknownKeys(RequiredObject(raw, "stages", "config"),
                          {"cf3d", "dependency", "projections_1d", "projections_2d", "ratios"},
                          "config.stages");
    }
    if (raw.contains("selection")) {
        ForbidUnknownKeys(RequiredObject(raw, "selection", "config"), {"charges", "centralities"},
                          "config.selection");
    }

    if (raw.contains("logging")) {
        ForbidUnknownKeys(RequiredObject(raw, "logging", "config"), {"level", "file"},
                          "config.logging");
    }
}

Config load(const std::string& path)
{
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open config: " + path);
    }

    nlohmann::json raw;
    try {
        in >> raw;
    } catch (const nlohmann::json::parse_error& e) {
        throw std::runtime_error("failed to parse config " + path + ": " + std::string(e.what()));
    }

    const auto& machine = RequiredObject(raw, "machine", "config");
    const auto& vars = RequiredObject(raw, "vars", "config");
    const auto& input = RequiredObject(vars, "input", "config.vars");
    const auto& output = RequiredObject(vars, "output", "config.vars");
    StrictConfig(raw);

    Config cfg;

    cfg.machine.base_input = RequiredString(machine, "base_input", "config.machine");
    cfg.machine.base_output = RequiredString(machine, "base_output", "config.machine");

    cfg.input.file =
        cfg.machine.base_input + "/" + RequiredString(input, "file", "config.vars.input");

    cfg.input.type = RequiredString(input, "type", "config.vars.input");

    cfg.output.dir = cfg.machine.base_output + "/" +
                     OptionalString(output, "dir", "results", "config.vars.output");
    if (raw.contains("fit")) {
        const auto& fit = RequiredObject(raw, "fit", "config");
        cfg.fit.q_max = OptionalDouble(fit, "q_max", cfg.fit.q_max, "config.fit");
        cfg.fit.use_default_ip =
            OptionalBool(fit, "use_default_ip", cfg.fit.use_default_ip, "config.fit");
        cfg.fit.options = OptionalString(fit, "options", cfg.fit.options, "config.fit");
        cfg.fit.minimizer = OptionalString(fit, "minimizer", cfg.fit.minimizer, "config.fit");
        cfg.fit.retry_with_defaults =
            OptionalBool(fit, "retry_with_defaults", cfg.fit.retry_with_defaults, "config.fit");
        cfg.fit.use_integral = OptionalBool(fit, "use_integral", cfg.fit.use_integral, "config.fit");
        cfg.fit.minos_errors = OptionalBool(fit, "minos_errors", cfg.fit.minos_errors, "config.fit");
        ParseLimits(fit, cfg.fit);
        ParseFreeze(fit, cfg.fit);
    }
    if (raw.contains("projections")) {
        const auto& proj = RequiredObject(raw, "projections", "config");
        cfg.projections.slice_1d =
            OptionalDouble(proj, "slice_1d", cfg.projections.slice_1d, "config.projections");
        cfg.projections.slice_2d =
            OptionalDouble(proj, "slice_2d", cfg.projections.slice_2d, "config.projections");
        cfg.projections.crop_2d =
            OptionalDouble(proj, "crop_2d", cfg.projections.crop_2d, "config.projections");
        cfg.projections.slice_ratio =
            OptionalDouble(proj, "slice_ratio", cfg.projections.slice_ratio, "config.projections");
        cfg.projections.fit_over_cf_range = OptionalDouble(
            proj, "fit_over_cf_range", cfg.projections.fit_over_cf_range, "config.projections");
        cfg.projections.axis_range_1d = OptionalDouble(
            proj, "axis_range_1d", cfg.projections.axis_range_1d, "config.projections");
        ParseOptionalRange(proj, "cf_y_range_1d", cfg.projections.cf_y_min_1d,
                           cfg.projections.cf_y_max_1d, "config.projections");
        ParseOptionalRange(proj, "fit_over_cf_y_range_1d", cfg.projections.fit_over_cf_y_min_1d,
                           cfg.projections.fit_over_cf_y_max_1d, "config.projections");
    }

    if (raw.contains("binning")) {
        const auto& binning = RequiredObject(raw, "binning", "config");
        cfg.binning.values = RequiredDoubleArray(binning, "values", "config.binning");
        cfg.binning.names = RequiredStringArray(binning, "names", "config.binning");
        cfg.binning.file_names = OptionalStringArray(binning, "file_names", "config.binning");
    }

    if (raw.contains("images")) {
        const auto& images = RequiredObject(raw, "images", "config");
        cfg.general.images.need =
            OptionalBool(images, "need", cfg.general.images.need, "config.images");
        cfg.general.images.format =
            OptionalString(images, "format", cfg.general.images.format, "config.images");
    }

    if (raw.contains("stages")) {
        const auto& stages = RequiredObject(raw, "stages", "config");
        cfg.stages.cf3d = OptionalBool(stages, "cf3d", cfg.stages.cf3d, "config.stages");
        cfg.stages.dependency =
            OptionalBool(stages, "dependency", cfg.stages.dependency, "config.stages");
        cfg.stages.projections_1d =
            OptionalBool(stages, "projections_1d", cfg.stages.projections_1d, "config.stages");
        cfg.stages.projections_2d =
            OptionalBool(stages, "projections_2d", cfg.stages.projections_2d, "config.stages");
        cfg.stages.ratios = OptionalBool(stages, "ratios", cfg.stages.ratios, "config.stages");
    }

    if (raw.contains("selection")) {
        const auto& selection = RequiredObject(raw, "selection", "config");
        cfg.selection.charges =
            OptionalIntArray(selection, "charges", cfg.selection.charges, "config.selection");
        cfg.selection.centralities = OptionalIntArray(
            selection, "centralities", cfg.selection.centralities, "config.selection");
    }

    cfg.threads = OptionalInt(raw, "threads", cfg.threads, "config");

    if (raw.contains("logging")) {
        const auto& logging = RequiredObject(raw, "logging", "config");
        cfg.logging.level = OptionalString(logging, "level", cfg.logging.level, "config.logging");
        cfg.logging.file = OptionalString(logging, "file", cfg.logging.file, "config.logging");
    }

    Validate(cfg);
    build(cfg);

    ensure_dir(cfg.output.dir);

    const std::string snapshotPath = cfg.output.dir + "/run_config.json";
    std::ofstream snapshot(snapshotPath);
    if (snapshot) {
        snapshot << ToJson(cfg).dump(2);
    }

    return cfg;
}
