#include "config.h"

#include <fstream>
#include <set>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include "helpers.h"

namespace {

void CheckObject(const nlohmann::json& j, const std::string& path) {
    if (!j.is_object()) {
        throw std::runtime_error(path + " must be an object");
    }
}

void ForbidUnknownKeys(
    const nlohmann::json& j,
    const std::set<std::string>& allowed,
    const std::string& path
) {
    CheckObject(j, path);

    for (const auto& item : j.items()) {
        if (!allowed.count(item.key())) {
            throw std::runtime_error(
                "unknown config field: " + path + "." + item.key()
            );
        }
    }
}

const nlohmann::json& RequiredObject(
    const nlohmann::json& j,
    const char* key,
    const std::string& path
) {
    if (!j.contains(key)) {
        throw std::runtime_error("missing required field: " + path + "." + key);
    }

    if (!j.at(key).is_object()) {
        throw std::runtime_error("field must be object: " + path + "." + key);
    }

    return j.at(key);
}

std::string RequiredString(
    const nlohmann::json& j,
    const char* key,
    const std::string& path
) {
    if (!j.contains(key)) {
        throw std::runtime_error("missing required field: " + path + "." + key);
    }

    if (!j.at(key).is_string()) {
        throw std::runtime_error("field must be string: " + path + "." + key);
    }

    return j.at(key).get<std::string>();
}

std::string OptionalString(
    const nlohmann::json& j,
    const char* key,
    const std::string& fallback,
    const std::string& path
) {
    if (!j.contains(key)) {
        return fallback;
    }

    if (!j.at(key).is_string()) {
        throw std::runtime_error("field must be string: " + path + "." + key);
    }

    return j.at(key).get<std::string>();
}

void Validate(const Config& cfg) {
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
        throw std::runtime_error(
            "vars.input.type must be either \"kt\" or \"rapidity\""
        );
    }
}

} // namespace

void Build(Config& cfg) {
    if (cfg.input.type == "kt") {
        cfg.bining.count = Kt::kCount;
        cfg.bining.names = Kt::kNames;
        cfg.bining.values = Kt::kValues;
        cfg.bining.fileNames = Kt::kFileNames;
    } else if (cfg.input.type == "rapidity") {
        cfg.bining.count = Rapidity::kCount;
        cfg.bining.names = Rapidity::kNames;
        cfg.bining.values = Rapidity::kValues;
        cfg.bining.fileNames = Rapidity::kFileNames;
    }

    cfg.fitRes = FitGrid(
        Charge::kCount,
        std::vector<std::vector<FitResult>>(
            Centrality::kCount,
            std::vector<FitResult>(cfg.bining.count)
        )
    );
}

Config Load(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open config: " + path);
    }

    nlohmann::json j;
    try {
        in >> j;
    } catch (const nlohmann::json::parse_error& e) {
        throw std::runtime_error(
            "failed to parse config " + path + ": " + std::string(e.what())
        );
    }

    ForbidUnknownKeys(j, {"machine", "vars"}, "config");

    const auto& machine = RequiredObject(j, "machine", "config");
    const auto& vars = RequiredObject(j, "vars", "config");

    ForbidUnknownKeys(
        machine,
        {"base_input", "base_output"},
        "config.machine"
    );

    ForbidUnknownKeys(
        vars,
        {"input", "output"},
        "config.vars"
    );

    const auto& input = RequiredObject(vars, "input", "config.vars");
    const auto& output = RequiredObject(vars, "output", "config.vars");

    ForbidUnknownKeys(
        input,
        {"file", "type"},
        "config.vars.input"
    );

    ForbidUnknownKeys(
        output,
        {"dir"},
        "config.vars.output"
    );

    Config cfg;

    cfg.machine.base_input = RequiredString(
        machine,
        "base_input",
        "config.machine"
    );

    cfg.machine.base_output = RequiredString(
        machine,
        "base_output",
        "config.machine"
    );

    cfg.input.file =
        cfg.machine.base_input + "/" +
        RequiredString(input, "file", "config.vars.input");

    cfg.input.type = RequiredString(
        input,
        "type",
        "config.vars.input"
    );

    cfg.output.dir =
        cfg.machine.base_output + "/" +
        OptionalString(output, "dir", "results", "config.vars.output");

    Validate(cfg);
    Build(cfg);

    EnsureDir(cfg.output.dir);

    return cfg;
}