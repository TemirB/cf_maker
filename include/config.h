#pragma once

#include <string>

#include <helpers.h>

struct Machine {
    std::string base_input;
    std::string base_output;
};

struct Input {
    std::string type;
    std::string file;
};

struct Output {
    std::string dir;
};

struct Config {
    Machine machine;

    Input input;
    Output output;

    Bin bining;
    FitGrid fitRes;
};

Config Load(const std::string& path);
void Build(Config& cfg);