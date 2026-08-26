#pragma once

#include <array>
#include <optional>
#include <string>
#include <vector>

#include "core/binning.h"
#include "fit/types.h"

struct Machine
{
    std::string base_input;
    std::string base_output;
};

struct Input
{
    std::string type;
    std::string file;
};

struct Output
{
    std::string dir;
};

struct Images
{
    bool need = true;
    std::string format = "pdf";
};

struct General
{
    Images images;
};

struct FitConfig
{
    double q_max = 0.20;
    bool use_default_ip = false;
    double radius_sq_min = 0.0;
    double radius_sq_max = 100.0;
    double cross_min = -30.0;
    double cross_max = 30.0;
    double lambda_min = 0.3;
    double lambda_max = 1.0;
    std::string options = "RQS0";
    std::string minimizer = "Minuit2";
    bool retry_with_defaults = true;
    bool use_integral = false;
    bool minos_errors = false;
    std::array<std::optional<double>, 7> freeze = {std::nullopt, std::nullopt, std::nullopt, 0.0,
                                                   0.0,          0.0,          std::nullopt};
};

struct ProjectionConfig
{
    double slice_1d = 0.05;
    double slice_2d = 0.2;
    double crop_2d = 0.1;
    double slice_ratio = 0.05;
    double fit_over_cf_range = 0.08;
    double axis_range_1d = 0.2;
    double cf_y_min_1d = 0.9;
    double cf_y_max_1d = 1.7;
    double fit_over_cf_y_min_1d = 0.9;
    double fit_over_cf_y_max_1d = 1.1;
};

struct Stages
{
    bool cf3d = true;
    bool dependency = true;
    bool projections_1d = true;
    bool projections_2d = true;
    bool ratios = true;
};

struct Selection
{
    std::vector<int> charges = {0, 1};
    std::vector<int> centralities = {0, 1, 2, 3};
};

struct LoggingConfig
{
    std::string level = "info";
    std::string file = "run.log";
};

struct Config
{
    General general;
    Machine machine;

    Input input;
    Output output;

    FitConfig fit;
    ProjectionConfig projections;
    Stages stages;
    Selection selection;
    LoggingConfig logging;

    int threads = 0;

    Bin binning;
    FitGrid fit_results;
};

[[nodiscard]] Config load(const std::string& path);
void build(Config& cfg);
