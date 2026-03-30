#pragma once

#include <string>

#include "helpers.h"

struct Context {
    std::string inputFile;
    std::string outDir;

    Bin bining;
    FitGrid fitRes;
};

Context BuildContext(char** argv);