#pragma once

#include <functional>
#include <string>

#include "analysis/context.h"

using StageFn = std::function<void(RunContext&)>;

struct Stage {
    std::string name;
    StageFn fn;
};
