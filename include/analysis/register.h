#pragma once

#include "analysis/registry.h"

struct StageAutoRegister {
    StageAutoRegister(const char* name, StageFn fn) {
        GetRegistry().push_back({name, fn});
    }
};

#define REGISTER_STAGE(NAME, FN) \
    static StageAutoRegister _reg_##__COUNTER__(NAME, FN);
