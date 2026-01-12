#include "analysis/registry.h"

std::vector<Stage>& GetRegistry() {
    static std::vector<Stage> r;
    return r;
}
