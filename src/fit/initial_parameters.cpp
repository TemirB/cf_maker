#include <fit/initial_parameters.h>

InitialParameters::InitialParameters() {
    data_ = {
        {
            "out", {
                {0, {6.2, 5.9, 5.7, 5.4}},
                {1, {5.35, 5.15, 4.9, 4.75}},
                {2, {4.6, 4.4, 4.2, 4.0}},
                {3, {4.2, 4.0, 3.8, 3.6}}
            }
        },
        {
            "side", {
                {0, {4.8, 4.1, 3.8, 3.6}},
                {1, {4.1, 3.6, 3.15, 2.8}},
                {2, {3.4, 3.0, 2.8, 2.6}},
                {3, {2.8, 2.5, 2.3, 2.15}}
            }
        },
        {
            "long", {
                {0, {4.5, 3.5, 3.0, 2.8}},
                {1, {4.3, 3.35, 2.85, 2.6}},
                {2, {3.67, 2.86, 2.5, 2.3}},
                {3, {2.8, 2.5, 2.3, 2.15}}
            }
        },
        {
            "lambda", {
                {0, {0.835, 0.867, 0.944, 0.96}}, // 3
                {1, {0.819, 0.855, 0.900, 0.93}}, // 23
                {2, {0.799, 0.830, 0.870, 0.92}}, // 123
                {3, {0.785, 0.800, 0.850, 0.87}} // 0123
            }
        }
    };
}

double InitialParameters::get(std::string parameter, int centrality, int kt) {
    auto centralityMap = data_.at(parameter).find(centrality);
    const auto& vec = centralityMap->second;
    
    return vec[kt];
}