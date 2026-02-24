#include <fit/initial_parameters.h>

InitialParameters::InitialParameters() {
    // ПРИНЦИП: centrality 0 (most central) → largest radii
    //          centrality 3 (most peripheral) → smallest radii
    // Значения, отмеченные "// ESTIMATED", были скорректированы для соблюдения тренда
    
    // pData_ : charge = 0 (positive)
    pData_ = {
        {
            "out", {
                {0, {6.76, 6.34, 6.05, 6.05, 5.93, 5.99, 6.00, 6.20, 6.28, 6.79}},  // centrality 0: reference
                {1, {6.39, 5.76, 5.40, 5.19, 5.08, 5.12, 5.22, 5.34, 5.74, 6.33}},  // centrality 1: ~10% smaller
                {2, {5.20, 4.90, 4.60, 4.40, 4.30, 4.35, 4.45, 4.70, 5.10, 5.50}},  // ESTIMATED: centrality 2: ~20% smaller, smoothed
                {3, {4.20, 4.00, 3.85, 3.75, 3.70, 3.75, 3.85, 4.05, 4.30, 4.60}}   // ESTIMATED: centrality 3: ~30% smaller, smoothed
            }
        },
        {
            "side", {
                {0, {4.40, 4.43, 4.44, 4.51, 4.55, 4.52, 4.50, 4.48, 4.41, 4.45}},  // centrality 0: reference
                {1, {3.99, 3.94, 3.92, 3.88, 3.87, 3.87, 3.86, 3.89, 3.92, 3.90}},  // centrality 1: ~10% smaller
                {2, {3.60, 3.50, 3.45, 3.40, 3.42, 3.50, 3.48, 3.55, 3.65, 3.75}},  // ESTIMATED: smoothed, removed outliers (9.99, 8.35)
                {3, {3.20, 3.15, 3.10, 3.05, 3.08, 3.15, 3.20, 3.30, 3.40, 3.50}}   // ESTIMATED: peripheral, smooth trend
            }
        },
        {
            "long", {
                {0, {5.50, 4.80, 4.29, 4.06, 3.91, 3.93, 4.06, 4.34, 4.75, 5.61}},  // centrality 0: reference
                {1, {5.18, 4.55, 4.10, 3.86, 3.74, 3.72, 3.87, 4.12, 4.53, 5.13}},  // centrality 1: ~10% smaller
                {2, {4.50, 4.10, 3.75, 3.55, 3.45, 3.48, 3.60, 3.80, 4.15, 4.70}},  // ESTIMATED: removed outliers (9.95, 9.99)
                {3, {3.80, 3.60, 3.45, 3.35, 3.30, 3.35, 3.45, 3.65, 3.90, 4.20}}   // ESTIMATED: smooth peripheral trend
            }
        },
        {
            "lambda", {
                {0, {0.886, 0.831, 0.803, 0.801, 0.799, 0.801, 0.801, 0.819, 0.821, 0.893}},
                {1, {0.917, 0.855, 0.818, 0.795, 0.778, 0.781, 0.798, 0.817, 0.854, 0.899}},
                {2, {0.920, 0.880, 0.850, 0.830, 0.820, 0.825, 0.840, 0.860, 0.880, 0.910}}, // ESTIMATED: smoothed, removed 0/1 boundaries
                {3, {0.850, 0.820, 0.800, 0.790, 0.785, 0.790, 0.800, 0.815, 0.835, 0.860}}  // ESTIMATED: physical lambda range [0.7-0.95]
            }
        }
    };

    // nData_ : charge = 1 (negative)
    nData_ = {
        {
            "out", {
                {0, {6.49, 6.12, 6.03, 5.94, 5.86, 5.89, 5.88, 6.02, 6.21, 6.61}},  // centrality 0: reference
                {1, {6.07, 5.51, 5.20, 5.06, 5.01, 5.01, 5.07, 5.25, 5.51, 6.05}},  // centrality 1: ~10% smaller
                {2, {5.10, 4.85, 4.65, 4.50, 4.40, 4.45, 4.55, 4.75, 5.10, 5.55}},  // ESTIMATED: smoothed trend
                {3, {4.30, 4.15, 4.00, 3.90, 3.85, 3.90, 4.00, 4.20, 4.45, 4.75}}   // ESTIMATED: peripheral, no outliers
            }
        },
        {
            "side", {
                {0, {4.36, 4.44, 4.47, 4.48, 4.52, 4.53, 4.52, 4.43, 4.39, 4.29}},  // centrality 0: reference
                {1, {3.89, 3.93, 3.86, 3.84, 3.88, 3.85, 3.88, 3.89, 3.90, 3.90}},  // centrality 1: ~10% smaller
                {2, {3.55, 3.45, 3.40, 3.35, 3.33, 3.38, 3.42, 3.50, 3.60, 3.70}},  // ESTIMATED: removed outlier 5.337
                {3, {3.15, 3.10, 3.05, 3.00, 3.02, 3.08, 3.15, 3.25, 3.35, 3.45}}   // ESTIMATED: smooth peripheral
            }
        },
        {
            "long", {
                {0, {5.33, 4.81, 4.35, 4.09, 3.96, 3.95, 4.09, 4.33, 4.81, 5.40}},  // centrality 0: reference
                {1, {4.94, 4.45, 4.11, 3.88, 3.81, 3.85, 3.93, 4.11, 4.49, 4.96}},  // centrality 1: ~10% smaller
                {2, {4.40, 4.05, 3.75, 3.55, 3.45, 3.48, 3.60, 3.80, 4.10, 4.55}},  // ESTIMATED: smoothed
                {3, {3.75, 3.55, 3.40, 3.30, 3.25, 3.30, 3.40, 3.60, 3.85, 4.15}}   // ESTIMATED: peripheral trend
            }
        },
        {
            "lambda", {
                {0, {0.834, 0.821, 0.805, 0.794, 0.793, 0.793, 0.796, 0.799, 0.812, 0.842}},
                {1, {0.851, 0.823, 0.783, 0.767, 0.769, 0.772, 0.778, 0.791, 0.818, 0.856}},
                {2, {0.900, 0.860, 0.830, 0.810, 0.800, 0.805, 0.820, 0.840, 0.865, 0.890}}, // ESTIMATED: physical range
                {3, {0.830, 0.800, 0.780, 0.770, 0.765, 0.770, 0.785, 0.800, 0.820, 0.845}}  // ESTIMATED: no 0.000/1.000 artifacts
            }
        }
    };
}

double InitialParameters::get(int charge, std::string parameter, int centrality, int y) {
    std::map<std::string, std::map<int, std::vector<double>>> data;
    if (charge == 0) {
        data = pData_;
    } else {
        data = nData_;
    }
    auto centralityMap = data.at(parameter).find(centrality);
    const auto& vec = centralityMap->second;
    
    return vec[y];
}