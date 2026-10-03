#include <cmath>
#include <iostream>
#include <stdexcept>

#include <TH1.h>
#include <TH3.h>
#include <TMemFile.h>
#include <TROOT.h>

#include "analysis/ratio.h"
#include "io/input.h"

namespace
{
void require_close(double actual, double expected)
{
    if (std::abs(actual - expected) > 1e-12) {
        throw std::runtime_error("Unexpected normalized charge-ratio projection");
    }
}

void check_ratio(double slice, bool sparse)
{
    Config cfg;
    cfg.input.type = "kt";
    cfg.selection.centralities = {0};
    cfg.binning.count = 1;
    cfg.binning.names = {"test"};
    cfg.projections.slice_ratio = slice;
    TMemFile input("charge_cf.root", "RECREATE");
    for (int charge = 0; charge < 2; ++charge) {
        TH3D histogram("cf", "", 8, -.2, .2, 8, -.2, .2, 8, -.2, .2);
        for (int x = 1; x <= 8; ++x) {
            for (int y = 1; y <= 8; ++y) {
                for (int z = 1; z <= 8; ++z) {
                    if (!sparse || (y == 4 && z == 4)) {
                        histogram.SetBinContent(x, y, z, 1);
                        histogram.SetBinError(x, y, z, charge == 0 ? .1 : .2);
                    }
                }
            }
        }
        histogram.Write(get_cf_name(charge, 0, "kt", "test").c_str());
    }
    TMemFile ratio_project("ratio_project.root", "RECREATE");
    TMemFile project_ratio("project_ratio.root", "RECREATE");
    do_cf_ratios(cfg, &input, &ratio_project, &project_ratio);
    auto* result = dynamic_cast<TH1*>(project_ratio.Get("proj_of_ratios_0_0_out"));
    if (!result) {
        throw std::runtime_error("Missing charge-ratio projection");
    }
    require_close(result->GetBinContent(4), 1);
    // For width .05 the slice has four cells; empty cells carry no weight.
    if (slice == .05) {
        require_close(result->GetBinError(4), std::sqrt(.05 / (sparse ? 1 : 4)));
    }
}
} // namespace

int main()
{
    gROOT->SetBatch(true);
    TH1::AddDirectory(false);
    check_ratio(.05, false);
    check_ratio(.15, false);
    check_ratio(.05, true);
    std::cout << "Charge-ratio normalization tests passed\n";
}
