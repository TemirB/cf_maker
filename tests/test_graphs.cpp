#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH3D.h>
#include <TList.h>
#include <TMemFile.h>
#include <TMultiGraph.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TText.h>

#include "analysis/dependency.h"
#include "analysis/graphs.h"
#include "io/input.h"

namespace
{
void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

Config make_config()
{
    Config cfg;
    cfg.input.type = "kt";
    cfg.output.dir = "graphs_test_output";
    cfg.general.images.need = false;
    cfg.selection.charges = {0};
    cfg.selection.centralities = {0};
    cfg.binning.count = 7;
    cfg.binning.values = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    cfg.binning.names = {"0", "1", "2", "3", "4", "5", "6"};
    cfg.fit_results = FitGrid(1, std::vector<std::vector<FitResult>>(1, std::vector<FitResult>(7)));

    FitResult successful;
    successful.ok = true;
    successful.status = 0;
    successful.cov_status = 3;
    successful.attempts = 1;
    // These values deliberately fall outside the old is_valid() heuristic.
    successful.r = {12.0, 0.5, 6.0, 0.0, 0.0, 0.0};
    successful.e_r.fill(0.1);
    successful.lambda = 1.0;
    successful.e_lambda = 0.01;
    successful.chi2 = 4.0;
    successful.ndf = 2;
    successful.p_value = 0.2;
    auto& fits = cfg.fit_results[0][0];
    fits[1] = successful;
    fits[2] = successful;
    fits[2].ok = false;
    fits[2].status = 4;
    fits[3] = successful;
    fits[3].e_r[4] = std::numeric_limits<double>::infinity();
    fits[4] = successful;
    fits[4].at_limit = true;
    fits[5] = successful;
    fits[5].ndf = 0;
    fits[6] = successful;
    fits[6].r[0] = 13.0;
    fits[6].chi2 = 8.0;
    fits[6].p_value = 0.4;
    return cfg;
}

void check_two_points(const TGraphErrors& graph, double first_y, double second_y)
{
    require(graph.GetN() == 2, "Unavailable fits were retained or valid fits were rejected");
    double x = 0.0;
    double y = 0.0;
    graph.GetPoint(0, x, y);
    require(x == 1.5 && std::abs(y - first_y) < 1e-12, "Wrong first point or placeholder point");
    graph.GetPoint(1, x, y);
    require(x == 6.5 && std::abs(y - second_y) < 1e-12, "Wrong second point or placeholder point");
}

TMultiGraph& read_multigraph(TMemFile& file, const char* name)
{
    auto* canvas = dynamic_cast<TCanvas*>(file.Get(name));
    require(canvas != nullptr, std::string("Missing dependency canvas: ") + name);
    for (TObject* object : *canvas->GetListOfPrimitives()) {
        if (auto* graph = dynamic_cast<TMultiGraph*>(object)) {
            return *graph;
        }
    }
    throw std::runtime_error("Dependency canvas does not contain a multigraph");
}

void check_empty_canvas(TMemFile& file, const char* name)
{
    auto* canvas = dynamic_cast<TCanvas*>(file.Get(name));
    require(canvas != nullptr, std::string("Missing empty dependency canvas: ") + name);
    bool message_found = false;
    for (TObject* object : *canvas->GetListOfPrimitives()) {
        if (auto* graph = dynamic_cast<TMultiGraph*>(object)) {
            for (TObject* component : *graph->GetListOfGraphs()) {
                const auto* points = dynamic_cast<TGraph*>(component);
                require(!points || points->GetN() == 0, "Empty canvas contains physical points");
            }
        }
        const auto* message = dynamic_cast<TPaveText*>(object);
        if (!message) {
            continue;
        }
        for (TObject* line : *message->GetListOfLines()) {
            const auto* text = dynamic_cast<TText*>(line);
            message_found |= text && std::string(text->GetTitle()) == "No usable fits";
        }
    }
    require(message_found, "Empty graph has no unavailable-fit annotation");
}
} // namespace

int main()
{
    gROOT->SetBatch(true);
    TH1::AddDirectory(false);
    Config cfg = make_config();
    for (int bin = 0; bin < cfg.binning.count; ++bin) {
        require(is_usable_fit(cfg.fit_results[0][0][bin]) == (bin == 1 || bin == 6),
                "Unexpected fit eligibility policy");
    }
    TMemFile cf_file("graphs_cf.root", "RECREATE");
    for (int bin = 0; bin < cfg.binning.count; ++bin) {
        TH3D cf(get_cf_name(0, 0, "kt", cfg.binning.names[bin]).c_str(), "", 1, -0.01, 0.01, 1,
                -0.01, 0.01, 1, -0.01, 0.01);
        cf.SetBinContent(1, 1, 1, 2.0);
        cf.Write();
    }
    {
        std::unique_ptr<TGraphErrors> chi2(build_chi2_ndf_graph(cfg, 0, 0));
        std::unique_ptr<TGraphErrors> pvalue(build_pvalue_graph(cfg, 0, 0));
        std::unique_ptr<TGraphErrors> ratio(build_fit_over_cf_graph(cfg, &cf_file, 0, 0));
        check_two_points(*chi2, 2.0, 4.0);
        check_two_points(*pvalue, 0.2, 0.4);
        check_two_points(*ratio, 1.0, 1.0);
    }
    {
        TMemFile output("graphs_output.root", "RECREATE");
        make_dependency(cfg, &cf_file, &output);
        auto& radii = read_multigraph(output, "mg_R_out_pos");
        auto* graph = dynamic_cast<TGraphErrors*>(radii.GetListOfGraphs()->At(0));
        require(graph != nullptr, "Missing radius graph");
        check_two_points(*graph, 12.0, 13.0);
        auto& lambda = read_multigraph(output, "mg_L_pos");
        graph = dynamic_cast<TGraphErrors*>(lambda.GetListOfGraphs()->At(0));
        require(graph != nullptr, "Missing lambda graph");
        check_two_points(*graph, 1.0, 1.0);
        require(cfg.fit_results[0][0][4].at_limit, "Boundary fit metadata was discarded");
        require(!cfg.fit_results[0][0][2].ok, "Failure metadata was discarded");
    }
    for (FitResult& result : cfg.fit_results[0][0]) {
        result.ok = false;
    }
    {
        std::unique_ptr<TGraphErrors> chi2(build_chi2_ndf_graph(cfg, 0, 0));
        std::unique_ptr<TGraphErrors> pvalue(build_pvalue_graph(cfg, 0, 0));
        std::unique_ptr<TGraphErrors> ratio(build_fit_over_cf_graph(cfg, &cf_file, 0, 0));
        require(chi2->GetN() == 0 && pvalue->GetN() == 0 && ratio->GetN() == 0,
                "Entirely unavailable data produced graph points");
        TMemFile output("graphs_empty_output.root", "RECREATE");
        make_dependency(cfg, &cf_file, &output);
        check_empty_canvas(output, "mg_R_out_pos");
        check_empty_canvas(output, "mg_L_pos");
        check_empty_canvas(output, "mg_chi2_ndf_pos");
        check_empty_canvas(output, "mg_FitOverCF_pos");
    }
    std::filesystem::remove_all(cfg.output.dir);
    std::cout << "All graph tests passed\n";
}
