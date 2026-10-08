#include "analysis/pipeline.h"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>

#include <TFile.h>
#include <TString.h>

#include "analysis/cf3d.h"
#include "analysis/dependency.h"
#include "analysis/projections1d.h"
#include "analysis/projections2d.h"
#include "analysis/ratio.h"
#include "core/log.h"
#include "io/fit_results.h"
#include "io/output.h"

void prepare_stage_dependencies(Config& cfg)
{
    if (cfg.stages.ratios && (std::find(cfg.selection.charges.begin(), cfg.selection.charges.end(),
                                        0) == cfg.selection.charges.end() ||
                              std::find(cfg.selection.charges.begin(), cfg.selection.charges.end(),
                                        1) == cfg.selection.charges.end())) {
        throw std::runtime_error("ratios require both selected charges");
    }
    if (cfg.stages.cf3d ||
        !(cfg.stages.dependency || cfg.stages.projections_1d || cfg.stages.ratios)) {
        return;
    }
    const std::string name = cfg.output.dir + "/cf3d.root";
    TFile previous(name.c_str(), "READ");
    if (previous.IsZombie()) {
        throw std::runtime_error("cannot resume: missing or unreadable " + name +
                                 "; rerun with stages.cf3d=true");
    }
    read_fit_results(previous, cfg);
    logging::info("Restored complete fit results from " + name);
}

void stage_cf3d(Config& cfg)
{
    logging::ScopedTimer timer("Stage cf3d");
    const std::string name = cfg.output.dir + "/cf3d.root";
    logging::info("Stage cf3d: output = " + name);
    auto f = create_output_file(name);

    build_and_fit_3d_correlation_functions(cfg, f.get());
    write_fit_results(*f, cfg);
    finish_output_file(*f);
}

void stage_dependency(Config& cfg)
{
    logging::ScopedTimer timer("Stage dependency");
    const std::string name = Form("%s/%s.root", cfg.output.dir.c_str(), cfg.input.type.c_str());
    logging::info("Stage dependency: output = " + name);

    const std::string cf3dName = cfg.output.dir + "/cf3d.root";
    TFile cf3d(cf3dName.c_str(), "READ");
    if (cf3d.IsZombie()) {
        throw std::runtime_error("cannot open saved 3D CF: " + cf3dName);
    }
    auto f = create_output_file(name);

    make_dependency(cfg, &cf3d, f.get());
    finish_output_file(*f);
}

void stage_1d_projections(Config& cfg, TFile* input)
{
    logging::ScopedTimer timer("Stage projections_1d");
    const std::string name = cfg.output.dir + "/1d.root";
    logging::info("Stage projections_1d: output = " + name);
    auto f = create_output_file(name);

    make_lcms_1d_projections(cfg, input, f.get());

    finish_output_file(*f);
}

void stage_2d_projections(Config& cfg, TFile* input)
{
    logging::ScopedTimer timer("Stage projections_2d");
    const std::string name = cfg.output.dir + "/2d.root";
    logging::info("Stage projections_2d: output = " + name);
    auto f = create_output_file(name);

    make_lcms_2d_projections(cfg, input, f.get());

    finish_output_file(*f);
}

void stage_ratios(Config& cfg)
{
    logging::ScopedTimer timer("Stage ratios");
    const std::string name1 = cfg.output.dir + "/ratio_projs.root";
    const std::string name2 = cfg.output.dir + "/proj_ratios.root";
    const std::string cf3d = cfg.output.dir + "/cf3d.root";
    logging::info("Stage ratios: outputs = " + name1 + ", " + name2);

    auto fCF3D = std::make_unique<TFile>(cf3d.c_str(), "READ");
    if (fCF3D->IsZombie()) {
        throw std::runtime_error("cannot open saved 3D CF: " + cf3d);
    }
    auto f1 = create_output_file(name1);
    auto f2 = create_output_file(name2);
    do_cf_ratios(cfg, fCF3D.get(), f1.get(), f2.get());

    finish_output_file(*f1);
    finish_output_file(*f2);
}
