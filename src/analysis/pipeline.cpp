#include "analysis/pipeline.h"

#include <memory>
#include <string>

#include <TFile.h>
#include <TString.h>

#include "analysis/cf3d.h"
#include "analysis/dependency.h"
#include "analysis/projections1d.h"
#include "analysis/projections2d.h"
#include "analysis/ratio.h"
#include "core/log.h"

void stage_cf3d(Config& cfg)
{
    log::ScopedTimer timer("Stage cf3d");
    const std::string name = cfg.output.dir + "/cf3d.root";
    log::Info("Stage cf3d: output = " + name);
    TFile f(name.c_str(), "RECREATE");

    build_and_fit_3d_correlation_functions(cfg, &f);

    f.Write();
}

void stage_dependency(Config& cfg)
{
    log::ScopedTimer timer("Stage dependency");
    const std::string name = Form("%s/%s.root", cfg.output.dir.c_str(), cfg.input.type.c_str());
    log::Info("Stage dependency: output = " + name);

    const std::string cf3dName = cfg.output.dir + "/cf3d.root";
    TFile f(name.c_str(), "RECREATE");
    TFile cf3d(cf3dName.c_str(), "READ");

    make_dependency(cfg, &cf3d, &f);
}

void stage_1d_projections(Config& cfg, TFile* input)
{
    log::ScopedTimer timer("Stage projections_1d");
    const std::string name = cfg.output.dir + "/1d.root";
    log::Info("Stage projections_1d: output = " + name);
    auto f = std::make_unique<TFile>(name.c_str(), "RECREATE");

    make_lcms_1d_projections(cfg, input, f.get());

    f->Write();
}

void stage_2d_projections(Config& cfg, TFile* input)
{
    log::ScopedTimer timer("Stage projections_2d");
    const std::string name = cfg.output.dir + "/2d.root";
    log::Info("Stage projections_2d: output = " + name);
    auto f = std::make_unique<TFile>(name.c_str(), "RECREATE");

    make_lcms_2d_projections(cfg, input, f.get());

    f->Write();
}

void stage_ratios(Config& cfg)
{
    log::ScopedTimer timer("Stage ratios");
    const std::string name1 = cfg.output.dir + "/ratio_projs.root";
    const std::string name2 = cfg.output.dir + "/proj_ratios.root";
    const std::string cf3d = cfg.output.dir + "/cf3d.root";
    log::Info("Stage ratios: outputs = " + name1 + ", " + name2);

    auto f1 = std::make_unique<TFile>(name1.c_str(), "RECREATE");
    auto f2 = std::make_unique<TFile>(name2.c_str(), "RECREATE");
    auto fCF3D = std::make_unique<TFile>(cf3d.c_str(), "READ");
    do_cf_ratios(cfg, fCF3D.get(), f1.get(), f2.get());

    f1->Write();
    f2->Write();
}
