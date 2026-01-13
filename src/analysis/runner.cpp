#include "analysis/runner.h"

#include <TROOT.h>

void RunAllStages(RunContext& ctx) {
    for (auto& stage : GetRegistry()) {
        gROOT->cd();              // ← КРИТИЧНО
        stage.fn(ctx);

        for (auto& [_, f] : ctx.outputs) {
            f->Write();
            f->Close();
            delete f;
        }
        ctx.outputs.clear();
    }
}
