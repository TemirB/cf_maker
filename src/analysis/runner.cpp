#include "analysis/runner.h"

void RunAllStages(RunContext& ctx) {
    for (auto& stage : GetRegistry()) {
        // std::cout << "Running " << stage.name << "\n";
        stage.fn(ctx);

        // закрываем все файлы, которые открыла эта stage
        for (auto& [_, f] : ctx.outputs) {
            f->Write();
            f->Close();
        }
        ctx.outputs.clear();
    }
}
