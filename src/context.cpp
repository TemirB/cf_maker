#include "context.h"

#include <iostream>
#include <limits.h>
#include <libgen.h>


#include "helpers.h"

std::string GetExeDir() {
    char buf[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf)-1);
    buf[len] = '\0';
    return std::string(dirname(buf));
}

Context BuildContext(char** argv) {
    Context ctx;

    std::string resDir  = GetExeDir() + "/../results";

    EnsureDir(resDir);

    ctx.inputFile = argv[1];
    ctx.outDir    = resDir + "/" + argv[2];
    EnsureDir(ctx.outDir);


    ctx.bining.type = argv[3];
    if (std::strcmp(ctx.bining.type, "kt") == 0) {
        ctx.bining.count = Kt::kCount;
        ctx.bining.names = Kt::kNames;
        ctx.bining.values = Kt::kValues;
    } else if (std::strcmp(ctx.bining.type, "rapidity") == 0) {
        ctx.bining.count = Rapidity::kCount;
        ctx.bining.names = Rapidity::kNames;
        ctx.bining.values = Rapidity::kValues;
    }


    ctx.fitRes = FitGrid(
        Charge::kCount,
        std::vector<std::vector<FitResult>>(
            Centrality::kCount,
            std::vector<FitResult>(ctx.bining.count)
        )
    );
    return ctx;
}