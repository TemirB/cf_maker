#include "context.h"

#include "helpers.h"

Context BuildContext(char** argv) {
    Context ctx;

    ctx.exeDir      = GetExeDir();
    ctx.projectRoot = ctx.exeDir + "/..";
    ctx.resDir  = ctx.projectRoot + "/results";

    EnsureDir(ctx.resDir);

    ctx.inputFile = argv[1];
    ctx.outDir    = ctx.resDir + "/" + argv[2];
    EnsureDir(ctx.outDir);

    ctx.cf3dFile  = ctx.outDir + "/cf3d.root";
    ctx.fitJsonFile   = ctx.outDir + "/fit3d.json";

    return ctx;
}