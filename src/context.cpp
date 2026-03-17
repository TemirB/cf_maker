#include "context.h"

#include <limits.h>
#include <libgen.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "helpers.h"

std::string GetExeDir() {
    char buf[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf)-1);
    buf[len] = '\0';
    return std::string(dirname(buf));
}

void EnsureDir(const std::string& dir) {
    struct stat st;
    if (stat(dir.c_str(), &st) != 0) {
        mkdir(dir.c_str(), 0755);
    }
}

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