#pragma once

#include <string>

struct Context {
    std::string exeDir;
    std::string projectRoot;
    std::string resDir;
    std::string inputFile;
    std::string outDir;
    std::string cf3dFile;
    std::string fitJsonFile;
};

Context BuildContext(char** argv);