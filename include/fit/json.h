#pragma once

#include <string>
#include <nlohmann/json.hpp>

#include "fit/types.h"
#include "helpers.h"

nlohmann::json BuildMeta(const std::string& inputFile);
nlohmann::json BuildDataPoint(int ch,int c,int k,int y,const FitResult&);

void WriteFitJson(
    const std::string& inputFile,
    const std::string& outDir,
    FitResult fitRes[chargeSize][centralitySize][ktSize]
);

FitResult FitResultFromJson(const nlohmann::json& j);