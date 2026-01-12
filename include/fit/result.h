#pragma once

#include <nlohmann/json.hpp>
#include <string>

#include "fit/types.h"
#include "helpers.h"

nlohmann::json BuildMeta(const std::string& inputFile);

nlohmann::json FitResultToJson(const FitResult&);

nlohmann::json BuildDataPoint(int ch,int c,int k,int y,const FitResult&);

void WriteFitJson(const std::string& inputFile,
                  const std::string& outDir,
                  FitResult res[chargeSize][centralitySize][ktSize][rapiditySize]);
