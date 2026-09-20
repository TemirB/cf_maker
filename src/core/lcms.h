#pragma once

#include <string>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>

enum class LCMSAxis : char
{
    Out,
    Side,
    Long
};

[[nodiscard]] std::string axis_name(LCMSAxis a);

[[nodiscard]] char projection_char(LCMSAxis a);

[[nodiscard]] LCMSAxis third_axis(LCMSAxis a, LCMSAxis b);

void reset_ranges(TH3D& h);

void slice_projection(TH3D& h, LCMSAxis axis, double w = 0.05);

void slice_freeze(TH3D& h, LCMSAxis freeze, double w = 0.2);

[[nodiscard]] TH1D* project_1d(TH3D& h, LCMSAxis axis, double w = 0.05);

[[nodiscard]] TH2D* project_2d(TH3D& h, LCMSAxis ax1, LCMSAxis ax2, double w = 0.2);
