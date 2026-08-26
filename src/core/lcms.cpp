#include "core/lcms.h"

std::string axis_name(LCMSAxis a)
{
    switch (a) {
    case LCMSAxis::Out:
        return "out";
    case LCMSAxis::Side:
        return "side";
    case LCMSAxis::Long:
        return "long";
    }
    return "";
}

char projection_char(LCMSAxis a)
{
    switch (a) {
    case LCMSAxis::Out:
        return 'x';
    case LCMSAxis::Side:
        return 'y';
    case LCMSAxis::Long:
        return 'z';
    }
    return '\0';
}

LCMSAxis third_axis(LCMSAxis a, LCMSAxis b)
{
    if (a != LCMSAxis::Out && b != LCMSAxis::Out) {
        return LCMSAxis::Out;
    }
    if (a != LCMSAxis::Side && b != LCMSAxis::Side) {
        return LCMSAxis::Side;
    }
    return LCMSAxis::Long;
}

void reset_ranges(TH3D& h)
{
    h.GetXaxis()->SetRange(0, 0);
    h.GetYaxis()->SetRange(0, 0);
    h.GetZaxis()->SetRange(0, 0);
}

void slice_projection(TH3D& h, LCMSAxis axis, double w)
{
    reset_ranges(h);
    if (axis != LCMSAxis::Out) {
        h.GetXaxis()->SetRangeUser(-w, w);
    }
    if (axis != LCMSAxis::Side) {
        h.GetYaxis()->SetRangeUser(-w, w);
    }
    if (axis != LCMSAxis::Long) {
        h.GetZaxis()->SetRangeUser(-w, w);
    }
}

void slice_freeze(TH3D& h, LCMSAxis freeze, double w)
{
    reset_ranges(h);
    if (freeze == LCMSAxis::Out) {
        h.GetXaxis()->SetRangeUser(-w, w);
    }
    if (freeze == LCMSAxis::Side) {
        h.GetYaxis()->SetRangeUser(-w, w);
    }
    if (freeze == LCMSAxis::Long) {
        h.GetZaxis()->SetRangeUser(-w, w);
    }
}

TH1D* project_1d(TH3D& h, LCMSAxis axis, double w)
{
    slice_projection(h, axis, w);
    char proj[3] = {projection_char(axis), 'e', '\0'};
    TH1D* out = static_cast<TH1D*>(h.Project3D(proj));
    out->SetDirectory(nullptr);
    return out;
}

TH2D* project_2d(TH3D& h, LCMSAxis ax1, LCMSAxis ax2, double w)
{
    slice_freeze(h, third_axis(ax1, ax2), w);
    std::string proj;
    proj += projection_char(ax1);
    proj += projection_char(ax2);
    proj += 'e';
    TH2D* out = static_cast<TH2D*>(h.Project3D(proj.c_str()));
    out->SetDirectory(nullptr);
    return out;
}
