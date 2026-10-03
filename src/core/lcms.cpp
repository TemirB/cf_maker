#include "core/lcms.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <TAxis.h>

namespace
{
void set_slice_range(TAxis& axis, double half_width)
{
    if (!std::isfinite(half_width) || half_width <= 0) {
        throw std::invalid_argument("projection slice width must be finite and positive");
    }

    std::vector<double> edges;
    edges.reserve(static_cast<std::size_t>(axis.GetNbins()) + 1);
    for (int edge = 1; edge <= axis.GetNbins() + 1; ++edge) {
        edges.push_back(axis.GetBinLowEdge(edge));
    }
    const auto snap_to_edge = [&](double boundary) {
        const auto next = std::lower_bound(edges.begin(), edges.end(), boundary);
        double snapped = boundary;
        double distance = std::numeric_limits<double>::infinity();
        const double scale =
            std::max({std::abs(axis.GetXmin()), std::abs(axis.GetXmax()), std::abs(boundary)});
        // Nominal decimal cuts can differ from stored edges by a few ULPs.
        // Limit snapping to a quarter of the interval so narrow cuts stay open.
        const double tolerance =
            std::min(8 * std::numeric_limits<double>::epsilon() * scale, half_width / 2);
        const auto consider = [&](double edge) {
            const double difference = std::abs(edge - boundary);
            if (difference <= tolerance && difference < distance) {
                snapped = edge;
                distance = difference;
            }
        };
        if (next != edges.end()) {
            consider(*next);
        }
        if (next != edges.begin()) {
            consider(*(next - 1));
        }
        return snapped;
    };

    const double lower = snap_to_edge(-half_width);
    const double upper = snap_to_edge(half_width);
    // Include complete bins with positive overlap, excluding bins which merely
    // touch a cut. Searching actual edges avoids ROOT's uniform-bin arithmetic,
    // whose rounding differs between architectures and compiler settings.
    const int first =
        static_cast<int>(std::upper_bound(edges.begin(), edges.end(), lower) - edges.begin());
    const int last =
        static_cast<int>(std::lower_bound(edges.begin(), edges.end(), upper) - edges.begin());
    // Project3D treats an underflow-only (0, 0) range as all regular bins,
    // even when kAxisRange is set. Reject it instead of projecting other data.
    if (first == 0 && last == 0) {
        throw std::invalid_argument(
            "projection slice selects only underflow; ROOT Project3D cannot represent this range");
    }
    axis.SetRange(first, last);
}
} // namespace

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
        set_slice_range(*h.GetXaxis(), w);
    }
    if (axis != LCMSAxis::Side) {
        set_slice_range(*h.GetYaxis(), w);
    }
    if (axis != LCMSAxis::Long) {
        set_slice_range(*h.GetZaxis(), w);
    }
}

void slice_freeze(TH3D& h, LCMSAxis freeze, double w)
{
    reset_ranges(h);
    if (freeze == LCMSAxis::Out) {
        set_slice_range(*h.GetXaxis(), w);
    }
    if (freeze == LCMSAxis::Side) {
        set_slice_range(*h.GetYaxis(), w);
    }
    if (freeze == LCMSAxis::Long) {
        set_slice_range(*h.GetZaxis(), w);
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
    // ROOT lists the vertical axis first and the horizontal axis second.
    proj += projection_char(ax2);
    proj += projection_char(ax1);
    proj += 'e';
    TH2D* out = static_cast<TH2D*>(h.Project3D(proj.c_str()));
    out->SetDirectory(nullptr);
    return out;
}
