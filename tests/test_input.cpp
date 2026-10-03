#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include <TH3F.h>
#include <TMemFile.h>
#include <TNamed.h>
#include <TROOT.h>

#include "io/input.h"

namespace
{
void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void compare_axis(const TAxis& source, const TAxis& copy)
{
    require(source.GetNbins() == copy.GetNbins(), "Axis bin count changed");
    require(source.GetXbins()->GetSize() == copy.GetXbins()->GetSize(),
            "Uniform/variable binning changed");
    require(std::string(source.GetTitle()) == copy.GetTitle(), "Axis title changed");
    require(source.GetFirst() == copy.GetFirst() && source.GetLast() == copy.GetLast(),
            "Axis range changed");
    for (int bin = 0; bin <= source.GetNbins() + 1; ++bin) {
        require(source.GetBinLowEdge(bin) == copy.GetBinLowEdge(bin), "Axis edges changed");
        require(std::string(source.GetBinLabel(bin)) == copy.GetBinLabel(bin),
                "Axis labels changed");
    }
}

void compare_histogram(const TH3& source, const TH3D& copy)
{
    require(copy.GetDirectory() == nullptr, "Input copy is still owned by a directory");
    require(copy.IsA() == TH3D::Class(), "Input was not converted to TH3D");
    require(std::string(source.GetName()) == copy.GetName(), "Histogram name changed");
    require(std::string(source.GetTitle()) == copy.GetTitle(), "Histogram title changed");
    compare_axis(*source.GetXaxis(), *copy.GetXaxis());
    compare_axis(*source.GetYaxis(), *copy.GetYaxis());
    compare_axis(*source.GetZaxis(), *copy.GetZaxis());
    require(source.GetNcells() == copy.GetNcells(), "Cell count changed");
    require(source.GetSumw2N() == copy.GetSumw2N(), "Sumw2 storage changed");
    for (int bin = 0; bin < source.GetNcells(); ++bin) {
        require(source.GetBinContent(bin) == copy.GetBinContent(bin),
                "Contents changed, including a flow bin");
        require(source.GetBinError(bin) == copy.GetBinError(bin), "Bin errors changed");
        if (source.GetSumw2N() != 0) {
            require(source.GetSumw2()->At(bin) == copy.GetSumw2()->At(bin),
                    "Squared weights changed");
        }
    }
    require(source.GetEntries() == copy.GetEntries(), "Entry count changed");
    require(source.GetStatOverflows() == copy.GetStatOverflows(), "Overflow statistics changed");
    double source_stats[TH1::kNstat]{};
    double copy_stats[TH1::kNstat]{};
    source.GetStats(source_stats);
    copy.GetStats(copy_stats);
    for (int index = 0; index < TH1::kNstat; ++index) {
        require(source_stats[index] == copy_stats[index], "Cached statistics changed");
    }
}

void populate(TH3& histogram, double scale)
{
    histogram.Sumw2();
    histogram.SetStatOverflows(TH1::kConsider);
    for (int x = 0; x <= histogram.GetNbinsX() + 1; ++x) {
        for (int y = 0; y <= histogram.GetNbinsY() + 1; ++y) {
            for (int z = 0; z <= histogram.GetNbinsZ() + 1; ++z) {
                const double weight = scale * (1.25 + x + 2.0 * y + 3.0 * z);
                // Fill off-centre to exercise cached moments, not moments reconstructed
                // from bin centres during conversion.
                histogram.Fill(histogram.GetXaxis()->GetBinCenter(x) + 0.001,
                               histogram.GetYaxis()->GetBinCenter(y) + 0.001,
                               histogram.GetZaxis()->GetBinCenter(z) + 0.001, weight);
            }
        }
    }
    histogram.GetXaxis()->SetTitle("q_out [GeV/c]");
    histogram.GetYaxis()->SetTitle("q_side [GeV/c]");
    histogram.GetZaxis()->SetTitle("q_long [GeV/c]");
    histogram.GetXaxis()->SetBinLabel(1, "out bin");
    histogram.SetEntries(123.5);
}

void check_float_input()
{
    TMemFile input("float_input.root", "RECREATE");
    const double x_edges[] = {-0.3, -0.1, 0.2};
    const double y_edges[] = {-0.4, -0.2, 0.0, 0.4};
    const double z_edges[] = {-0.5, 0.1, 0.6};
    TH3F den("bp_0_0_num_0", "Pair counts", 2, x_edges, 3, y_edges, 2, z_edges);
    TH3F num("bp_0_0_num_wei_0", "Pair weights", 2, x_edges, 3, y_edges, 2, z_edges);
    den.SetDirectory(nullptr);
    num.SetDirectory(nullptr);
    populate(den, 1.0);
    populate(num, 1.5);
    den.Write();
    num.Write();

    const auto [den_ptr, num_ptr] = get_hists(&input, 0, 0, 0);
    std::unique_ptr<TH3D> den_copy(den_ptr), num_copy(num_ptr);
    require(den_copy && num_copy, "TH3F input was rejected");
    compare_histogram(den, *den_copy);
    compare_histogram(num, *num_copy);

    auto* stored_den = dynamic_cast<TH3F*>(input.Get(den.GetName()));
    require(stored_den != nullptr, "Source histogram disappeared after conversion");
    std::unique_ptr<TH3F> stored_den_owner;
    if (stored_den->GetDirectory() == nullptr) {
        stored_den_owner.reset(stored_den);
    }
    const auto [second_den_ptr, second_num_ptr] = get_hists(&input, 0, 0, 0);
    std::unique_ptr<TH3D> second_den(second_den_ptr), second_num(second_num_ptr);
    require(second_den && second_num, "Repeated input read failed");
    require(second_den.get() != den_copy.get() && second_num.get() != num_copy.get(),
            "Repeated reads returned shared ownership");
    den_copy->SetBinContent(1, 1, 1, -42.0);
    require(second_den->GetBinContent(1, 1, 1) == den.GetBinContent(1, 1, 1),
            "Mutating one copy changed another copy");
    require(stored_den->GetBinContent(1, 1, 1) == den.GetBinContent(1, 1, 1),
            "Mutating a copy changed the input file's object");
    input.Close();
    compare_histogram(num, *num_copy);
    compare_histogram(den, *second_den);
}

void check_double_input(bool sumw2)
{
    TMemFile input(sumw2 ? "double_input.root" : "counts_input.root", "RECREATE");
    TH3D den("bp_0_0_num_0", "", 2, -0.2, 0.2, 3, -0.3, 0.3, 4, -0.4, 0.4);
    TH3D num("bp_0_0_num_wei_0", "", 2, -0.2, 0.2, 3, -0.3, 0.3, 4, -0.4, 0.4);
    den.SetDirectory(nullptr);
    num.SetDirectory(nullptr);
    if (sumw2) {
        populate(den, 1.0);
        populate(num, 1.5);
    } else {
        den.Fill(0.01, 0.01, 0.01);
        num.Fill(0.01, 0.01, 0.01);
        den.Sumw2(false);
        num.Sumw2(false);
    }
    den.GetYaxis()->SetRange(2, 3);
    num.GetZaxis()->SetRange(1, 2);
    den.Write();
    num.Write();
    const auto [den_ptr, num_ptr] = get_hists(&input, 0, 0, 0);
    std::unique_ptr<TH3D> den_copy(den_ptr), num_copy(num_ptr);
    require(den_copy && num_copy, "TH3D input was rejected");
    compare_histogram(den, *den_copy);
    compare_histogram(num, *num_copy);
    den.GetYaxis()->SetRange(0, 0);
    num.GetZaxis()->SetRange(0, 0);
    den_copy->GetYaxis()->SetRange(0, 0);
    num_copy->GetZaxis()->SetRange(0, 0);
    compare_histogram(den, *den_copy);
    compare_histogram(num, *num_copy);
}

void check_invalid_input()
{
    require(get_hists(nullptr, 0, 0, 0) == std::pair<TH3D*, TH3D*>{nullptr, nullptr},
            "Null input did not fail cleanly");
    TMemFile input("invalid_input.root", "RECREATE");
    require(get_hists(&input, 0, 0, 0) == std::pair<TH3D*, TH3D*>{nullptr, nullptr},
            "Missing denominator did not fail cleanly");
    TH3F den("bp_0_0_num_0", "", 1, 0, 1, 1, 0, 1, 1, 0, 1);
    den.SetDirectory(nullptr);
    den.Write();
    require(get_hists(&input, 0, 0, 0) == std::pair<TH3D*, TH3D*>{nullptr, nullptr},
            "Missing numerator did not fail cleanly");
    TNamed unsupported("bp_0_0_num_wei_0", "not a histogram");
    unsupported.Write();
    require(get_hists(&input, 0, 0, 0) == std::pair<TH3D*, TH3D*>{nullptr, nullptr},
            "Unsupported object class was accepted");
}
} // namespace

int main()
{
    gROOT->SetBatch(true);
    // Exercise both ROOT ownership policies; get_hists must be independent of them.
    for (const bool register_histograms : {false, true}) {
        TH1::AddDirectory(register_histograms);
        check_float_input();
        check_double_input(true);
        check_double_input(false);
        check_invalid_input();
    }
    TH1::AddDirectory(false);
    std::cout << "All input tests passed\n";
}
