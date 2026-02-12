#include "plots.h"

#include <TFile.h>
#include <TDirectory.h>
#include <TH3D.h>
#include <TH1D.h>
#include <TLine.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TCanvas.h>

#include "fit/types.h"
#include "helpers.h"

template<class T>
using RootPtr = std::unique_ptr<T>;

std::string GetCentralityLabel(int centIdx) {
    static const std::vector<std::string> labels = {
        "0-5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-80%"
    };
    return (centIdx >= 0 && centIdx < (int)labels.size()) 
        ? labels[centIdx] : Form("cent %d", centIdx);
}

std::string GetKTRangeLabel(int ktIdx) {
    static const std::vector<std::string> labels = {
        "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.8", "0.8-1.0"
    };
    return (ktIdx >= 0 && ktIdx < (int)labels.size()) 
        ? Form("k_{T} = %s GeV/c", labels[ktIdx].c_str()) 
        : Form("k_{T} bin %d", ktIdx);
}

std::string FormatRatioName(
    const char* prefix, 
    int centIdx, 
    int ktIdx, 
    LCMSAxis axis
) {
    return Form("%s_%d_%d_%s", prefix, centIdx, ktIdx, ToString(axis).c_str());
}

void ApplyRatioStyle(TH1D* h, LCMSAxis axis) {
    if (!h) return;
    
    // Подписи осей
    const char* axisLabel = 
        axis == LCMSAxis::Out  ? "q_{out} (GeV/c)" :
        axis == LCMSAxis::Side ? "q_{side} (GeV/c)" : "q_{long} (GeV/c)";
    
    h->GetXaxis()->SetTitle(axisLabel);
    h->GetYaxis()->SetTitle("CF^{-} / CF^{+}");
    h->SetTitle("");
    
    // Стиль осей
    h->GetXaxis()->SetTitleSize(0.055);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->CenterTitle();
    
    h->GetYaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.0);
    h->GetYaxis()->CenterTitle();
    
    // Стиль маркеров (частицы/античастицы)
    h->SetMarkerStyle(kFullCircle);   // 20 = filled circle
    h->SetMarkerSize(1.2);
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->SetLineWidth(1);
    
    // Адаптивный диапазон Y с защитой от выбросов
    double ymin = h->GetMinimum(0.0);  // игнорируем бины с 0
    double ymax = h->GetMaximum();
    
    if (ymin < 0.2) ymin = 0.2;
    if (ymax > 2.0) ymax = 2.0;
    
    double range = ymax - ymin;
    h->GetYaxis()->SetRangeUser(
        std::max(0.0, ymin - 0.15 * range),
        ymax + 0.15 * range
    );
    
    // Устанавливаем разумные границы по умолчанию если гистограмма пустая
    if (ymin == 0.0 && ymax == 0.0) {
        h->GetYaxis()->SetRangeUser(0.5, 1.5);
    }
}

// void DrawAndSaveRatio(
//     TH1D* hRatio, 
//     LCMSAxis axis, 
//     int centIdx, 
//     int ktIdx,
//     const std::string& outputDir = "."
// ) {
//     if (!hRatio) return;
    
//     // Создаём директорию если не существует
//     gSystem->mkdir(outputDir.c_str(), kTRUE);
    
//     // Имя файла без расширения
//     std::string baseName = FormatRatioName("ratio", centIdx, ktIdx, axis);
    
//     // Холст
//     TCanvas* c = new TCanvas(
//         Form("c_%s", baseName.c_str()),
//         Form("Ratio %s cent=%d kt=%d", ToString(axis).c_str(), centIdx, ktIdx),
//         800, 600
//     );
    
//     // Поля для читаемых подписей
//     c->SetBottomMargin(0.13);
//     c->SetLeftMargin(0.14);
//     c->SetRightMargin(0.05);
//     c->SetTopMargin(0.08);
    
//     // Рисуем гистограмму с ошибками
//     hRatio->Draw("P E");
    
//     // Горизонтальная линия на уровне 1.0
//     TLine* line = new TLine(
//         hRatio->GetXaxis()->GetXmin(), 1.0,
//         hRatio->GetXaxis()->GetXmax(), 1.0
//     );
//     line->SetLineColor(kRed+1);
//     line->SetLineStyle(7);  // пунктир
//     line->SetLineWidth(2);
//     line->Draw();
    
//     // Подписи на холсте
//     TLatex latex;
//     latex.SetNDC();
//     latex.SetTextFont(42);
//     latex.SetTextSize(0.045);
    
//     // Левый верхний угол: физические параметры
//     latex.SetTextAlign(11);
//     latex.DrawLatex(0.18, 0.88, ToString(axis).c_str());
//     latex.DrawLatex(0.18, 0.82, GetCentralityLabel(centIdx).c_str());
//     latex.DrawLatex(0.18, 0.76, GetKTRangeLabel(ktIdx).c_str());
    
//     // Правый верхний угол: коллаборация
//     latex.SetTextAlign(31);
//     latex.SetTextSize(0.042);
//     latex.SetTextColor(kGray+2);
//     latex.DrawLatex(0.95, 0.92, "ALICE Performance");
    
//     c->Update();
    
//     // Сохраняем в разные форматы
//     c->SaveAs(Form("%s/%s.png", outputDir.c_str(), baseName.c_str()));
//     c->SaveAs(Form("%s/%s.pdf", outputDir.c_str(), baseName.c_str()));
    
//     // Очищаем память
//     delete c;
// }

TH1D* Project1D(TH3D& h, LCMSAxis axis, double w = 0.05) {
    SetSlice1D(h, axis, w);

    const char* proj =
        axis==LCMSAxis::Out  ? "x" :
        axis==LCMSAxis::Side ? "y" : "z";

    TH1D* out = (TH1D*)h.Project3D(proj);
    out->SetDirectory(nullptr);
    return out;
}

void Style(TH1D* hRatio, LCMSAxis axis, int centIdx, int ktIdx) {
    if (!hRatio) return;

    TCanvas* c = new TCanvas(
        Form(
            "c_%s_%d_%d",
            ToString(axis).c_str(),
            centIdx,
            ktIdx
        ),
        Form(
            "Ratio %s cent=%d kt=%d", 
            ToString(axis).c_str(), 
            centIdx, 
            ktIdx
        ),
        800, 600
    );

    TPad* mainPad = new TPad("main", "main", 0.0, 0.0, 1.0, 1.0);
    mainPad->SetBottomMargin(0.13);
    mainPad->SetLeftMargin(0.14);
    mainPad->SetRightMargin(0.05);
    mainPad->SetTopMargin(0.08);
    mainPad->Draw();
    mainPad->cd();
    
    // Рисуем гистограмму
    hRatio->Draw("P E");  // маркеры с ошибками
    
    // Горизонтальная линия на уровне 1.0
    TLine* line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0,
                           hRatio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kRed);
    line->SetLineStyle(7);  // пунктир
    line->SetLineWidth(2);
    line->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.045);
    latex.SetTextFont(42);
    
    // Верхний левый угол: тип оси и диапазон kT
    std::string ktRange = GetKTRangeLabel(ktIdx);
    std::string centRange = GetCentralityLabel(centIdx);
    
    latex.DrawLatex(0.18, 0.88, Form("%s", ToString(axis).c_str()));
    latex.DrawLatex(0.18, 0.82, Form("%s", centRange.c_str()));
    latex.DrawLatex(0.18, 0.76, Form("%s", ktRange.c_str()));
    
    // Верхний правый угол: статистика (если нужно)
    // latex.SetTextAlign(31);
    // latex.DrawLatex(0.95, 0.88, "prelim");
    
    c->Update();
    // c->SaveAs(Form("ratio_%s_cent%d_kt%d.png", ToString(axis).c_str(), centIdx, ktIdx));
    // c->SaveAs(Form("ratio_%s_cent%d_kt%d.pdf", ToString(axis).c_str(), centIdx, ktIdx));
}

TH1D* ProjectRatio(TH3D& ratio3D, LCMSAxis axis, const char* name) {
    TH1D* r = Project1D(ratio3D, axis);
    
    r->SetName(name);

    return r;
}

TH1D* RatioProject(
    TH3D& Neg, TH3D& Pos,
    LCMSAxis axis, std::string name
) {
    TH1D* n = Project1D(Neg, axis);
    TH1D* p = Project1D(Pos, axis);

    auto ratio = RootPtr<TH1D>((TH1D*)n->Clone(name.c_str()));
    ratio->Divide(p, n);

    ApplyRatioStyle(ratio.get(), axis);
    return ratio.release();
}

void do_CF_ratios(TFile* fCF3D, TFile* fRatioProj, TFile* fProjRatio) {
    for (int centIdx = 0; centIdx < centralitySize; centIdx++)
    for (int ktIdx = 0; ktIdx < ktSize; ktIdx++) {
        TString nPos = getCFName(1, centIdx, ktIdx);
        TString nNeg = getCFName(0, centIdx, ktIdx);

        TH3D* hNeg = (TH3D*) fCF3D->Get(nNeg);
        TH3D* hPos = (TH3D*) fCF3D->Get(nPos);

        if (!hNeg || !hPos) continue;

        TH3D* neg = (TH3D*)hNeg->Clone(); neg->SetDirectory(nullptr);
        TH3D* pos = (TH3D*)hPos->Clone(); pos->SetDirectory(nullptr);

        TH3D* ratio = (TH3D*)neg->Clone(); ratio->SetDirectory(nullptr);
        ratio->Divide(pos, neg);

        for (auto ax : {LCMSAxis::Out, LCMSAxis::Side, LCMSAxis::Long}) {
            std::string rp_name = FormatRatioName("ratio_proj", centIdx, ktIdx, ax);
            TH1D* ratio_of_proj = RatioProject(*neg, *pos, ax, rp_name);

            fRatioProj->cd();
            ratio_of_proj->Write();

            std::string proj_of_ratio_name = Form("proj_of_ratios_%d_%d_%s", centIdx, ktIdx, ToString(ax).c_str());
            TH1D* proj_of_ratio = ProjectRatio(*ratio, ax, proj_of_ratio_name.c_str());

            fProjRatio->cd();
            proj_of_ratio->Write();
        }
    }
}