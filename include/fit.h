#pragma once

#include <TF3.h>
#include <TH3.h>

// фит-функция
Double_t CF_fit_3d(Double_t* q, Double_t* par);

// создаёт и настраивает
TF3* CreateCF3DFit();


// результат фита
struct FitResult {
    std::array<double,3> R{};
    std::array<double,3> eR{};
    double lambda{};
    double elambda{};
    bool ok{false};
};

// Фит + получение парамметров 
FitResult FitCF3D(TH3D* hCF, TF3* fit3d);
