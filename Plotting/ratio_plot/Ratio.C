#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"

void Ratio() {

    TFile * f_PL = TFile::Open("BParking_SF_PL_2024.root");
    TFile * f_UL = TFile::Open("BParking_SF_UL.root");
    TH2D * PL_2D = (TH2D * ) f_PL -> Get("BParking_trigger_efficiency");
    TH2D * UL_2D = (TH2D * ) f_UL -> Get("BParking_trigger_efficiency");
    
    Float_t probe_pt_range[] = {
    6.0,
    7.0,
    8.0,
    8.5,
    9.0,
    10.0,
    10.5,
    11.0,
    12.0,
    20.0,
    100.0
    };

    Float_t displacement_significance[] = {
    0.0,
    3.0,
    3.5,
    4.0,
    5.0,
    6.0,
    8.0,
    10.0,
    20.0,
    500.0
    };
    
    const int x_bins = sizeof(probe_pt_range) / sizeof(probe_pt_range[0]) - 1;
    const int y_bins = sizeof(displacement_significance) / sizeof(displacement_significance[0]) - 1;
    
    gStyle->SetOptStat(0000);
    
    TH2D * Ratio_2D = new TH2D("PU 2024 / UL ratio", "PU 2024 / UL ratio", 10, 0, 10, 9, 0, 9); // x_bins, probe_pt_range, y_bins, displacement_significance);
    TCanvas * c = new TCanvas("PU 2024 / UL ratio", "PU 2024 / UL ratio", 1000, 700);
    
    for (int i = 0; i <= x_bins; i++) {
        for (int j = 0; j <= y_bins; j++) {
            double num   = PL_2D->GetBinContent(i, j);
            double denom = UL_2D->GetBinContent(i, j);
            double num_err   = PL_2D->GetBinError(i, j);
            double denom_err = UL_2D->GetBinError(i, j);
            cout<<i<<"   "<<j<<endl;
            if (num != 0.0 && denom != 0.0) {
                double ratio = (num / denom);
                double error = ratio * (pow(num_err/num, 2) + pow(denom_err/denom, 2));
                cout<<ratio<<endl;
                Ratio_2D->SetBinContent(i, j, ratio);
                Ratio_2D->SetBinError(i, j, error);
            } else {
                Ratio_2D->SetBinContent(i, j, 0.0);
                Ratio_2D->SetBinError(i, j, 0.0);
            }
        }
    }
    
    for (int i = 0; i <= x_bins; i++) {
        for (int j = 0; j <= y_bins; j++) {
            char label[5];
            sprintf(label, "%.2f", displacement_significance[j]);
            Ratio_2D->GetYaxis()->ChangeLabel(j + 1, 0, 0.03, -1, -1, -1, label);
        }
        char label[5];
        sprintf(label, "%.2f", probe_pt_range[i]);
        Ratio_2D->GetXaxis()->ChangeLabel(i + 1, 60, 0.03, -1, -1, -1, label);
    }
    
    Ratio_2D->GetZaxis()->SetRangeUser(0.0, 10.0);
    Ratio_2D->SetTitle("Ratio of the scale factors PL / UL");
    Ratio_2D->Draw("colz text e");
    c->SaveAs("PL_2024_UL_SF_ratio.png");

}