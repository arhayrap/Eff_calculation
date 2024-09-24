#include "HttStyles.cc"
#include "TStyle.h"
#include "TGaxis.h"
#include "TRandom.h"

#include <iostream>
#include <math.h>
#include <TF1.h>
#include <TH1D.h>
#include "TCanvas.h"

void TwoD_Compare_Eff()
{

// *******************************************************************

  TFile *f = new TFile("PT_IP_2D_Eff.root","RECREATE");

  TH2D *_PT_IP_2D_Eff    = new TH2D("PT_IP_2D_Eff","PT_IP_2D_Eff",100,0,100,10,0,100);

  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0);
  SetStyle();
  gStyle->SetPaintTextFormat("2.3f");

  TString PlotTitle="";

  // TCanvas *c1 = new TCanvas(PlotTitle, PlotTitle,100,52,900,686);
  TCanvas *c1 = new TCanvas(PlotTitle, PlotTitle,1500,1100);
//  c1->SetGridy();
  c1->Range(0,0,1,1);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(10);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetFrameFillStyle(0);
  c1->SetFrameLineStyle(0);
  c1->SetFrameLineWidth(2);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderSize(10);
  c1->SetBottomMargin(0.2);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.15);

  // TFile* f_MC     = new TFile("Histograms_BParking.root");
  // TFile* f_MC     = new TFile("Histograms_MC_UL.root");
  // TFile* f_DATA   = new TFile("Histograms_DATA_UL.root");

  // TFile* f_MC     = new TFile("Histograms_MC_UL_111.root");
  // TFile* f_DATA   = new TFile("Histograms_DATA_UL_111.root");
  // TFile* f_DATA   = new TFile("Histograms.root");

  TFile* f_MC     = new TFile("./tmp_7/Histograms_MC_UL_111.root");
  TFile* f_DATA   = new TFile("./tmp_7/Histograms_DATA_UL_111.root");

  TH2D * _MC_Probe           = (TH2D*)f_MC->Get("Mu_Probe_2D_PT_IP");
  TH2D * _MC_Probe_Matched   = (TH2D*)f_MC->Get("Mu_Probe_Matched_2D_PT_IP");
  TH2D * _DATA_Probe         = (TH2D*)f_DATA->Get("Mu_Probe_2D_PT_IP");
  TH2D * _DATA_Probe_Matched = (TH2D*)f_DATA->Get("Mu_Probe_Matched_2D_PT_IP");

  Float_t probe_pt_range[]            = {6.0, 7.0, 8.0, 8.5, 9.0, 10.0, 10.5, 11.0, 12.0, 20.0, 100.0};
  Float_t displacement_significance[] = {0.0, 3.0, 3.5, 4.0, 5.0,  6.0,  8.0, 10.0, 20.0, 500.0};
  
  //Float_t probe_pt_range[]            = {6.0, 7.0, 8.0, 8.5, 9.0, 10.0, 10.5, 11.0, 12.0, 20.0, 100.0};
  //  Float_t displacement_significance[] = {0.0, 3.0, 3.5, 4.0, 5.0,  6.0,  8.0, 10.0, 20.0, 500.0};
  //Float_t displacement_significance[] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5}; //, 500.0};

  const int x_bins = sizeof(probe_pt_range) / sizeof(probe_pt_range[0]) - 1;
  const int y_bins = sizeof(displacement_significance) / sizeof(displacement_significance[0]) - 1;
  cout<<x_bins<<"   "<<y_bins<<endl;
  
  // TH2D * _MC_Probe           = (TH2D*)f_MC->Get("Mu_Probe_2D_PT_IP");
  // TH2D * _MC_Probe_Matched   = (TH2D*)f_MC->Get("Mu_Probe_Matched_2D_PT_IP");
  
  TH2F * Num_MC = new TH2F("MC_numerator", "MC_numerator", x_bins, 0, x_bins, y_bins, 0, y_bins);
  TH2F * Nom_MC = new TH2F("MC_denominator", "MC_denominator", x_bins, 0, x_bins, y_bins, 0, y_bins);
  TH2D * Eff_MC = new TH2D("Efficiency MC PL", "Efficiency MC PL", x_bins, 0, x_bins, y_bins, 0, y_bins);

  TH2F * Num_DATA = new TH2F("DATA_numerator", "DATA_numerator", x_bins, 0, x_bins, y_bins, 0, y_bins);
  TH2F * Nom_DATA = new TH2F("DATA_denominator", "DATA_denominator", x_bins, 0, x_bins, y_bins, 0, y_bins);
  TH2D * Eff_DATA = new TH2D("Efficiency BParking PL", "Efficiency BParking PL", x_bins, 0, x_bins, y_bins, 0, y_bins);

  TH2D * Trigger_SF = new TH2D("Trigger_SF", "Trigger_SF", x_bins, 0, x_bins, y_bins, 0, y_bins);

  gStyle->SetPalette(kBird);
  gStyle->SetMarkerSize(0.5);
  // gStyle->SetTextSize(0.01);
// ******************************************************************* MC
  _MC_Probe_Matched->SetMarkerStyle(23);
  _MC_Probe_Matched->GetXaxis()->SetLabelFont(42);
  _MC_Probe_Matched->GetXaxis()->SetLabelOffset(0.02);
  _MC_Probe_Matched->GetXaxis()->SetTitle("Mu - PT");
  _MC_Probe_Matched->GetXaxis()->SetTitleSize(0.045);
  _MC_Probe_Matched->GetXaxis()->SetTitleOffset(1.2);
  _MC_Probe_Matched->GetXaxis()->SetTitleFont(42);
  _MC_Probe_Matched->GetYaxis()->SetNdivisions(510);
  _MC_Probe_Matched->GetXaxis()->SetNdivisions(505);
  _MC_Probe_Matched->GetYaxis()->SetLabelFont(42);
  _MC_Probe_Matched->GetYaxis()->SetLabelSize(0.04);
  _MC_Probe_Matched->GetYaxis()->SetLabelOffset(0.01);
  _MC_Probe_Matched->GetYaxis()->SetTitle("Mu - IP");
  _MC_Probe_Matched->GetYaxis()->SetTitleOffset(1.5);
  _MC_Probe_Matched->GetYaxis()->SetTitleSize(0.045);
  _MC_Probe_Matched->GetYaxis()->SetTitleFont(42);
  _MC_Probe_Matched->SetTitle(PlotTitle);
  _MC_Probe_Matched->GetZaxis()->SetTitle(" ");

  _MC_Probe->Sumw2();
  _MC_Probe_Matched->Sumw2();
  
  cout<<"Before the loop"<<endl;
  int N_bins_x = Num_MC->GetNbinsX();
  int N_bins_y = Nom_MC->GetNbinsY();
  float step_x = 0.5; // 100.0 / N_bins_x;
  float step_y = 0.5; // 0.025; // 500.0 / N_bins_y;
  for (int i = 0; i < x_bins; i++) { // PT
    for (int j = 0; j < y_bins; j++) { // Displacement significance
      for (int k = 0; k < 200; k++) { // PT
       for (int l = 0; l < 1000; l++) { // Displacement significance
        //for (int k = 0; k < 100; k++) { // PT
          //for (int l = 0; l < 50; l++) { // Displacement significance
            // cout<<k*step_x<<"   "<<l*step_y<<endl;
          if ((k*step_x >= probe_pt_range[i]) && (k*step_x < probe_pt_range[i+1]) && (l*step_y >= displacement_significance[j]) && (l*step_y < displacement_significance[j+1])) {
            Num_MC->SetBinContent(i+1, j+1, Num_MC->GetBinContent(i+1, j+1) + _MC_Probe_Matched->GetBinContent(k+1, l+1));
            Nom_MC->SetBinContent(i+1, j+1, Nom_MC->GetBinContent(i+1, j+1) + _MC_Probe        ->GetBinContent(k+1, l+1));
            Num_DATA->SetBinContent(i+1, j+1, Num_DATA->GetBinContent(i+1, j+1) + _DATA_Probe_Matched->GetBinContent(k+1, l+1));
            Nom_DATA->SetBinContent(i+1, j+1, Nom_DATA->GetBinContent(i+1, j+1) + _DATA_Probe        ->GetBinContent(k+1, l+1));
          }
        }
      }
      cout<<i+1<<"  "<<j+1<<")    "<<Num_MC->GetBinContent(i+1, j+1)<<"    "<<Nom_MC->GetBinContent(i+1, j+1)<<"    "<<Num_MC->GetBinContent(i+1, j+1) / Nom_MC->GetBinContent(i+1, j+1)<<endl;
      cout<<i+1<<"  "<<j+1<<")    "<<Num_DATA->GetBinContent(i+1, j+1)<<"    "<<Nom_DATA->GetBinContent(i+1, j+1)<<"    "<<Num_DATA->GetBinContent(i+1, j+1) / Nom_DATA->GetBinContent(i+1, j+1)<<endl;
      // Eff_MC->SetBinContent(i+1, j+1, Num_MC->GetBinContent(i+1, j+1) / Nom_MC->GetBinContent(i+1, j+1));
      float sigma_num = sqrt(Num_MC->GetBinContent(i+1, j+1));
      float sigma_nom = sqrt(Nom_MC->GetBinContent(i+1, j+1));
      float N_num = Num_MC->GetBinContent(i+1, j+1);
      float N_nom = Nom_MC->GetBinContent(i+1, j+1);
      Eff_MC->SetBinContent(i+1, j+1, (N_num / N_nom));
      Eff_MC->SetBinError(i+1, j+1, (N_num / N_nom) * sqrt(pow(sigma_num / N_num, 2) + pow(sigma_nom / N_nom, 2)));
      sigma_num = sqrt(Num_DATA->GetBinContent(i+1, j+1));
      sigma_nom = sqrt(Nom_DATA->GetBinContent(i+1, j+1));
      N_num = Num_DATA->GetBinContent(i+1, j+1);
      N_nom = Nom_DATA->GetBinContent(i+1, j+1);
      Eff_DATA->SetBinContent(i+1, j+1, (N_num / N_nom));
      Eff_DATA->SetBinError(i+1, j+1, (N_num / N_nom) * sqrt(pow(sigma_num / N_num, 2) + pow(sigma_nom / N_nom, 2)));
      float eff_MC = Eff_MC->GetBinContent(i+1, j+1);
      float eff_DATA = Eff_DATA->GetBinContent(i+1, j+1);
      float eff_MC_err = Eff_MC->GetBinError(i+1, j+1);
      float eff_DATA_err = Eff_DATA->GetBinError(i+1, j+1);
      float ratio = (eff_DATA / eff_MC);
      cout<<"SF: "<<eff_DATA<<"    "<<eff_MC<<"    "<<(ratio)<<endl;
      Trigger_SF->SetBinContent(i+1, j+1, ratio);
      Trigger_SF->SetBinError(i+1, j+1, (ratio) * sqrt(pow(eff_MC_err / eff_MC, 2) + pow(eff_DATA_err / eff_DATA, 2)));
    }
  }

  cout<<"After the loop"<<endl;

  TCanvas * c_eff_MC = new TCanvas("Efficiency MC", "Efficiency MC", 1000, 1000);
  
  for (int i = 0; i <= x_bins; i++) {
    for (int j = 0; j <= y_bins; j++) {
      char label[5];
      sprintf(label, "%.2f", displacement_significance[j]);
      Eff_MC->GetYaxis()->ChangeLabel(j + 1, 0, 0.03, -1, -1, -1, label);
    }
    char label[5];
    sprintf(label, "%.2f", probe_pt_range[i]);
    Eff_MC->GetXaxis()->ChangeLabel(i + 1, 60, 0.03, -1, -1, -1, label);
  }

  // Eff_MC->GetZaxis()->SetRangeUser(0, 1.0);
  c_eff_MC->SetLeftMargin(0.15);
  c_eff_MC->SetRightMargin(0.15);
  Eff_MC->SetMarkerSize(1.0);
  Eff_MC->Draw("colz text e");
  c_eff_MC->SaveAs("EFF_MC.png");

  TCanvas * c_eff_DATA = new TCanvas("Efficiency DATA", "Efficiency DATA", 1000, 1000);
  
  for (int i = 0; i <= x_bins; i++) {
    for (int j = 0; j <= y_bins; j++) {
      char label[5];
      sprintf(label, "%.2f", displacement_significance[j]);
      Eff_DATA->GetYaxis()->ChangeLabel(j + 1, 0, 0.03, -1, -1, -1, label);
    }
    char label[5];
    sprintf(label, "%.2f", probe_pt_range[i]);
    Eff_DATA->GetXaxis()->ChangeLabel(i + 1, 60, 0.03, -1, -1, -1, label);
  }

  // Eff_DATA->GetZaxis()->SetRangeUser(0, 1.0);
  c_eff_DATA->SetLeftMargin(0.15);
  c_eff_DATA->SetRightMargin(0.15);
  Eff_DATA->SetMarkerSize(1.0);
  Eff_DATA->Draw("colz text e");
  c_eff_DATA->SaveAs("EFF_Data.png");

  TFile * BParking_triggers = new TFile("BParking_triggers_PL.root");

  TCanvas * C_trigger_SF = new TCanvas("Trigger scale factors PL", "Trigger scale factors PL", 1000, 1000);
  
  for (int i = 0; i <= x_bins; i++) {
    for (int j = 0; j <= y_bins; j++) {
      char label[5];
      sprintf(label, "%.2f", displacement_significance[j]);
      Trigger_SF->GetYaxis()->ChangeLabel(j + 1, 0, 0.03, -1, -1, -1, label);
    }
    char label[5];
    sprintf(label, "%.2f", probe_pt_range[i]);
    Trigger_SF->GetXaxis()->ChangeLabel(i + 1, 60, 0.03, -1, -1, -1, label);
  }
  
  // C_trigger_SF->GetZaxis()->SetRangeUser(0, 1.0);
  C_trigger_SF->SetLeftMargin(0.15);
  C_trigger_SF->SetRightMargin(0.15);
  Trigger_SF->SetMarkerSize(1.0);
  Trigger_SF->Draw("colz text e");
  C_trigger_SF->SaveAs("Trigger_SF.png");
  
  BParking_triggers->Close();
  
  _MC_Probe_Matched->Divide(_MC_Probe);
  _MC_Probe_Matched->SetMaximum(1.0);
  _MC_Probe_Matched->SetMinimum(0.0);
  _MC_Probe_Matched->Draw("text colz");

   c1->Modified();
   c1->cd();


// ************************************* Fill SFs Histo 

  f->Write();

}