#include <iostream>
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iomanip>
#include <math.h>
#include "TStyle.h"

using namespace std;

// enum ParIndex_t {
//    Bkg0=0, Bkg1=1, Bkg2,
//    SigScale, SigSigma, SigMean,
//    N_PAR};

// const std::map<ParIndex_t,std::string> parNames{
//    {SigScale, "Gauss scale"}, {SigSigma, "Gauss #sigma"}, {SigMean, "Gauss #mu"}
// };

// double gaussian(double *x, double *par) {
//     return par[0] * TMath::Gaus(x[0], par[SigMean], par[SigSigma], true) + x[0];
// }

// double background(double *x, double *par) {
//     return par[Bkg0]
// }

void SF_calculation() {

    TFile * file_data = TFile::Open("../Condor/BParking_skimmed_A_PL_111.root");
    TFile * file_mc   = TFile::Open("../Condor/JPsi_PL_list.root");
    system("mkdir -p FitPlots");
    
            
    // for (int i = 5; i <= 5; i++) {
    //     for (int j = 4; j <= 4; j++) {

    for (int i = 1; i <= 10; i++) {
        for (int j = 1; j <= 9; j++) {
            
            TString name = "Minv_" + to_string(i) + "_" + to_string(j) + "_num";
            TH1D * file_data_num_inv_mass   = (TH1D * ) file_data -> Get(name);
            name = "Minv_" + to_string(i) + "_" + to_string(j) + "_denom";
            TH1D * file_data_denom_inv_mass = (TH1D * ) file_data -> Get(name);
            name = "Minv_" + to_string(i) + "_" + to_string(j) + "_num";
            TH1D * file_mc_num_inv_mass   = (TH1D * ) file_mc -> Get(name);
            name = "Minv_" + to_string(i) + "_" + to_string(j) + "_denom";
            TH1D * file_mc_denom_inv_mass = (TH1D * ) file_mc -> Get(name);
            
            file_data_num_inv_mass->GetXaxis()->SetRangeUser(2.9, 3.3);
            file_data_denom_inv_mass->GetXaxis()->SetRangeUser(2.9, 3.3);
            file_mc_num_inv_mass->GetXaxis()->SetRangeUser(2.9, 3.3);
            file_mc_denom_inv_mass->GetXaxis()->SetRangeUser(2.9, 3.3);
            
            cout<<"Denominator mean and STD: "<<file_data_denom_inv_mass->GetMean()<<"   "<<file_data_denom_inv_mass->GetRMS()<<endl;
            cout<<"Numerator mean and STD: "<<file_data_num_inv_mass->GetMean()<<"   "<<file_data_num_inv_mass->GetRMS()<<endl;
            
            TF1 * gaussianFit_num  = new TF1("gaus_num", "[0]*TMath::Exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x", 2.9, 3.3); // [0]*TMath::Exp(-0.5*((x-[1])/[2])^2)
            gaussianFit_num->SetParLimits(0, file_data_num_inv_mass->GetMaximum()*0.5, file_data_num_inv_mass->GetMaximum()*1.5);
            gaussianFit_num->SetParLimits(1, 3.1 - 1*file_data_num_inv_mass->GetRMS(), 3.1 + 1*file_data_num_inv_mass->GetRMS());
            gaussianFit_num->SetParLimits(2, file_data_num_inv_mass->GetRMS()*0.5, file_data_num_inv_mass->GetRMS()*1.5);
            gaussianFit_num->SetParameter(0, file_data_num_inv_mass->GetMaximum());
            gaussianFit_num->SetParameter(1, file_data_num_inv_mass->GetMean());
            gaussianFit_num->SetParameter(2, file_data_num_inv_mass->GetRMS());
            
            TF1 * gaussianFit_denom  = new TF1("gaus_denom", "[0]*TMath::Exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x", 2.9, 3.3);
            gaussianFit_denom->SetParLimits(0, file_data_denom_inv_mass->GetMaximum()*0.5, file_data_denom_inv_mass->GetMaximum()*1.5);
            gaussianFit_denom->SetParLimits(1, 3.1 - 1*file_data_denom_inv_mass->GetRMS(), 3.1 + 1*file_data_denom_inv_mass->GetRMS());
            gaussianFit_denom->SetParLimits(2, file_data_denom_inv_mass->GetRMS()*0.5, file_data_denom_inv_mass->GetRMS()*1.5);
            gaussianFit_denom->SetParameter(0, file_data_denom_inv_mass->GetMaximum());
            gaussianFit_denom->SetParameter(1, file_data_denom_inv_mass->GetMean());
            gaussianFit_denom->SetParameter(2, file_data_denom_inv_mass->GetRMS());

            file_data_denom_inv_mass -> Fit (gaussianFit_denom, "R");
            file_data_num_inv_mass -> Fit (gaussianFit_num, "R");

            // double c_denom = gaussianFit_denom->GetParameter(3);
            // double c_num = gaussianFit_num->GetParameter(3);

            TF1 * chebyshevFit_num = new TF1("chebyshevFit_num", "[0]+[1]*x", 2.9, 3.3);
            cout<<gaussianFit_num->GetParameter(3)<<"   "<<gaussianFit_num->GetParameter(4)<<endl;
            chebyshevFit_num -> SetParameter(0, gaussianFit_num->GetParameter(3));
            chebyshevFit_num -> SetParameter(1, gaussianFit_num->GetParameter(4));
            cout<<endl;
            TF1 * chebyshevFit_denom = new TF1("chebyshevFit_denom", "[0]+[1]*x", 2.9, 3.3);
            cout<<gaussianFit_denom->GetParameter(3)<<"   "<<gaussianFit_denom->GetParameter(4)<<endl;
            chebyshevFit_denom -> SetParameter(0, gaussianFit_denom->GetParameter(3));
            chebyshevFit_denom -> SetParameter(1, gaussianFit_denom->GetParameter(4));
            
            // file_data_denom_inv_mass -> Fit (chebyshevFit_denom, "R");
            // file_data_num_inv_mass -> Fit (chebyshevFit_num, "R");
            
            TString   s = "Bin_" + to_string(i) + "_" + to_string(j);
            TString output_name = s + ".png";
            TCanvas * c = new TCanvas(s, s, 2400, 2000);
            c->Divide(2,2);
            
            c->cd(1);
            file_data_num_inv_mass->SetMarkerColor(kBlack);
            file_data_num_inv_mass->SetMarkerStyle(20);
            file_data_num_inv_mass->Draw("o E");
            gaussianFit_num->SetLineColor(kGreen);
            gaussianFit_num->SetLineWidth(4);
            gaussianFit_num->Draw("SAME");
            chebyshevFit_num->SetLineWidth(4);
            chebyshevFit_num->SetLineColor(kGreen);
            chebyshevFit_num->Draw("SAME");
            // c->cd(2);
            
            c->cd(3);
            file_data_denom_inv_mass->SetMarkerColor(kBlack);
            file_data_denom_inv_mass->SetMarkerStyle(20);
            file_data_denom_inv_mass->Draw("o E");
            gaussianFit_denom->SetLineColor(kBlue);
            gaussianFit_denom->SetLineWidth(4);
            gaussianFit_denom->Draw("SAME");
            chebyshevFit_denom->SetLineWidth(4);
            chebyshevFit_denom->SetLineColor(kBlue);
            chebyshevFit_denom->Draw("SAME");
            // c->cd(4);
            
            c->SaveAs("./FitPlots/" + output_name);
        }
    }

}