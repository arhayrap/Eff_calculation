#include <cstdlib>

#include <vector>

#include <iostream>

#include <map>

#include <string>

#include <math.h>

#include "TFile.h"

#include "TTree.h"

#include "TString.h"

#include "TSystem.h"

#include "TStyle.h"

#include "TROOT.h"

#include "TStopwatch.h"

#include "TMVA/Tools.h"

#include "TMVA/Reader.h"

#include "TMVA/MethodCuts.h"

#include <chrono>

#include <filesystem>

namespace fs = std::filesystem;

typedef std::vector<std::string> stringvec;

// D-phi
double phi_dist(double a, double b) {
  if (fabs(a - b) > 3.14159265) {
    return 6.2831853 - fabs(a - b);
  }
  return fabs(a - b);
}

// D-R
double dr_Func(double eta1, double eta2, double phidist12) {
  return sqrt((eta1 - eta2) * (eta1 - eta2) + phidist12 * phidist12);
}

using namespace std;
using namespace std::chrono;

void Make_Histograms() {
  TString dataset = "DATA_PL";
  TFile * f = new TFile("Histograms_" + dataset + ".root", "RECREATE");

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

  Float_t cutflow_bins[] = {
    0.0,
    1.0,
    2.0,
    3.0,
    4.0,
    5.0,
    6.0,
    7.0,
    8.0,
    9.0
  };
  
  auto start = high_resolution_clock::now();

  int N_cutflow_bins = sizeof(cutflow_bins) / sizeof(cutflow_bins[0]);
  const int x_bins = sizeof(probe_pt_range) / sizeof(probe_pt_range[0]) - 1;
  const int y_bins = sizeof(displacement_significance) / sizeof(displacement_significance[0]) - 1;

  TH1D * _Mu_Tag_pt = new TH1D("Mu_Tag_pt", "Mu_Tag_pt", 100, 0, 100);
  TH1D * _Mu_Tag_eta = new TH1D("Mu_Tag_eta", "Mu_Tag_eta", 100, -5, 5);
  TH1D * _Mu_Probe_pt = new TH1D("Mu_Probe_pt", "Mu_Probe_pt", 200, 0, 100);
  TH1D * _Mu_Probe_eta = new TH1D("Mu_Probe_eta", "Mu_Probe_eta", 100, -5, 5);
  TH1D * _Mu_Probe_IP = new TH1D("Mu_Probe_IP", "Mu_Probe_IP", 100, 0, 100);
  TH1D * _Inv_Mass = new TH1D("Inv_Mass", "Inv_Mass", 100, 2.5, 3.5);
  TH1D * _DR_Tag_Probe = new TH1D("DR_Tag_Probe", "DR_Tag_Probe", 30, 0, 1);
  TH1D * _DPHI_Tag_Probe = new TH1D("DPHI_Tag_Probe", "DPHI_Tag_Probe", 30, 0, 4);
  TH1D * IP = new TH1D("Probe_IP", "Probe_IP", 150, 0, 500);

  TH1D * PROB_SV = new TH1D("PROB_SV", "PROB_SV", 30, 0, 1);
  TH1D * TAG_DR_SV = new TH1D("TAG_DR_SV", "TAG_DR_SV", 30, -150, 150);
  TH1D * PROBE_DR_SV = new TH1D("PROBE_DR_SV", "PROBE_DR_SV", 30, -150, 150);
  TH1D * cutflow_table = new TH1D("cutflow_table", "cutflow_table", 20, 0, 10); // cutflow_bins);
  cutflow_table -> GetXaxis() -> SetNdivisions(10);

  TH2D * _Mu_Probe_2D_PT_IP = new TH2D("Mu_Probe_2D_PT_IP", "Mu_Probe_2D_PT_IP", 200, 0, 100, 1000, 0, 500);
  TH2D * _Mu_Probe_2D_PT_IP_GRID = new TH2D("Mu_Probe_2D_PT_IP_GRID", "Mu_Probe_2D_PT_IP_GRID", x_bins, probe_pt_range, y_bins, displacement_significance);
  TH2D * _Tag_Probe_DR_SV = new TH2D("_Tag_Probe_DR_SV", "_Tag_Probe_DR_SV", 30, 0, 1.5, 30, 0, 1.5);
  // TH2D * _Tag_Probe_DR_SV = new TH2D("_Tag_Probe_DR_SV", "_Tag_Probe_DR_SV", 30, 0, 1.5, 30, 0, 1.5);
  TH2D * _Tag_Probe_charge = new TH2D("_Tag_Probe_charge", "_Tag_Probe_charge", 10, -5, 5, 10, -5, 5);

  TH1D * _Mu_Probe_Matched_pt = new TH1D("Mu_Probe_Matched_pt", "Mu_Probe_Matched_pt", 100, 0, 100);
  TH1D * _Mu_Probe_Matched_IP = new TH1D("Mu_Probe_Matched_IP", "Mu_Probe_Matched_IP", 500, 0, 500);
  TH2D * _Mu_Probe_Matched_2D_PT_IP = new TH2D("Mu_Probe_Matched_2D_PT_IP", "Mu_Probe_Matched_2D_PT_IP", 200, 0, 100, 1000, 0, 500);
  TH2D * _Mu_Probe_Matched_2D_MAP_NUM = new TH2D("Mu_Probe_Matched_2D_MAP_NUM", "Mu_Probe_Matched_2D_MAP_NUM", 10, 0, 10, 9, 0, 9);
  TH2D * _Mu_Probe_Matched_2D_MAP_NOM = new TH2D("Mu_Probe_Matched_2D_MAP_NOM", "Mu_Probe_Matched_2D_MAP_NUM", 10, 0, 10, 9, 0, 9);
  TH2D * _Mu_Probe_Matched_2D_MAP = new TH2D("Mu_Probe_Matched_2D_MAP_NOM", "Mu_Probe_Matched_2D_MAP_NUM", 10, 0, 10, 9, 0, 9);
  TH2D * _Mu_Probe_Matched_2D_PT_IP_GRID = new TH2D("Mu_Probe_Matched_2D_PT_IP_GRID", "Mu_Probe_Matched_2D_PT_IP_GRID", x_bins, probe_pt_range, y_bins, displacement_significance);
  TH2D * _Mu_Probe_Matched_2D_PT_IP_GRID_EFF = new TH2D("Mu_Probe_Matched_2D_PT_IP_GRID_EFF", "Mu_Probe_Matched_2D_PT_IP_GRID_EFF", x_bins, probe_pt_range, y_bins, displacement_significance);

  TFile * f0;
  TChain * ch;
  std::string base_folder = "/eos/user/a/arhayrap/";
  std::string base = "";

  if (dataset == "MC_PL") {
    f0 = TFile::Open("/eos/user/a/arhayrap/JPsi_output_data_PL_Tight_Loose_ID_test_old_version/JPsi_output_data_PL_Tight_Loose_ID_test_old_version.root");
  } else if (dataset == "MC_UL") {
    f0 = TFile::Open("/eos/user/a/arhayrap/JPsi_output_data_UL_Tight_Loose_ID_test_old_version/JPsi_output_data_UL_Tight_Loose_ID_test_old_version.root");
  } else {
    f0 = TFile::Open("/eos/user/a/arhayrap/BParking_skimmed_A_1_UL/BParking_skimmed_A_1_UL.root");
  }
  
  /*std::string bases[] = {
      // "BParking_skimmed_A_1_PL",
      // "BParking_skimmed_A_2_PL",
      // "BParking_skimmed_A_3_PL",
      // "BParking_skimmed_A_4_PL",
      // "BParking_skimmed_A_5_PL",
      // "BParking_skimmed_A_6_PL"
  };*/
  
  std::string bases[] = {
      "BParking_skimmed_A_1_UL",
      "BParking_skimmed_A_2_UL",
      "BParking_skimmed_A_3_UL",
      "BParking_skimmed_A_4_UL",
      "BParking_skimmed_A_5_UL",
      "BParking_skimmed_A_6_UL"
  };
  
  TChain * T_Tag = new TChain("Training/Reco_Jets/Tag_Candidate");
  TChain * T_Probe = new TChain("Training/Reco_Jets/Probe_Candidate");
  TChain * T_Tag_tr = new TChain("Training/Reco_Jets/Tag_Candidate_Track");
  TChain * T_Probe_tr = new TChain("Training/Reco_Jets/Probe_Candidate_Track");
  TChain * T_L3_0 = new TChain("Training/Reco_Jets/L3_0");
  TChain * T_L3_1 = new TChain("Training/Reco_Jets/L3_1");
  TChain * T_L3_2 = new TChain("Training/Reco_Jets/L3_2");
  TChain * T_L3_3 = new TChain("Training/Reco_Jets/L3_3");
  TChain * HLT_result = new TChain("Training/Reco_Jets/HLT_result");

  for (int b = 0; b < sizeof(bases) / sizeof(bases[0]); b++) {
      std::string path = base_folder + "/" + bases[b] + "/";
      for (const auto & entry : fs::directory_iterator(path)) {
          cout<<entry.path().string()<<endl;
          TString entry_path = entry.path().string();
          T_Tag->AddFile(entry_path);
          T_Probe->AddFile(entry_path);
          T_Tag_tr->AddFile(entry_path);
          T_Probe_tr->AddFile(entry_path);
          T_L3_0->AddFile(entry_path);
          T_L3_1->AddFile(entry_path);
          T_L3_2->AddFile(entry_path);
          T_L3_3->AddFile(entry_path);
          HLT_result->AddFile(entry_path);
      }
  }

  /*
  TTree * T_Tag = (TTree * ) f0 -> Get("Training/Reco_Jets/Tag_Candidate");
  TTree * T_Probe = (TTree * ) f0 -> Get("Training/Reco_Jets/Probe_Candidate");
  TTree * T_Tag_tr = (TTree * ) f0 -> Get("Training/Reco_Jets/Tag_Candidate_Track");
  TTree * T_Probe_tr = (TTree * ) f0 -> Get("Training/Reco_Jets/Probe_Candidate_Track");
  TTree * T_L3_0 = (TTree * ) f0 -> Get("Training/Reco_Jets/L3_0");
  TTree * T_L3_1 = (TTree * ) f0 -> Get("Training/Reco_Jets/L3_1");
  TTree * T_L3_2 = (TTree * ) f0 -> Get("Training/Reco_Jets/L3_2");
  TTree * T_L3_3 = (TTree * ) f0 -> Get("Training/Reco_Jets/L3_3");
  TTree * HLT_result = (TTree * ) f0 -> Get("Training/Reco_Jets/HLT_result");
  */

  Float_t weight;
  Float_t NPV;
  Float_t Tag_pt, Tag_e, Tag_eta, Tag_phi, Tag_charge, Tag_dxy, Tag_dxyErr, Tag_dz, Tag_LooseID, Tag_TightID, Tag_Global, Tag_PV_theta;
  Float_t Probe_pt, Probe_e, Probe_eta, Probe_phi, Probe_charge, Probe_dxy, Probe_dxyErr, Probe_dz, Probe_LooseID, Probe_TightID, Probe_Global, Probe_PV_theta;

  Float_t Tag_nstrips, Tag_ntracks, Tag_npixels;
  Float_t Tag_track_purity;
  Float_t Tag_probSV;
  Float_t Tag_dr_SV;

  Float_t Probe_nstrips, Probe_ntracks, Probe_npixels;
  Float_t Probe_track_purity;
  Float_t Probe_probSV;
  Float_t Probe_dr_SV;

  Float_t Prob_SV;

  Float_t L3_0_pt, L3_0_eta, L3_0_phi;
  Float_t L3_1_pt, L3_1_eta, L3_1_phi;
  Float_t L3_2_pt, L3_2_eta, L3_2_phi;
  Float_t L3_3_pt, L3_3_eta, L3_3_phi;

  Float_t L3_0_HLT_Mu12_IP6, L3_0_HLT_Mu10p5_IP3p5, L3_0_HLT_Mu9_IP6, L3_0_HLT_Mu9_IP5, L3_0_HLT_Mu9_IP4, L3_0_HLT_Mu8p5_IP3p5, L3_0_HLT_Mu8_IP5, L3_0_HLT_Mu8_IP6, L3_0_HLT_Mu8_IP3, L3_0_HLT_Mu7_IP4;
  Float_t L3_1_HLT_Mu12_IP6, L3_1_HLT_Mu10p5_IP3p5, L3_1_HLT_Mu9_IP6, L3_1_HLT_Mu9_IP5, L3_1_HLT_Mu9_IP4, L3_1_HLT_Mu8p5_IP3p5, L3_1_HLT_Mu8_IP5, L3_1_HLT_Mu8_IP6, L3_1_HLT_Mu8_IP3, L3_1_HLT_Mu7_IP4;
  Float_t L3_2_HLT_Mu12_IP6, L3_2_HLT_Mu10p5_IP3p5, L3_2_HLT_Mu9_IP6, L3_2_HLT_Mu9_IP5, L3_2_HLT_Mu9_IP4, L3_2_HLT_Mu8p5_IP3p5, L3_2_HLT_Mu8_IP5, L3_2_HLT_Mu8_IP6, L3_2_HLT_Mu8_IP3, L3_2_HLT_Mu7_IP4;
  Float_t L3_3_HLT_Mu12_IP6, L3_3_HLT_Mu10p5_IP3p5, L3_3_HLT_Mu9_IP6, L3_3_HLT_Mu9_IP5, L3_3_HLT_Mu9_IP4, L3_3_HLT_Mu8p5_IP3p5, L3_3_HLT_Mu8_IP5, L3_3_HLT_Mu8_IP6, L3_3_HLT_Mu8_IP3, L3_3_HLT_Mu7_IP4;

  T_Tag -> SetBranchAddress("Gweight", & weight);
  T_Tag -> SetBranchAddress("Tnpv", & NPV);

  T_Tag -> SetBranchAddress("mu_pt", & Tag_pt);
  T_Tag -> SetBranchAddress("mu_energy", & Tag_e);
  T_Tag -> SetBranchAddress("mu_eta", & Tag_eta);
  T_Tag -> SetBranchAddress("mu_phi", & Tag_phi);
  T_Tag -> SetBranchAddress("mu_charge", & Tag_charge);
  T_Tag -> SetBranchAddress("mu_dxy", & Tag_dxy);
  T_Tag -> SetBranchAddress("mu_dxyErr", & Tag_dxyErr);
  T_Tag -> SetBranchAddress("mu_dz", & Tag_dz);
  T_Tag -> SetBranchAddress("mu_LooseID", & Tag_LooseID);
  T_Tag -> SetBranchAddress("mu_TightID", & Tag_TightID);
  T_Tag -> SetBranchAddress("mu_isGlobal", & Tag_Global);

  T_Tag_tr -> SetBranchAddress("n_strip_layers", & Tag_nstrips);
  T_Tag_tr -> SetBranchAddress("n_pixel_layers", & Tag_npixels);
  T_Tag_tr -> SetBranchAddress("n_track_layers", & Tag_ntracks);
  T_Tag_tr -> SetBranchAddress("mu_PV_theta", & Tag_PV_theta);
  T_Tag_tr -> SetBranchAddress("track_purity", & Tag_track_purity);
  T_Tag_tr -> SetBranchAddress("probSV", & Tag_probSV);
  T_Tag_tr -> SetBranchAddress("probSV", & Tag_dr_SV);

  T_Probe -> SetBranchAddress("mu_pt", & Probe_pt);
  T_Probe -> SetBranchAddress("mu_energy", & Probe_e);
  T_Probe -> SetBranchAddress("mu_eta", & Probe_eta);
  T_Probe -> SetBranchAddress("mu_phi", & Probe_phi);
  T_Probe -> SetBranchAddress("mu_charge", & Probe_charge);
  T_Probe -> SetBranchAddress("mu_dxy", & Probe_dxy);
  T_Probe -> SetBranchAddress("mu_dxyErr", & Probe_dxyErr);
  T_Probe -> SetBranchAddress("mu_dz", & Probe_dz);
  T_Probe -> SetBranchAddress("mu_LooseID", & Probe_LooseID);
  T_Probe -> SetBranchAddress("mu_TightID", & Probe_TightID);
  T_Probe -> SetBranchAddress("mu_isGlobal", & Probe_Global);

  T_Probe_tr -> SetBranchAddress("n_strip_layers", & Probe_nstrips);
  T_Probe_tr -> SetBranchAddress("n_pixel_layers", & Probe_npixels);
  T_Probe_tr -> SetBranchAddress("n_track_layers", & Probe_ntracks);
  T_Probe_tr -> SetBranchAddress("mu_PV_theta", & Probe_PV_theta);
  T_Probe_tr -> SetBranchAddress("track_purity", & Probe_track_purity);
  T_Probe_tr -> SetBranchAddress("probSV", & Probe_probSV);
  T_Probe_tr -> SetBranchAddress("probSV", & Probe_dr_SV);

  T_L3_0 -> SetBranchAddress("mu_pt", & L3_0_pt);
  T_L3_0 -> SetBranchAddress("mu_eta", & L3_0_eta);
  T_L3_0 -> SetBranchAddress("mu_phi", & L3_0_phi);
  T_L3_0 -> SetBranchAddress("HLT_Mu12_IP6", & L3_0_HLT_Mu12_IP6);
  T_L3_0 -> SetBranchAddress("HLT_Mu10p5_IP3p5", & L3_0_HLT_Mu10p5_IP3p5);
  T_L3_0 -> SetBranchAddress("HLT_Mu9_IP6", & L3_0_HLT_Mu9_IP6);
  T_L3_0 -> SetBranchAddress("HLT_Mu9_IP5", & L3_0_HLT_Mu9_IP5);
  T_L3_0 -> SetBranchAddress("HLT_Mu9_IP4", & L3_0_HLT_Mu9_IP4);
  T_L3_0 -> SetBranchAddress("HLT_Mu8p5_IP3p5", & L3_0_HLT_Mu8p5_IP3p5);
  T_L3_0 -> SetBranchAddress("HLT_Mu8_IP5", & L3_0_HLT_Mu8_IP5);
  T_L3_0 -> SetBranchAddress("HLT_Mu8_IP6", & L3_0_HLT_Mu8_IP6);
  T_L3_0 -> SetBranchAddress("HLT_Mu8_IP3", & L3_0_HLT_Mu8_IP3);
  T_L3_0 -> SetBranchAddress("HLT_Mu7_IP4", & L3_0_HLT_Mu7_IP4);

  T_L3_1 -> SetBranchAddress("mu_pt", & L3_1_pt);
  T_L3_1 -> SetBranchAddress("mu_eta", & L3_1_eta);
  T_L3_1 -> SetBranchAddress("mu_phi", & L3_1_phi);
  T_L3_1 -> SetBranchAddress("HLT_Mu12_IP6", & L3_1_HLT_Mu12_IP6);
  T_L3_1 -> SetBranchAddress("HLT_Mu10p5_IP3p5", & L3_1_HLT_Mu10p5_IP3p5);
  T_L3_1 -> SetBranchAddress("HLT_Mu9_IP6", & L3_1_HLT_Mu9_IP6);
  T_L3_1 -> SetBranchAddress("HLT_Mu9_IP5", & L3_1_HLT_Mu9_IP5);
  T_L3_1 -> SetBranchAddress("HLT_Mu9_IP4", & L3_1_HLT_Mu9_IP4);
  T_L3_1 -> SetBranchAddress("HLT_Mu8p5_IP3p5", & L3_1_HLT_Mu8p5_IP3p5);
  T_L3_1 -> SetBranchAddress("HLT_Mu8_IP5", & L3_1_HLT_Mu8_IP5);
  T_L3_1 -> SetBranchAddress("HLT_Mu8_IP6", & L3_1_HLT_Mu8_IP6);
  T_L3_1 -> SetBranchAddress("HLT_Mu8_IP3", & L3_1_HLT_Mu8_IP3);
  T_L3_1 -> SetBranchAddress("HLT_Mu7_IP4", & L3_1_HLT_Mu7_IP4);

  T_L3_2 -> SetBranchAddress("mu_pt", & L3_2_pt);
  T_L3_2 -> SetBranchAddress("mu_eta", & L3_2_eta);
  T_L3_2 -> SetBranchAddress("mu_phi", & L3_2_phi);
  T_L3_2 -> SetBranchAddress("HLT_Mu12_IP6", & L3_2_HLT_Mu12_IP6);
  T_L3_1 -> SetBranchAddress("HLT_Mu10p5_IP3p5", & L3_1_HLT_Mu10p5_IP3p5);
  T_L3_2 -> SetBranchAddress("HLT_Mu9_IP6", & L3_2_HLT_Mu9_IP6);
  T_L3_2 -> SetBranchAddress("HLT_Mu9_IP5", & L3_2_HLT_Mu9_IP5);
  T_L3_2 -> SetBranchAddress("HLT_Mu9_IP4", & L3_2_HLT_Mu9_IP4);
  T_L3_2 -> SetBranchAddress("HLT_Mu8p5_IP3p5", & L3_2_HLT_Mu8p5_IP3p5);
  T_L3_2 -> SetBranchAddress("HLT_Mu8_IP5", & L3_2_HLT_Mu8_IP5);
  T_L3_2 -> SetBranchAddress("HLT_Mu8_IP6", & L3_2_HLT_Mu8_IP6);
  T_L3_2 -> SetBranchAddress("HLT_Mu8_IP3", & L3_2_HLT_Mu8_IP3);
  T_L3_2 -> SetBranchAddress("HLT_Mu7_IP4", & L3_2_HLT_Mu7_IP4);

  T_L3_3 -> SetBranchAddress("mu_pt", & L3_3_pt);
  T_L3_3 -> SetBranchAddress("mu_eta", & L3_3_eta);
  T_L3_3 -> SetBranchAddress("mu_phi", & L3_3_phi);
  T_L3_3 -> SetBranchAddress("HLT_Mu12_IP6", & L3_3_HLT_Mu12_IP6);
  T_L3_3 -> SetBranchAddress("HLT_Mu10p5_IP3p5", & L3_3_HLT_Mu10p5_IP3p5);
  T_L3_3 -> SetBranchAddress("HLT_Mu9_IP6", & L3_3_HLT_Mu9_IP6);
  T_L3_3 -> SetBranchAddress("HLT_Mu9_IP5", & L3_3_HLT_Mu9_IP5);
  T_L3_3 -> SetBranchAddress("HLT_Mu9_IP4", & L3_3_HLT_Mu9_IP4);
  T_L3_3 -> SetBranchAddress("HLT_Mu8p5_IP3p5", & L3_3_HLT_Mu8p5_IP3p5);
  T_L3_3 -> SetBranchAddress("HLT_Mu8_IP5", & L3_3_HLT_Mu8_IP5);
  T_L3_3 -> SetBranchAddress("HLT_Mu8_IP6", & L3_3_HLT_Mu8_IP6);
  T_L3_3 -> SetBranchAddress("HLT_Mu8_IP3", & L3_3_HLT_Mu8_IP3);
  T_L3_3 -> SetBranchAddress("HLT_Mu7_IP4", & L3_3_HLT_Mu7_IP4);

  HLT_result -> SetBranchAddress("Prob_SV", & Prob_SV);

    TString cuts[] = {
        "No cuts",
        "Charge sum",
        "Minv cut",
        "SV_matching",
        "Loose Prob(SV)",
        "Soft ID",
        "cos(#eta)"
    };
    
    const int N = sizeof(cuts) / sizeof(cuts[0]);

  TH2D * cutflow_plots[N];
  TH2D * cutflow_plots_denominator[N];

  for (int c = 0; c < N; c++) {
    cutflow_plots[c] = new TH2D(cuts[c], cuts[c], 10, 0, 10, 9, 0, 9);
    cutflow_plots_denominator[c] = new TH2D(cuts[c], cuts[c], 10, 0, 10, 9, 0, 9);
  }

  double n_matched = 0.0;
  double matched_denominator = 0.0;
  double test_num = 0;
  double test_nom = 0;
  for (int i = 0; i < T_Tag -> GetEntries(); i++) {
    T_Tag -> GetEntry(i);
    T_Probe -> GetEntry(i);
    T_Tag_tr -> GetEntry(i);
    T_Probe_tr -> GetEntry(i);
    T_L3_0 -> GetEntry(i);
    T_L3_1 -> GetEntry(i);
    T_L3_2 -> GetEntry(i);
    T_L3_3 -> GetEntry(i);
    HLT_result -> GetEntry(i);

    // if (i > 1000000) break;

    if (i % 100000 == 0) {
      // if (i % 100 == 0) {
      cout << i << " / " << T_Tag -> GetEntries() << endl;
    }

    // ***********************************  Selection 
    // if(Tag_pt < 12) continue; // TAG muon PT selection.

    Bool_t L3_0_HLT = (L3_0_HLT_Mu12_IP6 || L3_0_HLT_Mu10p5_IP3p5 || L3_0_HLT_Mu9_IP6 || L3_0_HLT_Mu9_IP5 || L3_0_HLT_Mu9_IP4 || L3_0_HLT_Mu8p5_IP3p5 || L3_0_HLT_Mu8_IP5 || L3_0_HLT_Mu8_IP6 || L3_0_HLT_Mu8_IP3 || L3_0_HLT_Mu7_IP4);
    Bool_t L3_1_HLT = (L3_1_HLT_Mu12_IP6 || L3_1_HLT_Mu10p5_IP3p5 || L3_1_HLT_Mu9_IP6 || L3_1_HLT_Mu9_IP5 || L3_1_HLT_Mu9_IP4 || L3_1_HLT_Mu8p5_IP3p5 || L3_1_HLT_Mu8_IP5 || L3_1_HLT_Mu8_IP6 || L3_1_HLT_Mu8_IP3 || L3_1_HLT_Mu7_IP4);

    // Bool_t L3_0_HLT = (L3_0_HLT_Mu12_IP6 || L3_0_HLT_Mu9_IP6 || L3_0_HLT_Mu9_IP5 || L3_0_HLT_Mu9_IP4 || L3_0_HLT_Mu8_IP5 || L3_0_HLT_Mu8_IP6 || L3_0_HLT_Mu8_IP3 || L3_0_HLT_Mu7_IP4);
    // Bool_t L3_1_HLT = (L3_1_HLT_Mu12_IP6 || L3_1_HLT_Mu9_IP6 || L3_1_HLT_Mu9_IP5 || L3_1_HLT_Mu9_IP4 || L3_1_HLT_Mu8_IP5 || L3_1_HLT_Mu8_IP6 || L3_1_HLT_Mu8_IP3 || L3_1_HLT_Mu7_IP4);
    
    bool matchet_to_leading = false;
    bool matchet_to_subleading = false;
    double Probe_IP = abs(Probe_dxy / Probe_dxyErr);

    if (dr_Func(Tag_eta, L3_0_eta, phi_dist(Tag_phi, L3_0_phi)) < 0.2) matchet_to_leading = true;
    if (dr_Func(Tag_eta, L3_1_eta, phi_dist(Tag_phi, L3_1_phi)) < 0.2) matchet_to_subleading = true;

    Bool_t L3_0_HLT_tag = (L3_0_HLT_Mu12_IP6 || L3_0_HLT_Mu9_IP6);
    Bool_t L3_1_HLT_tag = (L3_1_HLT_Mu12_IP6 || L3_1_HLT_Mu9_IP6);

    TLorentzVector p4_Tag, p4_Probe, DiMu_System;
    Float_t Tag_m = sqrt(pow(Tag_e, 2) - (pow(Tag_pt * cos(Tag_phi), 2) + pow(Tag_pt * sin(Tag_phi), 2) + pow(Tag_pt * sinh(Tag_eta), 2)));
    Float_t Probe_m = sqrt(pow(Probe_e, 2) - (pow(Probe_pt * cos(Probe_phi), 2) + pow(Probe_pt * sin(Probe_phi), 2) + pow(Probe_pt * sinh(Probe_eta), 2)));

    p4_Tag.SetPtEtaPhiM(Tag_pt, Tag_eta, Tag_phi, Tag_m);
    p4_Probe.SetPtEtaPhiM(Probe_pt, Probe_eta, Probe_phi, Probe_m);
    DiMu_System = p4_Tag + p4_Probe;

    // ========================================================================================================================== NO CUTS
    // ==========================================================================================================================
    // ==========================================================================================================================
    // ==========================================================================================================================

    if (Tag_TightID != 1) continue;
    if (Probe_LooseID != 1) continue;
    // if (!((matchet_to_leading && L3_0_HLT_tag) || (matchet_to_subleading_1 && L3_1_HLT_tag) || (matchet_to_subleading_2 && L3_2_HLT_tag) || (matchet_to_subleading_3 && L3_3_HLT_tag))) continue; // The TAG muon had matched to a leading or subleading trigger object.
    if (!((matchet_to_leading && L3_0_HLT_tag) || (matchet_to_subleading && L3_1_HLT_tag))) continue; // The TAG muon had matched to a leading or subleading trigger object.
    // if (!((L3_0_HLT_tag) || (L3_1_HLT_tag))) continue; // The TAG muon had matched to a leading or subleading trigger object.
    // if (!(Tag_dr_SV >= 0 && Probe_dr_SV >= 0)) continue;
    
    // if (!((Tag_charge + Probe_charge) == 0)) continue;
    
    double dr_Tag_Probe = dr_Func(Tag_eta, Probe_eta, phi_dist(Tag_phi, Probe_phi));
    double dphi_Tag_Probe = phi_dist(Tag_phi, Probe_phi);

    _Tag_Probe_DR_SV -> Fill(Tag_dr_SV, Probe_dr_SV);
    _Tag_Probe_charge -> Fill(Tag_charge, Probe_charge);
    TAG_DR_SV -> Fill(Tag_dr_SV);
    PROBE_DR_SV -> Fill(Probe_dr_SV);
    
    bool selection[N] = {};
    bool sel = false;
    
    selection[0] = true;
    selection[1] = ((Tag_charge + Probe_charge) == 0);
    selection[2] = ((DiMu_System.M() > 2.9) && (DiMu_System.M() < 3.3));
    selection[3] = ((Tag_dr_SV < 0.4) && (Tag_dr_SV > 0.0) && (Probe_dr_SV < 0.4) && (Probe_dr_SV > 0.0));
    selection[4] = ((Prob_SV > 1e-5));
    // selection[4] = ((Prob_SV > 0.1));
    selection[5] = ((Tag_Global) && (Tag_ntracks > 5) && (Tag_npixels > 0) && (Tag_dxy < 0.3) && (Tag_dz < 20) && (Tag_track_purity == 1.0));
    selection[6] = ((cos(Tag_PV_theta) > 0.9) && (cos(Probe_PV_theta) > 0.9));
    // selection[6] = ((cos(Tag_PV_theta + Probe_PV_theta) > 0.9));

    bool all_selections = false;

    for (int s = 0; s < N; s++) {
        if (s == 0) {
            all_selections = selection[0];
        } else {
            all_selections = (all_selections && selection[s]);
        }
    }

    for (int iter = 0; iter < N; iter++) {
      // Fill the selection bins that correspond to the iteration number.
      if (iter == 0) {
        sel = selection[0];
      } else {
        for (int s = 0; s <= iter; s++) {
          sel = (sel && selection[s]);
        }
      }
      if (sel) {
        cutflow_table -> Fill(iter);
        for (int j = 0; j < x_bins; j++) { // PT
          for (int k = 0; k < y_bins; k++) { // Displacement significance
            if (Probe_pt >= probe_pt_range[j] && Probe_pt < probe_pt_range[j + 1] && Probe_IP >= displacement_significance[k] && Probe_IP < displacement_significance[k + 1]) {
              cutflow_plots_denominator[iter] -> SetBinContent(j + 1, k + 1, cutflow_plots_denominator[iter] -> GetBinContent(j + 1, k + 1) + 1.0);
            }
          }
        }
        if (matchet_to_leading && L3_1_HLT) {
          if (dr_Func(Probe_eta, L3_1_eta, phi_dist(Probe_phi, L3_1_phi)) < 0.2) {
            // cout<<"matched to a subleading trigger object"<<endl;
            for (int j = 0; j < x_bins; j++) { // PT
              for (int k = 0; k < y_bins; k++) { // Displacement significance
                if (Probe_pt >= probe_pt_range[j] && Probe_pt < probe_pt_range[j + 1] && Probe_IP >= displacement_significance[k] && Probe_IP < displacement_significance[k + 1]) {
                  cutflow_plots[iter] -> SetBinContent(j + 1, k + 1, cutflow_plots[iter] -> GetBinContent(j + 1, k + 1) + 1.0);
                }
              }
            }
          }
        }
        else if (matchet_to_subleading && L3_0_HLT) {
          if (dr_Func(Probe_eta, L3_0_eta, phi_dist(Probe_phi, L3_0_phi)) < 0.2) {
            // cout<<"matched to a leading trigger object"<<endl;
            for (int j = 0; j < x_bins; j++) { // PT
              for (int k = 0; k < y_bins; k++) { // Displacement significance
                if (Probe_pt >= probe_pt_range[j] && Probe_pt < probe_pt_range[j + 1] && Probe_IP >= displacement_significance[k] && Probe_IP < displacement_significance[k + 1]) {
                  cutflow_plots[iter] -> SetBinContent(j + 1, k + 1, cutflow_plots[iter] -> GetBinContent(j + 1, k + 1) + 1.0);
                }
              }
            }
          }
        }
      }
    }
    
    
    if (!all_selections) continue;

    // if (!(Tag_charge + Probe_charge == 0)) continue; // Summary charge of the TAG and PROBE muons is 0.
    // if (!((L3_0_HLT_tag && matchet_to_leading) || (L3_1_HLT_tag && matchet_to_subleading))) continue; // Matching to one of the leading or subleading HLT trigger objects.
    // if (!((Tag_Global) && (Tag_ntracks > 5) && (Tag_npixels > 0) && (Tag_dxy < 0.3) && (Tag_dz < 20) && (Tag_track_purity == 1.0))) continue; // Soft muon identification criteria
    // if (!(Prob_SV > 1e-5)) continue;
    // if (!(cos(Tag_PV_theta) > 0.9 && cos(Probe_PV_theta) > 0.9)) continue;
    // if (!(DiMu_System.M() > 2.9 && DiMu_System.M() < 3.3)) continue;
    // if (!(Tag_track_purity == 1.0 && Probe_track_purity == 1.0)) continue;
    // if (dr_Func(Tag_eta, Probe_eta, phi_dist(Tag_phi, Probe_phi))>0.6) continue;

    test_nom+=1.0;
    IP -> Fill(Probe_IP);
    _Inv_Mass -> Fill(DiMu_System.M());
    _Mu_Tag_pt -> Fill(Tag_pt);
    _Mu_Tag_eta -> Fill(Tag_eta);
    _Mu_Probe_pt -> Fill(Probe_pt);
    _Mu_Probe_eta -> Fill(Probe_eta);
    _DR_Tag_Probe -> Fill(dr_Tag_Probe);
    _DPHI_Tag_Probe -> Fill(dphi_Tag_Probe);
    PROB_SV -> Fill(Prob_SV);

    _Mu_Probe_IP -> Fill(Probe_IP);
    _Mu_Probe_2D_PT_IP -> Fill(Probe_pt, Probe_IP);
    _Mu_Probe_2D_PT_IP_GRID -> Fill(Probe_pt, Probe_IP);

    if (L3_1_pt > 0 && matchet_to_leading && L3_1_HLT) {
      test_num+=1.0;
      // cutflow_table->Fill(9);
      if (dr_Func(Probe_eta, L3_1_eta, phi_dist(Probe_phi, L3_1_phi)) < 0.2) {
        _Mu_Probe_Matched_pt -> Fill(Probe_pt);
        _Mu_Probe_Matched_IP -> Fill(Probe_IP);
        for (int i = 0; i < x_bins; i++) { // PT
          for (int j = 0; j < y_bins; j++) { // Displacement significance
            if (Probe_pt > probe_pt_range[i] && Probe_pt < probe_pt_range[i + 1] && Probe_IP > displacement_significance[j] && Probe_IP < displacement_significance[j + 1]) {
              _Mu_Probe_Matched_2D_MAP_NUM -> SetBinContent(i + 1, j + 1, _Mu_Probe_Matched_2D_MAP_NUM -> GetBinContent(i + 1, j + 1) + 1);
            }
          }
        }
        _Mu_Probe_Matched_2D_PT_IP -> Fill(Probe_pt, Probe_IP);
        _Mu_Probe_Matched_2D_PT_IP_GRID -> Fill(Probe_pt, Probe_IP);
       } else {
           cout<<matchet_to_leading<<"       "<<L3_1_HLT<<"       "<<phi_dist(Probe_phi, L3_1_phi)<<endl;
       }
    } else if (L3_0_pt > 0 && matchet_to_subleading && L3_0_HLT) {
      test_num+=1.0;
      // cutflow_table->Fill(9);
       if (dr_Func(Probe_eta, L3_0_eta, phi_dist(Probe_phi, L3_0_phi)) < 0.2) {
        _Mu_Probe_Matched_pt -> Fill(Probe_pt);
        _Mu_Probe_Matched_IP -> Fill(Probe_IP);
        for (int i = 0; i < x_bins; i++) { // PT
          for (int j = 0; j < y_bins; j++) { // Displacement significance
            if (Probe_pt > probe_pt_range[i] && Probe_pt < probe_pt_range[i + 1] && Probe_IP > displacement_significance[j] && Probe_IP < displacement_significance[j + 1]) {
              _Mu_Probe_Matched_2D_MAP_NUM -> SetBinContent(i + 1, j + 1, _Mu_Probe_Matched_2D_MAP_NUM -> GetBinContent(i + 1, j + 1) + 1);
            }
          }
        }
        _Mu_Probe_Matched_2D_PT_IP -> Fill(Probe_pt, Probe_IP);
        _Mu_Probe_Matched_2D_PT_IP_GRID -> Fill(Probe_pt, Probe_IP);
      } else {
        cout<<matchet_to_subleading<<"       "<<L3_0_HLT<<"      "<<phi_dist(Probe_phi, L3_0_phi)<<endl;
      }
    }
    for (int i = 0; i < x_bins; i++) { // PT
       for (int j = 0; j < y_bins; j++) { // Displacement significance
        if (Probe_pt > probe_pt_range[i] && Probe_pt < probe_pt_range[i + 1] && Probe_IP > displacement_significance[j] && Probe_IP < displacement_significance[j + 1]) {
          _Mu_Probe_Matched_2D_MAP_NOM -> SetBinContent(i + 1, j + 1, _Mu_Probe_Matched_2D_MAP_NOM -> GetBinContent(i + 1, j + 1) + 1);
        }
      }
    }

  }

  for (int j = 0; j < x_bins; j++) {
    for (int k = 0; k < y_bins; k++) {
      Float_t Nom = _Mu_Probe_2D_PT_IP_GRID -> GetBinContent(j + 1, k + 1);
      Float_t Nom_err = _Mu_Probe_2D_PT_IP_GRID -> GetBinError(j + 1, k + 1);
      for (int l = 0; l < N; l++) {
        Float_t Num_1 = cutflow_plots[l] -> GetBinContent(j + 1, k + 1);
        Float_t Num_err_1 = cutflow_plots[l] -> GetBinError(j + 1, k + 1);
        Float_t Nom_1 = cutflow_plots_denominator[l] -> GetBinContent(j + 1, k + 1);
        Float_t Nom_err_1 = cutflow_plots_denominator[l] -> GetBinError(j + 1, k + 1);
        cutflow_plots[l] -> SetBinContent(j + 1, k + 1, Num_1 / Nom_1);
        cutflow_plots[l] -> SetBinError(j + 1, k + 1, (Num_1 / Nom_1) * sqrt(pow(Num_err_1 / Num_1, 2) + pow(Nom_err_1 / Nom_1, 2)));
      }

      Float_t Num = _Mu_Probe_Matched_2D_PT_IP_GRID -> GetBinContent(j + 1, k + 1);
      Float_t Num_err = _Mu_Probe_Matched_2D_PT_IP_GRID -> GetBinError(j + 1, k + 1);

      _Mu_Probe_Matched_2D_PT_IP_GRID_EFF -> SetBinContent(j + 1, k + 1, Num / Nom);
      _Mu_Probe_Matched_2D_PT_IP_GRID_EFF -> SetBinError(j + 1, k + 1, (Num / Nom) * sqrt(pow(Num_err / Num, 2) + pow(Nom_err / Nom, 2)));

      Num = _Mu_Probe_Matched_2D_MAP_NUM -> GetBinContent(j + 1, k + 1);
      Nom = _Mu_Probe_Matched_2D_MAP_NOM -> GetBinContent(j + 1, k + 1);
      Nom_err = _Mu_Probe_Matched_2D_MAP_NOM -> GetBinError(j + 1, k + 1);
      Num_err = _Mu_Probe_Matched_2D_MAP_NUM -> GetBinError(j + 1, k + 1);

      _Mu_Probe_Matched_2D_MAP -> SetBinContent(j + 1, k + 1, Num / Nom);
      _Mu_Probe_Matched_2D_MAP -> SetBinError(j + 1, k + 1, (Num / Nom) * sqrt(pow(Num_err / Num, 2) + pow(Nom_err / Nom, 2)));
    }
  }

  for (int j = 0; j < N; j++) {
    cutflow_table -> GetXaxis() -> ChangeLabel(j + 1, 60, 0.025, -1, -1, -1, cuts[j]);
  }

  f -> Write();

  TCanvas * Minv = new TCanvas("Minv", "Minv");
  _Inv_Mass -> Draw();
  Minv -> SaveAs("Minv.png");

  gStyle -> SetOptFile(0);
  gStyle -> SetOptStat(0);
  gStyle -> SetPaintTextFormat("2.3f");
  int w = 1200;
  int h = 1100;
  TCanvas * HIST_2D_PT_IP = new TCanvas("HIST_2D_PT_IP", "HIST_2D_PT_IP", w, h);
  _Mu_Probe_Matched_2D_PT_IP -> Draw("COLZ");
  HIST_2D_PT_IP -> SaveAs("HIST_2D_PT_IP.png");

  TCanvas * HIST_Tag_Probe_DR_SV = new TCanvas("HIST_Tag_Probe_DR_SV", "HIST_Tag_Probe_DR_SV", w, h);
  _Tag_Probe_DR_SV -> Draw("COLZ");
  HIST_Tag_Probe_DR_SV -> SaveAs("HIST_2D_DR_SV.png");

  TCanvas * HIST_Tag_Probe_charge = new TCanvas("HIST_Tag_Probe_charge", "HIST_Tag_Probe_charge", w, h);
  _Tag_Probe_charge -> Draw("COLZ TEXT");
  HIST_Tag_Probe_charge -> SaveAs("HIST_Tag_Probe_charge.png");

  TCanvas * HIST_2D_PT_IP_GRID = new TCanvas("HIST_2D_PT_IP_GRID", "HIST_2D_PT_IP_GRID", w, h);
  _Mu_Probe_Matched_2D_PT_IP_GRID -> Draw("COLZ TEXT");
  HIST_2D_PT_IP_GRID -> SaveAs("HIST_2D_PT_IP_GRID.png");

  TCanvas * HIST_DR_Tag_Probe = new TCanvas("HIST_DR_Tag_Probe", "HIST_DR_Tag_Probe", w, h);
  _DR_Tag_Probe -> Draw();
  HIST_DR_Tag_Probe -> SaveAs("HIST_DR_Tag_Probe.png");

  TCanvas * HIST_DPHI_Tag_Probe = new TCanvas("HIST_DPHI_Tag_Probe", "HIST_DPHI_Tag_Probe", w, h);
  _DPHI_Tag_Probe -> Draw();
  HIST_DPHI_Tag_Probe -> SaveAs("HIST_DPHI_Tag_Probe.png");

  TCanvas * HIST_2D_PT_IP_GRID_EFF = new TCanvas("HIST_2D_PT_IP_GRID_EFF", "HIST_2D_PT_IP_GRID_EFF", w, h);
  _Mu_Probe_Matched_2D_PT_IP_GRID_EFF -> Draw("COLZ TEXT");
  HIST_2D_PT_IP_GRID_EFF -> SetLogx();
  HIST_2D_PT_IP_GRID_EFF -> SetLogy();
  HIST_2D_PT_IP_GRID_EFF -> SaveAs("HIST_2D_PT_IP_GRID_EFF.png");

  TCanvas * HIST_Mu_Probe_Matched_2D_MAP = new TCanvas("_Mu_Probe_Matched_2D_MAP", "_Mu_Probe_Matched_2D_MAP", w, h);
  _Mu_Probe_Matched_2D_MAP -> GetZaxis() -> SetRangeUser(0.0, 1.0);
  _Mu_Probe_Matched_2D_MAP -> Draw("COLZ TEXT");
  HIST_Mu_Probe_Matched_2D_MAP -> SaveAs("HIST_Mu_Probe_Matched_2D_MAP.png");

  TCanvas * HIST_PROB_SV = new TCanvas("HIST_PROB_SV", "HIST_PROB_SV", w, h);
  PROB_SV -> Draw();
  HIST_PROB_SV -> SaveAs("PROB_SV.png");

  TCanvas * HIST_TAG_DR_SV = new TCanvas("HIST_TAG_DR_SV", "HIST_TAG_DR_SV", w, h);
  TAG_DR_SV -> Draw();
  HIST_TAG_DR_SV -> SaveAs("HIST_TAG_DR_SV.png");

  TCanvas * HIST_PROBE_DR_SV = new TCanvas("HIST_PROBE_DR_SV", "HIST_PROBE_DR_SV", w, h);
  PROBE_DR_SV -> Draw();
  HIST_PROBE_DR_SV -> SaveAs("HIST_PROBE_DR_SV.png");

  TCanvas * HIST_IP = new TCanvas("HIST_IP", "HIST_IP", w, h);
  IP -> Draw();
  HIST_IP -> SaveAs("HIST_IP.png");

  TCanvas * CUTFLOW_TABLE = new TCanvas("CUTFLOW_TABLE", "CUTFLOW_TABLE", w, h);
  cutflow_table -> Draw("hist text");
  CUTFLOW_TABLE -> SaveAs("CUTFLOW_TABLE.png");

  for (int p = 0; p < N; p++) {
    TCanvas * cutflow_canvas = new TCanvas(cuts[p], cuts[p], w, h);
    cutflow_plots[p] -> GetZaxis() -> SetRangeUser(0, 1.0);
    cutflow_plots[p] -> Draw("COLZ TEXT E");
    cutflow_canvas -> SaveAs(cuts[p] + "_2D_HIST.png");
  }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  cout<<"The duration of the execution of the script: "<<duration.count() * 1000000<<" s."<<endl;
}