#include "../interface/TrigScaleFactors_OnBParking.h"
#include "../interface/SpecFunc.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <TROOT.h>
#include <TSystem.h>
#include "TFile.h"
#include <TStyle.h>
#include "TMath.h"
#include "TNtuple.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <math.h>

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/one/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include "FWCore/Framework/interface/Run.h"
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"


#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/Event.h>
#include "FWCore/Framework/interface/Run.h"
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/L1Trigger/interface/Jet.h"

#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

double TrigScaleFactors_OnBParking::calculateMuonSVProbability(const pat::Muon& muon, const std::vector<reco::Vertex>& vertices) {
    // Placeholder for probability calculation (replace with your implementation)
    double probability = 0.0;

    // Get the inner track of the muon
    const reco::TrackRef innerTrack = muon.innerTrack();

    if (innerTrack.isNonnull()) {
        // Check the association of the inner track with the vertices
        size_t numAssociatedVertices = 0;
        int index = 0;
        for (const auto& vertex : vertices) {
            // cout<<"Track ("<<index<<") "<<vertex.trackWeight(innerTrack)<<endl;
            if (vertex.trackWeight(innerTrack) > 0.5) { // Consider track with weight > 0.5 associated
                numAssociatedVertices++;
            }
        }

        // Probability calculation: Example only (replace with your implementation)
        probability = static_cast<double>(numAssociatedVertices) / static_cast<double>(vertices.size());
    }
    // cout<<"probability: "<<probability<<endl;
    return probability;
}

double TrigScaleFactors_OnBParking::calculateBackpointingAngle(const reco::Vertex& primaryVertex, const reco::VertexCompositePtrCandidate& secondaryVertex, const pat::Muon& muon) {
    // Vector from primary to secondary vertex
    // cout<<"1: just entered the function"<<endl;
    math::XYZVector displacement(
        secondaryVertex.vertex().x() - primaryVertex.x(),
        secondaryVertex.vertex().y() - primaryVertex.y(),
        0.0
    );
    // secondaryVertex.vertex().z() - primaryVertex.z()
    // cout<<"2: after making the displacement vertex"<<endl;

    // Momentum vector
    // math::XYZVector momentumVec(muon.px(), muon.py()); //, muon.pz());
    math::XYZVector momentumVec(muon.px(), muon.py(), 0.0); //, muon.pz());

    // cout<<"3: after making the momentum vector"<<endl;

    // Calculate dot product
    double dotProduct = displacement.Dot(momentumVec);
    // cout<<"4: dotProduct"<<endl;

    // Calculate magnitudes
    double displacementMag = displacement.R(); // sqrt(pow(displacement.X(),2) + pow(displacement.Y(),2) + pow(displacement.Z(),2));
    
    cout<<"Magnitude: "<<displacementMag<<endl;
    cout<<"Magnitude R: "<<displacement.R()<<endl;
    // cout<<"5: displacementMag"<<endl;

    double momentumMag = momentumVec.R();
    // cout<<"6: momentumMag"<<endl;

    // Calculate cosine of the angle
    double cosTheta = dotProduct / (displacementMag * momentumMag);
    // cout<<"7: cosTheta"<<endl;

    // Calculate the angle in radians
    double angle = std::acos(cosTheta);
    // cout<<"8: angle"<<endl;

    return angle;
}

double TrigScaleFactors_OnBParking::DeltaR(double eta1, double eta2, double phi1, double phi2) {
    return sqrt(pow(eta1 - eta2, 2) + pow(MyFunc::phi_dist(phi1, phi2), 2));
}

///////////////////////////////////////////////////////////////////

TrigScaleFactors_OnBParking::TrigScaleFactors_OnBParking(const edm::ParameterSet& iConfig):
   // genInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("GenInf"))),
   triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
   triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone> >(iConfig.getParameter<edm::InputTag>("objects"))),
   triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
   l1GtToken_(consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("l1GtSrc"))),
   l1jetToken_(consumes<BXVector<l1t::Jet>>(iConfig.getParameter<edm::InputTag>("l1jetSrc"))),
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("offlinePrimaryVertices"))),
   // vtxToken_(consumes<reco::Vertex>(iConfig.getParameter<edm::InputTag>("offlinePrimaryVertices"))),
   // secondary_vtxToken_(consumes<std::vector<reco::VertexCompositePtrCandidate>>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))),
   // secondary_vtxToken_(iConfig.getParameter<edm::InputTag>("secondary_vertices")),
   // consumes<std::vector<reco::VertexCompositePtrCandidate>>(iConfig.getParameter<edm::InputTag>("secondary_vertices")),
   secondaryVerticesToken_(consumes<std::vector<reco::VertexCompositePtrCandidate>>(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),
   puToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupInfo"))),
   recjetToken_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("recJet"))),
   recmuonToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("recMuon"))){
   runBegin = -1;
   lumibegin = 0;
   lumiend = 0;
   edm::Service<TFileService> fs;
   bookHistograms();
  }

TrigScaleFactors_OnBParking::~TrigScaleFactors_OnBParking(){}

// Function to calculate the backpointing angle


//void TrigScaleFactors_OnBParking::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {

//}

void TrigScaleFactors_OnBParking::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace pat;
  using namespace std;
  using namespace MyFunc;

  int lumi = iEvent.luminosityBlock();
  if (runBegin < 0) {lumibegin = lumiend = lumi;runBegin = iEvent.id().run();}
  if (lumi < lumibegin) lumibegin = lumi;
  if (lumi > lumiend)   lumiend = lumi;

// Generator weights *********************************************************
  // Handle<GenEventInfoProduct> genEvtInfo;
  // iEvent.getByToken(genInfoToken_, genEvtInfo);
  double Gweight=1.0; // genEvtInfo->weight();
  _Tot_Wgt->Fill(Gweight);

// HLT Control and Signal TriggerResults **************************************
   Handle<edm::TriggerResults> triggerBits;
   iEvent.getByToken(triggerBits_, triggerBits);
   const edm::TriggerNames & trigNames = iEvent.triggerNames(*triggerBits);
   string path_Mu9_IP5_part0   = "HLT_Mu9_IP5_part0_v", path_Mu8_IP5_part0   = "HLT_Mu8_IP5_part0_v";
   string path_Mu9_IP5_part1   = "HLT_Mu9_IP5_part1_v", path_Mu8_IP5_part1   = "HLT_Mu8_IP5_part1_v";
   string path_Mu9_IP5_part2   = "HLT_Mu9_IP5_part2_v", path_Mu8_IP5_part2   = "HLT_Mu8_IP5_part2_v";
   string path_Mu9_IP5_part3   = "HLT_Mu9_IP5_part3_v", path_Mu8_IP5_part3   = "HLT_Mu8_IP5_part3_v";
   string path_Mu9_IP5_part4   = "HLT_Mu9_IP5_part4_v", path_Mu8_IP5_part4   = "HLT_Mu8_IP5_part4_v";
   string path_Mu9_IP5_part5   = "HLT_Mu9_IP5_part5_v", path_Mu8_IP5_part5   = "HLT_Mu8_IP5_part5_v";

   string path_Mu9_IP4_part0   = "HLT_Mu9_IP4_part0_v", path_Mu8_IP6_part0   = "HLT_Mu8_IP6_part0_v";
   string path_Mu9_IP4_part1   = "HLT_Mu9_IP4_part1_v", path_Mu8_IP6_part1   = "HLT_Mu8_IP6_part1_v";
   string path_Mu9_IP4_part2   = "HLT_Mu9_IP4_part2_v", path_Mu8_IP6_part2   = "HLT_Mu8_IP6_part2_v";
   string path_Mu9_IP4_part3   = "HLT_Mu9_IP4_part3_v", path_Mu8_IP6_part3   = "HLT_Mu8_IP6_part3_v";
   string path_Mu9_IP4_part4   = "HLT_Mu9_IP4_part4_v", path_Mu8_IP6_part4   = "HLT_Mu8_IP6_part4_v";
   string path_Mu9_IP4_part5   = "HLT_Mu9_IP4_part5_v", path_Mu8_IP6_part5   = "HLT_Mu8_IP6_part5_v";

   string path_Mu9_IP6_part0   = "HLT_Mu9_IP6_part0_v", path_Mu8_IP3_part0   = "HLT_Mu8_IP3_part0_v";
   string path_Mu9_IP6_part1   = "HLT_Mu9_IP6_part1_v", path_Mu8_IP3_part1   = "HLT_Mu8_IP3_part1_v";
   string path_Mu9_IP6_part2   = "HLT_Mu9_IP6_part2_v", path_Mu8_IP3_part2   = "HLT_Mu8_IP3_part2_v";
   string path_Mu9_IP6_part3   = "HLT_Mu9_IP6_part3_v", path_Mu8_IP3_part3   = "HLT_Mu8_IP3_part3_v";
   string path_Mu9_IP6_part4   = "HLT_Mu9_IP6_part4_v", path_Mu8_IP3_part4   = "HLT_Mu8_IP3_part4_v";
   string path_Mu9_IP6_part5   = "HLT_Mu9_IP6_part5_v", path_Mu8_IP3_part5   = "HLT_Mu8_IP3_part5_v";

   string path_Mu12_IP6_part0  = "HLT_Mu12_IP6_part0_v";
   string path_Mu12_IP6_part1  = "HLT_Mu12_IP6_part1_v";
   string path_Mu12_IP6_part2  = "HLT_Mu12_IP6_part2_v";
   string path_Mu12_IP6_part3  = "HLT_Mu12_IP6_part3_v";
   string path_Mu12_IP6_part4  = "HLT_Mu12_IP6_part4_v";
   string path_Mu12_IP6_part5  = "HLT_Mu12_IP6_part5_v";

   string path_Mu7_IP4_part0   = "HLT_Mu7_IP4_part0_v";
   string path_Mu7_IP4_part1   = "HLT_Mu7_IP4_part1_v";
   string path_Mu7_IP4_part2   = "HLT_Mu7_IP4_part2_v";
   string path_Mu7_IP4_part3   = "HLT_Mu7_IP4_part3_v";
   string path_Mu7_IP4_part4   = "HLT_Mu7_IP4_part4_v";
   string path_Mu7_IP4_part5   = "HLT_Mu7_IP4_part5_v";
   
   string path_Mu8p5_IP3p5_part0  = "HLT_Mu8p5_IP3p5_part0_v";
   string path_Mu8p5_IP3p5_part1  = "HLT_Mu8p5_IP3p5_part1_v";
   string path_Mu8p5_IP3p5_part2  = "HLT_Mu8p5_IP3p5_part2_v";
   string path_Mu8p5_IP3p5_part3  = "HLT_Mu8p5_IP3p5_part3_v";
   string path_Mu8p5_IP3p5_part4  = "HLT_Mu8p5_IP3p5_part4_v";
   string path_Mu8p5_IP3p5_part5  = "HLT_Mu8p5_IP3p5_part5_v";

   string path_Mu10p5_IP3p5_part0  = "HLT_Mu10p5_IP3p5_part0_v";
   string path_Mu10p5_IP3p5_part1  = "HLT_Mu10p5_IP3p5_part1_v";
   string path_Mu10p5_IP3p5_part2  = "HLT_Mu10p5_IP3p5_part2_v";
   string path_Mu10p5_IP3p5_part3  = "HLT_Mu10p5_IP3p5_part3_v";
   string path_Mu10p5_IP3p5_part4  = "HLT_Mu10p5_IP3p5_part4_v";
   string path_Mu10p5_IP3p5_part5  = "HLT_Mu10p5_IP3p5_part5_v";

   int trgSize = trigNames.size();
   // cout<<trgSize<<endl;
   
   for(int i=0;i<trgSize;i++)
    {
     cout<<trigNames.triggerName(i)<<endl;
     if(trigNames.triggerName(i).find(path_Mu9_IP5_part0)  != string::npos)
       path_Mu9_IP5_part0  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu9_IP5_part1)  != string::npos)
       path_Mu9_IP5_part1  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu9_IP5_part2)  != string::npos)
       path_Mu9_IP5_part2  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu9_IP5_part3)  != string::npos)
       path_Mu9_IP5_part3  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu9_IP5_part4)  != string::npos)
       path_Mu9_IP5_part4  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu9_IP5_part5)  != string::npos)
       path_Mu9_IP5_part5  = trigNames.triggerName(i);
     
     if(trigNames.triggerName(i).find(path_Mu9_IP4_part0)  != string::npos)
       path_Mu9_IP4_part0  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu9_IP4_part1)  != string::npos)
       path_Mu9_IP4_part1  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu9_IP4_part2)  != string::npos)
       path_Mu9_IP4_part2  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu9_IP4_part3)  != string::npos)
       path_Mu9_IP4_part3  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu9_IP4_part4)  != string::npos)
       path_Mu9_IP4_part4  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu9_IP4_part5)  != string::npos)
       path_Mu9_IP4_part5  = trigNames.triggerName(i);
       
     if(trigNames.triggerName(i).find(path_Mu9_IP6_part0)  != string::npos)
       path_Mu9_IP6_part0  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu9_IP6_part1)  != string::npos)
       path_Mu9_IP6_part1  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu9_IP6_part2)  != string::npos)
       path_Mu9_IP6_part2  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu9_IP6_part3)  != string::npos)
       path_Mu9_IP6_part3  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu9_IP6_part4)  != string::npos)
       path_Mu9_IP6_part4  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu9_IP6_part5)  != string::npos)
       path_Mu9_IP6_part5  = trigNames.triggerName(i);

     if(trigNames.triggerName(i).find(path_Mu8_IP6_part0)  != string::npos)
       path_Mu8_IP6_part0  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8_IP6_part1)  != string::npos)
       path_Mu8_IP6_part1  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8_IP6_part2)  != string::npos)
       path_Mu8_IP6_part2  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8_IP6_part3)  != string::npos)
       path_Mu8_IP6_part3  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8_IP6_part4)  != string::npos)
       path_Mu8_IP6_part4  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8_IP6_part5)  != string::npos)
       path_Mu8_IP6_part5  = trigNames.triggerName(i);
       
     if(trigNames.triggerName(i).find(path_Mu8_IP5_part0)  != string::npos)
       path_Mu8_IP5_part0  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8_IP5_part1)  != string::npos)
       path_Mu8_IP5_part1  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8_IP5_part2)  != string::npos)
       path_Mu8_IP5_part2  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8_IP5_part3)  != string::npos)
       path_Mu8_IP5_part3  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8_IP5_part4)  != string::npos)
       path_Mu8_IP5_part4  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8_IP5_part5)  != string::npos)
       path_Mu8_IP5_part5  = trigNames.triggerName(i);
       
     if(trigNames.triggerName(i).find(path_Mu8_IP3_part0)  != string::npos)
       path_Mu8_IP3_part0  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8_IP3_part1)  != string::npos)
       path_Mu8_IP3_part1  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8_IP3_part2)  != string::npos)
       path_Mu8_IP3_part2  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8_IP3_part3)  != string::npos)
       path_Mu8_IP3_part3  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8_IP3_part4)  != string::npos)
       path_Mu8_IP3_part4  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8_IP3_part5)  != string::npos)
       path_Mu8_IP3_part5  = trigNames.triggerName(i);

     if(trigNames.triggerName(i).find(path_Mu12_IP6_part0)  != string::npos)
       path_Mu12_IP6_part0  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu12_IP6_part1)  != string::npos)
       path_Mu12_IP6_part1  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu12_IP6_part2)  != string::npos)
       path_Mu12_IP6_part2  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu12_IP6_part3)  != string::npos)
       path_Mu12_IP6_part3  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu12_IP6_part4)  != string::npos)
       path_Mu12_IP6_part4  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu12_IP6_part5)  != string::npos)
       path_Mu12_IP6_part5  = trigNames.triggerName(i);
       
     if(trigNames.triggerName(i).find(path_Mu7_IP4_part0)  != string::npos)
       path_Mu7_IP4_part0  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu7_IP4_part1)  != string::npos)
       path_Mu7_IP4_part1  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu7_IP4_part2)  != string::npos)
       path_Mu7_IP4_part2  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu7_IP4_part3)  != string::npos)
       path_Mu7_IP4_part3  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu7_IP4_part4)  != string::npos)
       path_Mu7_IP4_part4  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu7_IP4_part5)  != string::npos)
       path_Mu7_IP4_part5  = trigNames.triggerName(i);
       
     if(trigNames.triggerName(i).find(path_Mu10p5_IP3p5_part0)  != string::npos)
       path_Mu10p5_IP3p5_part0  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu10p5_IP3p5_part1)  != string::npos)
       path_Mu10p5_IP3p5_part1  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu10p5_IP3p5_part2)  != string::npos)
       path_Mu10p5_IP3p5_part2  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu10p5_IP3p5_part3)  != string::npos)
       path_Mu10p5_IP3p5_part3  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu10p5_IP3p5_part4)  != string::npos)
       path_Mu10p5_IP3p5_part4  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu10p5_IP3p5_part5)  != string::npos)
       path_Mu10p5_IP3p5_part5  = trigNames.triggerName(i);
       
     if(trigNames.triggerName(i).find(path_Mu8p5_IP3p5_part0)  != string::npos)
       path_Mu8p5_IP3p5_part0  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8p5_IP3p5_part1)  != string::npos)
       path_Mu8p5_IP3p5_part1  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8p5_IP3p5_part2)  != string::npos)
       path_Mu8p5_IP3p5_part2  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8p5_IP3p5_part3)  != string::npos)
       path_Mu8p5_IP3p5_part3  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8p5_IP3p5_part4)  != string::npos)
       path_Mu8p5_IP3p5_part4  = trigNames.triggerName(i);
     if(trigNames.triggerName(i).find(path_Mu8p5_IP3p5_part5)  != string::npos)
       path_Mu8p5_IP3p5_part5  = trigNames.triggerName(i);
       
    }
   bool pass_Mu9_IP5_part0 = false,pass_Mu9_IP5_part1 = false,pass_Mu9_IP5_part2 = false,pass_Mu9_IP5_part3 = false,pass_Mu9_IP5_part4 = false,pass_Mu9_IP5_part5 = false;
   bool pass_Mu9_IP4_part0 = false,pass_Mu9_IP4_part1 = false,pass_Mu9_IP4_part2 = false,pass_Mu9_IP4_part3 = false,pass_Mu9_IP4_part4 = false,pass_Mu9_IP4_part5 = false;
   bool pass_Mu9_IP6_part0 = false,pass_Mu9_IP6_part1 = false,pass_Mu9_IP6_part2 = false,pass_Mu9_IP6_part3 = false,pass_Mu9_IP6_part4 = false,pass_Mu9_IP6_part5 = false;
   bool pass_Mu8_IP6_part0 = false,pass_Mu8_IP6_part1 = false,pass_Mu8_IP6_part2 = false,pass_Mu8_IP6_part3 = false,pass_Mu8_IP6_part4 = false,pass_Mu8_IP6_part5 = false;
   bool pass_Mu8_IP5_part0 = false,pass_Mu8_IP5_part1 = false,pass_Mu8_IP5_part2 = false,pass_Mu8_IP5_part3 = false,pass_Mu8_IP5_part4 = false,pass_Mu8_IP5_part5 = false;
   bool pass_Mu8_IP3_part0 = false,pass_Mu8_IP3_part1 = false,pass_Mu8_IP3_part2 = false,pass_Mu8_IP3_part3 = false,pass_Mu8_IP3_part4 = false,pass_Mu8_IP3_part5 = false;
   bool pass_Mu12_IP6_part0 = false,pass_Mu12_IP6_part1 = false,pass_Mu12_IP6_part2 = false,pass_Mu12_IP6_part3 = false,pass_Mu12_IP6_part4 = false,pass_Mu12_IP6_part5 = false;
   bool pass_Mu7_IP4_part0 = false,pass_Mu7_IP4_part1 = false,pass_Mu7_IP4_part2 = false,pass_Mu7_IP4_part3 = false,pass_Mu7_IP4_part4 = false,pass_Mu7_IP4_part5 = false;
   bool pass_Mu8p5_IP3p5_part0 = false,pass_Mu8p5_IP3p5_part1 = false,pass_Mu8p5_IP3p5_part2 = false,pass_Mu8p5_IP3p5_part3 = false,pass_Mu8p5_IP3p5_part4 = false,pass_Mu8p5_IP3p5_part5 = false;
   bool pass_Mu10p5_IP3p5_part0 = false,pass_Mu10p5_IP3p5_part1 = false,pass_Mu10p5_IP3p5_part2 = false,pass_Mu10p5_IP3p5_part3 = false,pass_Mu10p5_IP3p5_part4 = false,pass_Mu10p5_IP3p5_part5 = false;
   
   if (trigNames.triggerIndex(path_Mu9_IP5_part0) != (unsigned)trgSize)  pass_Mu9_IP5_part0  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP5_part0));
   if (trigNames.triggerIndex(path_Mu9_IP5_part1) != (unsigned)trgSize)  pass_Mu9_IP5_part1  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP5_part1));
   if (trigNames.triggerIndex(path_Mu9_IP5_part2) != (unsigned)trgSize)  pass_Mu9_IP5_part2  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP5_part2));
   if (trigNames.triggerIndex(path_Mu9_IP5_part3) != (unsigned)trgSize)  pass_Mu9_IP5_part3  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP5_part3));
   if (trigNames.triggerIndex(path_Mu9_IP5_part4) != (unsigned)trgSize)  pass_Mu9_IP5_part4  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP5_part4));
   if (trigNames.triggerIndex(path_Mu9_IP5_part5) != (unsigned)trgSize)  pass_Mu9_IP5_part5  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP5_part5));

   if (trigNames.triggerIndex(path_Mu9_IP4_part0) != (unsigned)trgSize)  pass_Mu9_IP4_part0  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP4_part0));
   if (trigNames.triggerIndex(path_Mu9_IP4_part1) != (unsigned)trgSize)  pass_Mu9_IP4_part1  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP4_part1));
   if (trigNames.triggerIndex(path_Mu9_IP4_part2) != (unsigned)trgSize)  pass_Mu9_IP4_part2  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP4_part2));
   if (trigNames.triggerIndex(path_Mu9_IP4_part3) != (unsigned)trgSize)  pass_Mu9_IP4_part3  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP4_part3));
   if (trigNames.triggerIndex(path_Mu9_IP4_part4) != (unsigned)trgSize)  pass_Mu9_IP4_part4  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP4_part4));
   if (trigNames.triggerIndex(path_Mu9_IP4_part5) != (unsigned)trgSize)  pass_Mu9_IP4_part5  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP4_part5));

   if (trigNames.triggerIndex(path_Mu9_IP6_part0) != (unsigned)trgSize)  pass_Mu9_IP6_part0  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP6_part0));
   if (trigNames.triggerIndex(path_Mu9_IP6_part1) != (unsigned)trgSize)  pass_Mu9_IP6_part1  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP6_part1));
   if (trigNames.triggerIndex(path_Mu9_IP6_part2) != (unsigned)trgSize)  pass_Mu9_IP6_part2  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP6_part2));
   if (trigNames.triggerIndex(path_Mu9_IP6_part3) != (unsigned)trgSize)  pass_Mu9_IP6_part3  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP6_part3));
   if (trigNames.triggerIndex(path_Mu9_IP6_part4) != (unsigned)trgSize)  pass_Mu9_IP6_part4  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP6_part4));
   if (trigNames.triggerIndex(path_Mu9_IP6_part5) != (unsigned)trgSize)  pass_Mu9_IP6_part5  = triggerBits->accept(trigNames.triggerIndex(path_Mu9_IP6_part5));

   if (trigNames.triggerIndex(path_Mu8_IP6_part0) != (unsigned)trgSize)  pass_Mu8_IP6_part0  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP6_part0));
   if (trigNames.triggerIndex(path_Mu8_IP6_part1) != (unsigned)trgSize)  pass_Mu8_IP6_part1  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP6_part1));
   if (trigNames.triggerIndex(path_Mu8_IP6_part2) != (unsigned)trgSize)  pass_Mu8_IP6_part2  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP6_part2));
   if (trigNames.triggerIndex(path_Mu8_IP6_part3) != (unsigned)trgSize)  pass_Mu8_IP6_part3  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP6_part3));
   if (trigNames.triggerIndex(path_Mu8_IP6_part4) != (unsigned)trgSize)  pass_Mu8_IP6_part4  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP6_part4));
   if (trigNames.triggerIndex(path_Mu8_IP6_part5) != (unsigned)trgSize)  pass_Mu8_IP6_part5  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP6_part5));

   if (trigNames.triggerIndex(path_Mu8_IP5_part0) != (unsigned)trgSize)  pass_Mu8_IP5_part0  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP5_part0));
   if (trigNames.triggerIndex(path_Mu8_IP5_part1) != (unsigned)trgSize)  pass_Mu8_IP5_part1  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP5_part1));
   if (trigNames.triggerIndex(path_Mu8_IP5_part2) != (unsigned)trgSize)  pass_Mu8_IP5_part2  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP5_part2));
   if (trigNames.triggerIndex(path_Mu8_IP5_part3) != (unsigned)trgSize)  pass_Mu8_IP5_part3  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP5_part3));
   if (trigNames.triggerIndex(path_Mu8_IP5_part4) != (unsigned)trgSize)  pass_Mu8_IP5_part4  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP5_part4));
   if (trigNames.triggerIndex(path_Mu8_IP5_part5) != (unsigned)trgSize)  pass_Mu8_IP5_part5  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP5_part5));

   if (trigNames.triggerIndex(path_Mu8_IP3_part0) != (unsigned)trgSize)  pass_Mu8_IP3_part0  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP3_part0));
   if (trigNames.triggerIndex(path_Mu8_IP3_part1) != (unsigned)trgSize)  pass_Mu8_IP3_part1  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP3_part1));
   if (trigNames.triggerIndex(path_Mu8_IP3_part2) != (unsigned)trgSize)  pass_Mu8_IP3_part2  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP3_part2));
   if (trigNames.triggerIndex(path_Mu8_IP3_part3) != (unsigned)trgSize)  pass_Mu8_IP3_part3  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP3_part3));
   if (trigNames.triggerIndex(path_Mu8_IP3_part4) != (unsigned)trgSize)  pass_Mu8_IP3_part4  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP3_part4));
   if (trigNames.triggerIndex(path_Mu8_IP3_part5) != (unsigned)trgSize)  pass_Mu8_IP3_part5  = triggerBits->accept(trigNames.triggerIndex(path_Mu8_IP3_part5));

   if (trigNames.triggerIndex(path_Mu12_IP6_part0) != (unsigned)trgSize)  pass_Mu12_IP6_part0  = triggerBits->accept(trigNames.triggerIndex(path_Mu12_IP6_part0));
   if (trigNames.triggerIndex(path_Mu12_IP6_part1) != (unsigned)trgSize)  pass_Mu12_IP6_part1  = triggerBits->accept(trigNames.triggerIndex(path_Mu12_IP6_part1));
   if (trigNames.triggerIndex(path_Mu12_IP6_part2) != (unsigned)trgSize)  pass_Mu12_IP6_part2  = triggerBits->accept(trigNames.triggerIndex(path_Mu12_IP6_part2));
   if (trigNames.triggerIndex(path_Mu12_IP6_part3) != (unsigned)trgSize)  pass_Mu12_IP6_part3  = triggerBits->accept(trigNames.triggerIndex(path_Mu12_IP6_part3));
   if (trigNames.triggerIndex(path_Mu12_IP6_part4) != (unsigned)trgSize)  pass_Mu12_IP6_part4  = triggerBits->accept(trigNames.triggerIndex(path_Mu12_IP6_part4));
   if (trigNames.triggerIndex(path_Mu12_IP6_part5) != (unsigned)trgSize)  pass_Mu12_IP6_part5  = triggerBits->accept(trigNames.triggerIndex(path_Mu12_IP6_part5));

   if (trigNames.triggerIndex(path_Mu7_IP4_part0) != (unsigned)trgSize)  pass_Mu7_IP4_part0  = triggerBits->accept(trigNames.triggerIndex(path_Mu7_IP4_part0));
   if (trigNames.triggerIndex(path_Mu7_IP4_part1) != (unsigned)trgSize)  pass_Mu7_IP4_part1  = triggerBits->accept(trigNames.triggerIndex(path_Mu7_IP4_part1));
   if (trigNames.triggerIndex(path_Mu7_IP4_part2) != (unsigned)trgSize)  pass_Mu7_IP4_part2  = triggerBits->accept(trigNames.triggerIndex(path_Mu7_IP4_part2));
   if (trigNames.triggerIndex(path_Mu7_IP4_part3) != (unsigned)trgSize)  pass_Mu7_IP4_part3  = triggerBits->accept(trigNames.triggerIndex(path_Mu7_IP4_part3));
   if (trigNames.triggerIndex(path_Mu7_IP4_part4) != (unsigned)trgSize)  pass_Mu7_IP4_part4  = triggerBits->accept(trigNames.triggerIndex(path_Mu7_IP4_part4));
   if (trigNames.triggerIndex(path_Mu7_IP4_part5) != (unsigned)trgSize)  pass_Mu7_IP4_part5  = triggerBits->accept(trigNames.triggerIndex(path_Mu7_IP4_part5));

   if (trigNames.triggerIndex(path_Mu10p5_IP3p5_part0) != (unsigned)trgSize)  pass_Mu10p5_IP3p5_part0  = triggerBits->accept(trigNames.triggerIndex(path_Mu10p5_IP3p5_part0));
   if (trigNames.triggerIndex(path_Mu10p5_IP3p5_part1) != (unsigned)trgSize)  pass_Mu10p5_IP3p5_part1  = triggerBits->accept(trigNames.triggerIndex(path_Mu10p5_IP3p5_part1));
   if (trigNames.triggerIndex(path_Mu10p5_IP3p5_part2) != (unsigned)trgSize)  pass_Mu10p5_IP3p5_part2  = triggerBits->accept(trigNames.triggerIndex(path_Mu10p5_IP3p5_part2));
   if (trigNames.triggerIndex(path_Mu10p5_IP3p5_part3) != (unsigned)trgSize)  pass_Mu10p5_IP3p5_part3  = triggerBits->accept(trigNames.triggerIndex(path_Mu10p5_IP3p5_part3));
   if (trigNames.triggerIndex(path_Mu10p5_IP3p5_part4) != (unsigned)trgSize)  pass_Mu10p5_IP3p5_part4  = triggerBits->accept(trigNames.triggerIndex(path_Mu10p5_IP3p5_part4));
   if (trigNames.triggerIndex(path_Mu10p5_IP3p5_part5) != (unsigned)trgSize)  pass_Mu10p5_IP3p5_part5  = triggerBits->accept(trigNames.triggerIndex(path_Mu10p5_IP3p5_part5));

   if (trigNames.triggerIndex(path_Mu8p5_IP3p5_part0) != (unsigned)trgSize)  pass_Mu8p5_IP3p5_part0  = triggerBits->accept(trigNames.triggerIndex(path_Mu8p5_IP3p5_part0));
   if (trigNames.triggerIndex(path_Mu8p5_IP3p5_part1) != (unsigned)trgSize)  pass_Mu8p5_IP3p5_part1  = triggerBits->accept(trigNames.triggerIndex(path_Mu8p5_IP3p5_part1));
   if (trigNames.triggerIndex(path_Mu8p5_IP3p5_part2) != (unsigned)trgSize)  pass_Mu8p5_IP3p5_part2  = triggerBits->accept(trigNames.triggerIndex(path_Mu8p5_IP3p5_part2));
   if (trigNames.triggerIndex(path_Mu8p5_IP3p5_part3) != (unsigned)trgSize)  pass_Mu8p5_IP3p5_part3  = triggerBits->accept(trigNames.triggerIndex(path_Mu8p5_IP3p5_part3));
   if (trigNames.triggerIndex(path_Mu8p5_IP3p5_part4) != (unsigned)trgSize)  pass_Mu8p5_IP3p5_part4  = triggerBits->accept(trigNames.triggerIndex(path_Mu8p5_IP3p5_part4));
   if (trigNames.triggerIndex(path_Mu8p5_IP3p5_part5) != (unsigned)trgSize)  pass_Mu8p5_IP3p5_part5  = triggerBits->accept(trigNames.triggerIndex(path_Mu8p5_IP3p5_part5));

// Trigger-objects ************************************************************
   Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
   iEvent.getByToken(triggerObjects_, triggerObjects);

   Handle<pat::PackedTriggerPrescales> triggerPrescales;
   iEvent.getByToken(triggerPrescales_, triggerPrescales);

  // edm::Handle<reco::VertexCollection> vertices;
  Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  // Handle<std::vector<reco::VertexCompositePtrCandidate>> secondary_vertices;
  // iEvent.getByToken(secondary_vtxToken_, secondary_vertices);
  
  // edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> secondary_vertices;
  // iEvent.getByToken(secondary_vtxToken_, secondary_vertices);
  
  // edm::Handle< std::vector<reco::VertexCompositePtrCandidate> > secondary_vertices;
  // iEvent.getByLabel("slimmedSecondaryVertices", secondary_vertices);
  
  edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> secondaryVerticesHandle;
  iEvent.getByToken(secondaryVerticesToken_, secondaryVerticesHandle);
  
  // Handle<GenParticleCollection> genParticles;
  // iEvent.getByLabel("genParticles", genParticles);
  
// Trigger Object Collections
  double L3_Mu_pt[4];
  double L3_Mu_eta[4];
  double L3_Mu_phi[4];
  double L3_Mu_filters[4][10];
  int N_L3_Muons = 0;
  // cout<<"Before the trigger object collection"<<endl;
  cout<<endl;
  for(int i=0;i<4;i++)
    L3_Mu_pt[i]=L3_Mu_eta[i]=L3_Mu_phi[i]=-100;

    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    for(pat::TriggerObjectStandAlone obj : *triggerObjects)
       {
        obj.unpackPathNames(names);
        obj.unpackFilterLabels(iEvent,*triggerBits);

        //if(obj.collection()=="hltIterL3MuonCandidates::HLT")
        //  {
           bool condition = false;
           for(unsigned h=0; h<obj.filterLabels().size(); ++h) {
            // if (obj.filterLabels()[h].find("10p5") != std::string::npos) {
            //    std::cout << obj.filterLabels()[h] << '\n';
            // }
            // if (pass_Mu8p5_IP3p5_part0 || pass_Mu8p5_IP3p5_part1 || pass_Mu8p5_IP3p5_part2 || pass_Mu8p5_IP3p5_part3 || pass_Mu8p5_IP3p5_part4) {
            std::cout << obj.filterLabels()[h] << '\n';
            // }
            if(obj.filterLabels()[h]=="hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered12Q"){  // HLT_Mu12_IP6_part*
            L3_Mu_filters[N_L3_Muons][0]=1; condition = true;}
            if(obj.filterLabels()[h]=="hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered10p5Q"){  // HLT_Mu10p5_IP3p5 //_part*
            L3_Mu_filters[N_L3_Muons][1]=1; condition = true;}
            if(obj.filterLabels()[h]=="hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9Q") {   // HLT_Mu9_IP6_part*
            L3_Mu_filters[N_L3_Muons][2]=1; condition = true;}
            if(obj.filterLabels()[h]=="hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP5Q"){ // HLT_Mu9_IP5_part*
            L3_Mu_filters[N_L3_Muons][3]=1; condition = true;}
            if(obj.filterLabels()[h]=="hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP4Q"){ // HLT_Mu9_IP4_part*
            L3_Mu_filters[N_L3_Muons][4]=1; condition = true;}
            if(obj.filterLabels()[h]=="hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8p5Q") {// HLT_Mu8p5_IP3p5 // _part*
            L3_Mu_filters[N_L3_Muons][5]=1; condition = true;}
            if(obj.filterLabels()[h]=="hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8IP6Q") {// HLT_Mu8_IP6_part*
            L3_Mu_filters[N_L3_Muons][6]=1; condition = true;}
            if(obj.filterLabels()[h]=="hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8IP5Q") {// HLT_Mu8_IP5_part*
            L3_Mu_filters[N_L3_Muons][7]=1; condition = true;}
            if(obj.filterLabels()[h]=="hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8Q")   { // HLT_Mu8_IP3_part*
            L3_Mu_filters[N_L3_Muons][8]=1; condition = true;}
            if(obj.filterLabels()[h]=="hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered7IP4Q") {// HLT_Mu7_IP4_part*
            L3_Mu_filters[N_L3_Muons][9]=1; condition = true;}
            }
            if (condition) {
                L3_Mu_pt[N_L3_Muons]   = obj.pt();
                L3_Mu_eta[N_L3_Muons]  = obj.eta();
                L3_Mu_phi[N_L3_Muons]  = obj.phi();
                N_L3_Muons++;
            }
          // }
        }
   cout<<"After the trigger object collection"<<endl;
// Primary Vertex *************************************************************
  // bool norm_VTX = false;
  // double PV_x=5000;
  // double PV_y=5000;
  // double PV_z=5000;
  cout<<1<<endl;

  const reco::Vertex &PV = vertices->front();
  const reco::VertexCompositePtrCandidate &SV = secondaryVerticesHandle->front();
  double SV_probability = 10;
  // if (SV.isValid()) {
  cout<<2<<endl;
  if (secondaryVerticesHandle.isValid() && !(secondaryVerticesHandle->empty())) {
     cout<<SV.vertex().x()<<endl;
     cout<<SV.vertexChi2()<<endl;
     SV_probability = TMath::Prob(SV.vertexChi2(), SV.vertexNdof());
  }
  // }
  cout<<3<<endl;

  /*
  for (const reco::Vertex &vtx : *vertices) {
     // for (const reco::Vertex &snd_vtx : *secondary_vertices) {
        if( vtx.ndof()>4  &&  fabs(vtx.position().z())<=24  &&  fabs(vtx.position().rho()) <= 2)
        {
          // PV_x = vtx.position().x();
          // PV_y = vtx.position().y();
          // PV_z = vtx.position().z();
          norm_VTX=true;
          break;
        }
     // }
  }
  // if(!norm_VTX) return;
  */

// Offline Muons ***************************************************************

  double EtaMax=2.3;
  double MinPt=1.0;
  int    N_Off_Mu=0;

  double mupt[2];
  double mue[2];
  double mueta[2];
  double muphi[2];
  double mucharge[2];
  double mudxy[2];
  double mudxyErr[2];
  double mudz[2];
  bool muLooseID[2];
  bool muTightID[2];
  bool mu_isGlobal[2];
  int n_strip_layers[2];
  int n_tracker_layers[2];
  int n_pixel_layers[2];
  float mu_PV_theta[2];
  bool track_purity[2];
  float probSV[2];

  for(int i=0;i<2;i++){mupt[i]=mue[i]=mueta[i]=muphi[i]=mucharge[i]=mudxy[i]=mudxyErr[i]=mudz[i]=-100;
  muLooseID[i]=muTightID[i]=mu_isGlobal[i]=track_purity[i]=false;n_strip_layers[i]=n_tracker_layers[i]=n_pixel_layers[i]=-100;mu_PV_theta[i]=probSV[i]=-100.0;}
  
  edm::Handle<edm::View<pat::Muon>> recmuons;
  iEvent.getByToken(recmuonToken_, recmuons);
  Bool_t tight_id = false;
  Int_t tight_index = -10;
  Bool_t loose_id = false;
  Bool_t other = false;
  Int_t index = 0;
  cout<<"before the muon loop"<<endl;
  for (auto mu = recmuons->begin(); mu != recmuons->end(); mu++)
   {
    cout<<mu->pt()<<"    "<<fabs(mu->eta())<<"     "<<mu->isTightMuon(PV)<<endl;
    if(mu->pt() < MinPt || fabs(mu->eta()) > EtaMax) continue;
    if (mu->isTightMuon(PV) && !tight_id) { // The tag muon wasn't selected and is here.
        cout<<"Tight id"<<endl;
        tight_id = true;
        tight_index = index;
        N_Off_Mu=0;
        mupt[N_Off_Mu]=mu->pt();
        mueta[N_Off_Mu]=mu->eta();
        muphi[N_Off_Mu]=mu->phi();
        mue[N_Off_Mu]=mu->energy();
        mucharge[N_Off_Mu]=mu->charge();
        mudxy[N_Off_Mu]=mu->muonBestTrack()->dxy(PV.position());
        mudxyErr[N_Off_Mu]=mu->muonBestTrack()->d0Error();
        mudz[N_Off_Mu] = mu->muonBestTrack()->dz(PV.position());
        muLooseID[N_Off_Mu] = mu->isLooseMuon();
        muTightID[N_Off_Mu] = mu->isTightMuon(PV);
        mu_isGlobal[N_Off_Mu] = mu->isGlobalMuon();
        if (mu->isGlobalMuon() && mu->muonBestTrack().isNonnull()) {
            n_tracker_layers[N_Off_Mu] = mu->muonBestTrack()->hitPattern().trackerLayersWithMeasurement();
            n_strip_layers[N_Off_Mu] = mu->muonBestTrack()->hitPattern().numberOfValidStripHits();
            n_pixel_layers[N_Off_Mu] = mu->muonBestTrack()->hitPattern().pixelLayersWithMeasurement();
            track_purity[N_Off_Mu] = mu->muonBestTrack()->quality(reco::TrackBase::highPurity);
            probSV[N_Off_Mu] = TrigScaleFactors_OnBParking::calculateMuonSVProbability(*mu, *vertices);
            if (secondaryVerticesHandle.isValid() && !(secondaryVerticesHandle->empty())) {
                mu_PV_theta[N_Off_Mu] = TrigScaleFactors_OnBParking::calculateBackpointingAngle(PV, SV, *mu);
                cout<<mu_PV_theta[N_Off_Mu]<<endl;
            } else {
                mu_PV_theta[N_Off_Mu] = -10.0;
            }
        } else {
            n_tracker_layers[N_Off_Mu] = -1;
            n_strip_layers[N_Off_Mu] = -1;
            n_pixel_layers[N_Off_Mu] = -1;
            track_purity[N_Off_Mu] = false;
            probSV[N_Off_Mu] = -1;
            mu_PV_theta[N_Off_Mu] = -1.0;
        }
        break;
     }
     index++;
   }
   cout<<"before the muon loop"<<endl;
   index = 0;
   for (auto mu = recmuons->begin(); mu != recmuons->end(); mu++)
   {
    if(mu->pt() < MinPt || fabs(mu->eta()) > EtaMax) continue;
     if (mu->isLooseMuon() && !loose_id && index != tight_index) { // The tag muon wasn't selected and is here.
        cout<<"Loose id"<<endl;
        loose_id = true;
        N_Off_Mu=1;
        mupt[N_Off_Mu]=mu->pt();
        mueta[N_Off_Mu]=mu->eta();
        muphi[N_Off_Mu]=mu->phi();
        mue[N_Off_Mu]=mu->energy();
        mucharge[N_Off_Mu]=mu->charge();
        mudxy[N_Off_Mu]=mu->muonBestTrack()->dxy(PV.position());
        mudxyErr[N_Off_Mu]=mu->muonBestTrack()->d0Error();
        mudz[N_Off_Mu] = mu->muonBestTrack()->dz(PV.position());
        muLooseID[N_Off_Mu] = mu->isLooseMuon();
        muTightID[N_Off_Mu] = mu->isTightMuon(PV);
        mu_isGlobal[N_Off_Mu] = mu->isGlobalMuon();
        if (mu->isGlobalMuon() && mu->muonBestTrack().isNonnull()) {
            n_tracker_layers[N_Off_Mu] = mu->muonBestTrack()->hitPattern().trackerLayersWithMeasurement();
            n_strip_layers[N_Off_Mu] = mu->muonBestTrack()->hitPattern().numberOfValidStripHits();
            n_pixel_layers[N_Off_Mu] = mu->muonBestTrack()->hitPattern().pixelLayersWithMeasurement();
            track_purity[N_Off_Mu] = mu->muonBestTrack()->quality(reco::TrackBase::highPurity);
            probSV[N_Off_Mu] = TrigScaleFactors_OnBParking::calculateMuonSVProbability(*mu, *vertices);
            if (secondaryVerticesHandle.isValid() && !(secondaryVerticesHandle->empty())) {
                mu_PV_theta[N_Off_Mu] = TrigScaleFactors_OnBParking::calculateBackpointingAngle(PV, SV, *mu);
                cout<<mu_PV_theta[N_Off_Mu]<<endl;
            } else {
                mu_PV_theta[N_Off_Mu] = -10.0;
            }
        } else {
            n_tracker_layers[N_Off_Mu] = -1;
            n_strip_layers[N_Off_Mu] = -1;
            n_pixel_layers[N_Off_Mu] = -1;
            track_purity[N_Off_Mu] = false;
            probSV[N_Off_Mu] = -1;
            mu_PV_theta[N_Off_Mu] = -1.0;
        }
        break;
     }
     index++;
   }
   cout<<"the muons had been selected"<<endl;
   
   double dr_tag = 10;
   double dr_probe = 10;
   /*
   for(size_t i = 0; i < genParticles->size(); ++ i) {
     const GenParticle & p = (*genParticles)[i];
     int id = p.pdgId();
     int st = p.status();
     const Candidate * mom = p.mother();
     double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass = p.mass();
     double vx = p.vx(), vy = p.vy(), vz = p.vz();
     int charge = p.charge();
     if (abs(id) == 13 && mom->pdgId() == 443) {
       double dr_tag_0 = TrigScaleFactors_OnBParking::DeltaR(mueta[0], eta, muphi[0], phi);
       double dr_probe_0 = TrigScaleFactors_OnBParking::DeltaR(mueta[1], eta, muphi[1], phi);
       if (dr_tag > dr_tag_0) {
           dr_tag = dr_tag_0;
       }
       if (dr_probe > dr_probe_0) {
           dr_probe = dr_probe_0;
       }
     }
   }*/
   double dr_tag_track = 10;
   double dr_probe_track = 10;
   /*
   for (int t = 0; t < )
         if (mu->isGlobalMuon() && mu->innerTrack().isNonnull()) {
            // n_tracker_layers[N_Off_Mu] = mu->innerTrack()->hitPattern().numberOfValidTrackerHits();
            // n_strip_layers[N_Off_Mu] = mu->innerTrack()->hitPattern().numberOfValidStripHits();
            // n_pixel_layers[N_Off_Mu] = mu->innerTrack()->hitPattern().numberOfValidPixelHits();
            n_tracker_layers[N_Off_Mu] = mu->innerTrack()->hitPattern().trackerLayersWithMeasurement();
            n_strip_layers[N_Off_Mu] = mu->innerTrack()->hitPattern().numberOfValidStripHits();
            n_pixel_layers[N_Off_Mu] = mu->innerTrack()->hitPattern().pixelLayersWithMeasurement();
            // track_purity[N_Off_Mu] = mu->innerTrack()->hitPattern().trackHighPurity(); //->quality(reco::TrackBase::highPurity);
            track_purity[N_Off_Mu] = mu->innerTrack()->hitPattern()->quality("trackHighPurity");
            probSV[N_Off_Mu] = TrigScaleFactors_OnBParking::calculateMuonSVProbability(*mu, *vertices);
            if (secondaryVerticesHandle.isValid() && !(secondaryVerticesHandle->empty())) {
                mu_PV_theta[N_Off_Mu] = TrigScaleFactors_OnBParking::calculateBackpointingAngle(PV, SV, *mu);
                cout<<mu_PV_theta[N_Off_Mu]<<endl;
            } else {
                mu_PV_theta[N_Off_Mu] = -10.0;
            }
        } else {
            n_tracker_layers[N_Off_Mu] = -1;
            n_strip_layers[N_Off_Mu] = -1;
            n_pixel_layers[N_Off_Mu] = -1;
            track_purity[N_Off_Mu] = false;
            probSV[N_Off_Mu] = -1;
            mu_PV_theta[N_Off_Mu] = -1.0;
        }
   */
   // Matching with the tracks from the Secondary Vertex
   // auto tracks = SV.product();
   dr_tag = 10;
   dr_probe = 10;

   // for (auto itTrack = SV.tracks_begin(); itTrack != SV.tracks_end(); ++itTrack) {
   if (secondaryVerticesHandle.isValid() && !(secondaryVerticesHandle->empty())) {

   for (size_t d = 0; d < SV.numberOfDaughters(); d++) {
       const auto & track = SV.daughterPtr(d);
       double dr_tag_0 = TrigScaleFactors_OnBParking::DeltaR(mueta[0], track->eta(), muphi[0], track->phi());
       double dr_probe_0 = TrigScaleFactors_OnBParking::DeltaR(mueta[1], track->eta(), muphi[1], track->phi());
       if (dr_tag > dr_tag_0) {
           dr_tag = dr_tag_0;
       }
       if (dr_probe > dr_probe_0) {
           dr_probe = dr_probe_0;
       }
   }

   }
// ************* Selection

   // if(N_L3_Muons < 1) return; // || N_Off_Mu < 2
   // if(phi_dist(muphi[0],muphi[1]) > 1.5) return;

// ************* PU info for reweight
    /*
   edm::Handle<std::vector<PileupSummaryInfo>>  PupInfo;
   iEvent.getByToken(puToken_, PupInfo);
   std::vector<PileupSummaryInfo>::const_iterator PVI;
   int Tnpv = -1;
   for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) 
   {
      int BX = PVI->getBunchCrossing();
      if(BX == 0) {
         Tnpv = PVI->getTrueNumInteractions();
         continue;
      }
   }
   */
// *************************************** Ntuples ******************************************
         _Tag_Candidate   -> Fill (mupt[0],mue[0],mueta[0],muphi[0],mucharge[0],mudxy[0],mudxyErr[0],mudz[0],muLooseID[0],muTightID[0],mu_isGlobal[0],1.0,-10);
         _Tag_Candidate_Track -> Fill(n_tracker_layers[0],n_strip_layers[0],n_pixel_layers[0],mu_PV_theta[0],track_purity[0],probSV[0],dr_tag);
         _Probe_Candidate -> Fill (mupt[1],mue[1],mueta[1],muphi[1],mucharge[1],mudxy[1],mudxyErr[1],mudz[1],muLooseID[1],muTightID[1],mu_isGlobal[1]);
         _Probe_Candidate_Track -> Fill(n_tracker_layers[1],n_strip_layers[1],n_pixel_layers[1],mu_PV_theta[1],track_purity[1],probSV[1],dr_probe);

         _L3_0 -> Fill (L3_Mu_pt[0],L3_Mu_eta[0],L3_Mu_phi[0],L3_Mu_filters[0][0],L3_Mu_filters[0][1],L3_Mu_filters[0][2],L3_Mu_filters[0][3],L3_Mu_filters[0][4],L3_Mu_filters[0][5],L3_Mu_filters[0][6],L3_Mu_filters[0][7],L3_Mu_filters[0][8],L3_Mu_filters[0][9]);
         _L3_1 -> Fill (L3_Mu_pt[1],L3_Mu_eta[1],L3_Mu_phi[1],L3_Mu_filters[1][0],L3_Mu_filters[1][1],L3_Mu_filters[1][2],L3_Mu_filters[1][3],L3_Mu_filters[1][4],L3_Mu_filters[1][5],L3_Mu_filters[1][6],L3_Mu_filters[1][7],L3_Mu_filters[1][8],L3_Mu_filters[1][9]);
         _L3_2 -> Fill (L3_Mu_pt[2],L3_Mu_eta[2],L3_Mu_phi[2],L3_Mu_filters[2][0],L3_Mu_filters[2][1],L3_Mu_filters[2][2],L3_Mu_filters[2][3],L3_Mu_filters[2][4],L3_Mu_filters[2][5],L3_Mu_filters[2][6],L3_Mu_filters[2][7],L3_Mu_filters[2][8],L3_Mu_filters[2][9]);
         _L3_3 -> Fill (L3_Mu_pt[3],L3_Mu_eta[3],L3_Mu_phi[3],L3_Mu_filters[3][0],L3_Mu_filters[3][1],L3_Mu_filters[3][2],L3_Mu_filters[3][3],L3_Mu_filters[3][4],L3_Mu_filters[3][5],L3_Mu_filters[3][6],L3_Mu_filters[3][7],L3_Mu_filters[3][8],L3_Mu_filters[3][9]);

         _HLT_result -> Fill (SV_probability,(pass_Mu12_IP6_part0||pass_Mu12_IP6_part1||pass_Mu12_IP6_part2||pass_Mu12_IP6_part3||pass_Mu12_IP6_part4||pass_Mu12_IP6_part5), (pass_Mu10p5_IP3p5_part0||pass_Mu10p5_IP3p5_part1||pass_Mu10p5_IP3p5_part2||pass_Mu10p5_IP3p5_part3||pass_Mu10p5_IP3p5_part4||pass_Mu10p5_IP3p5_part5), (pass_Mu9_IP6_part0||pass_Mu9_IP6_part1||pass_Mu9_IP6_part2||pass_Mu9_IP6_part3||pass_Mu9_IP6_part4||pass_Mu9_IP6_part5), (pass_Mu9_IP5_part0||pass_Mu9_IP5_part1||pass_Mu9_IP5_part2||pass_Mu9_IP5_part3||pass_Mu9_IP5_part4||pass_Mu9_IP5_part5), (pass_Mu9_IP4_part0||pass_Mu9_IP4_part1||pass_Mu9_IP4_part2||pass_Mu9_IP4_part3||pass_Mu9_IP4_part4||pass_Mu9_IP4_part5), (pass_Mu8p5_IP3p5_part0||pass_Mu8p5_IP3p5_part1||pass_Mu8p5_IP3p5_part2||pass_Mu8p5_IP3p5_part3||pass_Mu8p5_IP3p5_part4||pass_Mu8p5_IP3p5_part5), (pass_Mu8_IP5_part0||pass_Mu8_IP5_part1||pass_Mu8_IP5_part2||pass_Mu8_IP5_part3||pass_Mu8_IP5_part4||pass_Mu8_IP5_part5), (pass_Mu8_IP6_part0||pass_Mu8_IP6_part1||pass_Mu8_IP6_part2||pass_Mu8_IP6_part3||pass_Mu8_IP6_part4||pass_Mu8_IP6_part5), (pass_Mu8_IP3_part0||pass_Mu8_IP3_part1||pass_Mu8_IP3_part2||pass_Mu8_IP3_part3||pass_Mu8_IP3_part4||pass_Mu8_IP3_part5), (pass_Mu7_IP4_part0||pass_Mu7_IP4_part1||pass_Mu7_IP4_part2||pass_Mu7_IP4_part3||pass_Mu7_IP4_part4||pass_Mu7_IP4_part5));
// *************************************************************************************************

}
// ------------ method called once each job just before starting event loop  ------------
void TrigScaleFactors_OnBParking::beginJob(){
}

    TH1F *TrigScaleFactors_OnBParking::book1DHistogram(TFileDirectory & fDir, const std::string & fName, const std::string & fTitle,
    int fNbins, double fXmin, double fXmax) const {
    char title[1024];

    sprintf(title, "%s [RUN:%ld/%ld]", fTitle.c_str(), runBegin, lumibegin);
    return fDir.make < TH1F > (fName.c_str(), title, fNbins, fXmin, fXmax);
}
    TProfile *TrigScaleFactors_OnBParking::bookProfileHistogram(TFileDirectory & fDir, const std::string & fName, const std::string & fTitle,
    int fNbins, double fXmin, double fXmax, double fYmin, double fYmax) const {
    char title[1024];
    sprintf(title, "%s [RUN:%ld/%ld]", fTitle.c_str(), runBegin, lumibegin);
    return fDir.make < TProfile > (fName.c_str(), title, fNbins, fXmin, fXmax, fYmin, fYmax); 
}

    TH2F *TrigScaleFactors_OnBParking::book2DHistogram(TFileDirectory & fDir, const std::string & fName, const std::string & fTitle,
    int fNbins, double fXmin, double fXmax, int fNYbins, double fYmin, double fYmax) const {
    char title[1024];

    sprintf(title, "%s [RUN:%ld/%ld]", fTitle.c_str(), runBegin, lumibegin);
    return fDir.make < TH2F > (fName.c_str(), title, fNbins, fXmin, fXmax,fNYbins,fYmin,fYmax);
}

   void TrigScaleFactors_OnBParking::bookHistograms() {
   TFileDirectory Dir2  = fs->mkdir("Reco_Jets");

//********************************** NTuples ************************************************

    _Tot_Wgt = book1DHistogram(Dir2,"Tot_Wgt","Tot_Wgt",20000,-10000,10000);

    _Tag_Candidate    = new TNtuple("Tag_Candidate","Tag_Candidate","mu_pt:mu_energy:mu_eta:mu_phi:mu_charge:mu_dxy:mu_dxyErr:mu_dz:mu_LooseID:mu_TightID:mu_isGlobal:Gweight:Tnpv");
    _Tag_Candidate_Track = new TNtuple("Tag_Candidate_Track", "Tag_Candidate_Track", "n_track_layers:n_strip_layers:n_pixel_layers:mu_PV_theta:track_purity:drSV:probSV:dr_gen");
    _Probe_Candidate  = new TNtuple("Probe_Candidate","Probe_Candidate","mu_pt:mu_energy:mu_eta:mu_phi:mu_charge:mu_dxy:mu_dxyErr:mu_dz:mu_LooseID:mu_TightID:mu_isGlobal");
    _Probe_Candidate_Track = new TNtuple("Probe_Candidate_Track", "Probe_Candidate_Track", "n_track_layers:n_strip_layers:n_pixel_layers:mu_PV_theta:track_purity:drSV:probSV:dr_gen");

    _L3_0    = new TNtuple("L3_0","L3_0","mu_pt:mu_eta:mu_phi:HLT_Mu12_IP6:HLT_Mu10p5_IP3p5:HLT_Mu9_IP6:HLT_Mu9_IP5:HLT_Mu9_IP4:HLT_Mu8p5_IP3p5:HLT_Mu8_IP6:HLT_Mu8_IP5:HLT_Mu8_IP3:HLT_Mu7_IP4");
    _L3_1    = new TNtuple("L3_1","L3_1","mu_pt:mu_eta:mu_phi:HLT_Mu12_IP6:HLT_Mu10p5_IP3p5:HLT_Mu9_IP6:HLT_Mu9_IP5:HLT_Mu9_IP4:HLT_Mu8p5_IP3p5:HLT_Mu8_IP6:HLT_Mu8_IP5:HLT_Mu8_IP3:HLT_Mu7_IP4");
    _L3_2    = new TNtuple("L3_2","L3_2","mu_pt:mu_eta:mu_phi:HLT_Mu12_IP6:HLT_Mu10p5_IP3p5:HLT_Mu9_IP6:HLT_Mu9_IP5:HLT_Mu9_IP4:HLT_Mu8p5_IP3p5:HLT_Mu8_IP6:HLT_Mu8_IP5:HLT_Mu8_IP3:HLT_Mu7_IP4");
    _L3_3    = new TNtuple("L3_3","L3_3","mu_pt:mu_eta:mu_phi:HLT_Mu12_IP6:HLT_Mu10p5_IP3p5:HLT_Mu9_IP6:HLT_Mu9_IP5:HLT_Mu9_IP4:HLT_Mu8p5_IP3p5:HLT_Mu8_IP6:HLT_Mu8_IP5:HLT_Mu8_IP3:HLT_Mu7_IP4");

    _HLT_result = new TNtuple("HLT_result","HLT_result","Prob_SV:HLT_Mu12_IP6:HLT_Mu10p5_IP3p5:HLT_Mu9_IP6:HLT_Mu9_IP5:HLT_Mu9_IP4:HLT_Mu8p5_IP3p5:HLT_Mu8_IP6:HLT_Mu8_IP5:HLT_Mu8_IP3:HLT_Mu7_IP4");

}

// ******************************************************************************************

void TrigScaleFactors_OnBParking::endJob() {
}

DEFINE_FWK_MODULE(TrigScaleFactors_OnBParking);
