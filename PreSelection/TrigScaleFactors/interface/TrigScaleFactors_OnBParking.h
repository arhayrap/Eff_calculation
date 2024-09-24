#include <memory>
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include <cmath>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/TH1Store.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include <TStyle.h>
#include "TMath.h"
#include <stdint.h>
#include "TTree.h"
#include <memory>
#include <fstream>
#include <math.h>
#include "TNtuple.h"
#include "TH1I.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TProfile.h"

// Trigger-Objects
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1Trigger/interface/Jet.h"

class TH1F;
class TH2F;
class TStyle;
class TTree;

class TrigScaleFactors_OnBParking : public edm::one::EDAnalyzer<> {
   public:
      explicit TrigScaleFactors_OnBParking(const edm::ParameterSet&);
      ~TrigScaleFactors_OnBParking();
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
//      virtual void endRun(const edm::Run&, const edm::EventSetup&) override;
      virtual void endJob() override;
      virtual double calculateBackpointingAngle(const reco::Vertex&, const reco::VertexCompositePtrCandidate&, const pat::Muon&);
      virtual double calculateMuonSVProbability(const pat::Muon&, const std::vector<reco::Vertex>&);
      virtual double DeltaR(double eta1, double eta2, double phi1, double phi2);

  void bookHistograms ();
  TH1F* book1DHistogram (TFileDirectory& fDir, const std::string& fName, const std::string& fTitle, 
			 int fNbins, double fXmin, double fXmax) const;
  TH1F* book1DHistogram ( const std::string& fName, const std::string& fTitle, 
			 int fNbins, double fXmin, double fXmax) const;

  TH2F* book2DHistogram (TFileDirectory& fDir, const std::string& fName, const std::string& fTitle, 
			 int fNbinsX, double fXmin, double fXmax,
			 int fNbinsY, double fYmin, double fYmax) const;
  TProfile* bookProfileHistogram (TFileDirectory & fDir, const std::string & fName, const std::string & fTitle, int fNbins, double fXmin, double fXmax, double fYmin, double fYmax) const;

  edm::Service<TFileService> fs;
  
//***************** NTuples *****************
      TH1F* _Tot_Wgt;

      TNtuple* _Tag_Candidate;
      TNtuple* _Probe_Candidate;
      TNtuple* _Tag_Candidate_Track;
      TNtuple* _Probe_Candidate_Track;

      TNtuple* _L1_0;
      TNtuple* _L3_0;
      TNtuple* _L3_1;
      TNtuple* _L3_2;
      TNtuple* _L3_3;

      TNtuple* _HLT_result;

//********************************************

      edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
      edm::EDGetToken l1GtToken_;
      edm::EDGetTokenT<BXVector<l1t::Jet>> l1jetToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      // edm::EDGetTokenT<reco::VertexCompositePtrCandidate> secondary_vtxToken_;
      edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate>> secondaryVerticesToken_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;
      edm::EDGetTokenT<edm::View<pat::Jet>> recjetToken_;
      edm::EDGetTokenT<edm::View<pat::Muon>> recmuonToken_;
      long runBegin,lumibegin,lumiend,evtNo, goodevtNo, trigg, trigg_tmp, flagECevtNo;
};
