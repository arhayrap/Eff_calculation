import FWCore.ParameterSet.Config as cms
process = cms.Process("newPAT")
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patAlgosToolsTask = getPatAlgosToolsTask(process)
import sys

from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v24', '')

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.EventContent.EventContentCosmics_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

############################################# Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     sys.argv[1],
     sys.argv[2],
     sys.argv[3],
     sys.argv[4],
     sys.argv[5]
     # "root://cms-xrd-global.cern.ch//store/data/Run2018D/ParkingBPH1/MINIAOD/PromptReco-v2/000/321/712/00000/02132230-F1A8-E811-8187-FA163E2C7D31.root"
    )
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


###########################  Analyser
process.Training = cms.EDAnalyzer('TrigScaleFactors_OnBParking',
                               # GenInf     = cms.InputTag("generator", "", ""),
                               bits       = cms.InputTag("TriggerResults","","HLT"),
                               prescales  = cms.InputTag("patTrigger"),
                               objects    = cms.InputTag("slimmedPatTrigger"),
                               # vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                               secondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
                               # vertices   = cms.InputTag("offlinePrimaryVertices", "", "RECO"),
                               offlinePrimaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices", "", "RECO"),  # Ensure this matches your input file
                               recMuon    = cms.InputTag("slimmedMuons"),
                               recJetAK8  = cms.InputTag("slimmedJetsAK8"),
                               recJet     = cms.InputTag("slimmedJets"),
                               pileupInfo = cms.InputTag("slimmedAddPileupInfo"),
                               l1jetSrc   = cms.InputTag("caloStage2Digis:Jet"),
                               l1GtSrc    = cms.InputTag("gtStage2Digis"),
)
###########################   Output
process.TFileService = cms.Service("TFileService",
    fileName  = cms.string("Trigger_ScaleFactors_OnBParking.root")
)
###########################   Schedule
process.p = cms.Path(process.Training, patAlgosToolsTask)
