import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
import sys

process = cms.Process("newPAT")

from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.Skimming.pwdgSkimBPark_cfi')
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_upgrade2018_realistic_v7', '')

patAlgosToolsTask = getPatAlgosToolsTask(process)

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.EventContent.EventContentCosmics_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

############################################# Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     sys.argv[1]
     # 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/N1_102X_upgrade2018_realistic_v15-v1/00000/000C9465-E6D9-C848-8211-E651470DA55A.root'
    )
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

###########################  Analyser
process.Training = cms.EDAnalyzer('TrigScaleFactors_OnMC',
                               GenInf     = cms.InputTag("generator", "", ""),
                               bits       = cms.InputTag("TriggerResults","","HLT"),
                               prescales  = cms.InputTag("patTrigger"),
                               objects    = cms.InputTag("slimmedPatTrigger"),
                               vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                               secondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
                               # vertices   = cms.InputTag("offlinePrimaryVertices"),
                               # primaryVertices = cms.InputTag("offlinePrimaryVertices")  # Ensure this matches your input file
                               recMuon    = cms.InputTag("slimmedMuons"),
                               recJetAK8  = cms.InputTag("slimmedJetsAK8"),
                               recJet     = cms.InputTag("slimmedJets"),
                               pileupInfo = cms.InputTag("slimmedAddPileupInfo"),
                               l1jetSrc   = cms.InputTag("caloStage2Digis:Jet"),
                               l1GtSrc    = cms.InputTag("gtStage2Digis"),
)
###########################   Output
process.TFileService = cms.Service("TFileService",
    fileName  = cms.string("Trigger_ScaleFactors_OnMC.root")
)
###########################   Schedule
process.p = cms.Path(process.Training,patAlgosToolsTask)
