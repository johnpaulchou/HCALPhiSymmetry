import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('USER',Run3)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32 (-1)
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.default.limit = 10

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryDB_cff")


#--- Global Tag conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


#process.GlobalTag.globaltag = '130X_dataRun3_Prompt_v3'
process.GlobalTag.globaltag = '140X_dataRun3_Prompt_v4'

process.phaseHF = cms.EDAnalyzer ("phiSym",
                                  textFile = cms.untracked.string('yhisto_6.txt'),
                                  hfreco =  cms.InputTag("hfreco"),
                                  triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                                  #triggerNamesSingleMu = cms.untracked.vstring("HLT_IsoMu24", "HLT_Mu50"),
                                 # triggerNamesDoubleMu = cms.untracked.vstring("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"),
                                  hbhereco = cms.InputTag("hbhereco")
 )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('phisym.root'),
)

#----------------------------
# Paths/Sequences Definitions
#----------------------------
process.p = cms.Path(
    process.phaseHF
)

process.source = cms.Source ("PoolSource" ,
                             fileNames=cms.untracked.vstring(
#'/store/data/Commissioning2021/MinimumBias9/ALCARECO/HcalCalIterativePhiSym-PromptReco-v1/000/346/169/00000/1e53514c-b551-40e7-b35d-5fb542ef8180.root'
#'/store/data/Run2024C/EGamma0/ALCARECO/HcalCalIterativePhiSym-PromptReco-v1/000/379/415/00000/08bee99c-d6c9-4ab8-a7b1-dfc4f2aaaf34.root',
#'/store/data/Run2024C/EGamma0/ALCARECO/HcalCalIterativePhiSym-PromptReco-v1/000/379/416/00000/01f03546-17e8-42b6-b745-24b30c7090ec.root',
#'/store/data/Run2024C/EGamma0/ALCARECO/HcalCalIterativePhiSym-PromptReco-v1/000/379/765/00000/94ce8d5e-721e-4332-bc2c-ac61745e8bce.root'
#'/store/data/Run2024C/EGamma0/ALCARECO/HcalCalIterativePhiSym-PromptReco-v1/000/379/416/00000/510a5706-91f2-491d-8c5b-c7f8cb606d9c.root'
'/store/data/Run2024I/Muon0/ALCARECO/HcalCalIterativePhiSym-PromptReco-v1/000/386/478/00000/00bf7983-375f-461c-b23d-9519fe075845.root'
#'/store/data/Run2024C/EGamma0/ALCARECO/HcalCalIterativePhiSym-PromptReco-v1/000/379/617/00000/1126e275-f0d9-4dd6-8662-a27fe3562950.root'
#'/store/data/Run2024C/EGamma0/ALCARECO/HcalCalIterativePhiSym-PromptReco-v1/000/379/530/00000/2036b45f-df3a-4c61-bad0-50af2f877e7d.root'
#'/store/data/Run2024C/EGamma0/ALCARECO/HcalCalIterativePhiSym-PromptReco-v1/000/379/416/00000/030fe9e1-9710-4aa3-b775-c1b7499cae15.root'
#'/SingleNeutrino_E-10-Egun/Run3Winter24Reco-NoPU_BPixV5_NZS_133X_mcRun3_2023_realistic_postBPix_v5-v2/GEN-SIM-RECO'
)
)
