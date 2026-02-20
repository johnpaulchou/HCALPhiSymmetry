import FWCore.ParameterSet.Config as cms
import glob

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('USER',Run3)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32 (-1)
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.default.limit = 100

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


#process.GlobalTag.globaltag = '130X_dataRun3_Prompt_v3'
process.GlobalTag.globaltag = '140X_dataRun3_Prompt_v4'

process.phisym = cms.EDAnalyzer("phiSymTree",
                                hbheRecHits = cms.InputTag("hbhereco"),
                                hfRecHits = cms.InputTag("hfreco"),
                                triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                                triggerSummary = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                                undoRespCorr = cms.bool(True),
                                HBHEthreshold = cms.double(1.0),
                                HFthreshold = cms.double(5.0),

                                hltPaths = cms.vstring("HLT_Ele*_WPTight_Gsf_v*",
                                                       "HLT_Photon*_TightID_TightIso_v*")
                                )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('/cms/johnpaul/HCAL/phisymtree.root'),
)

process.p = cms.Path(
    process.phisym
    )


files = glob.glob(os.path.join('/cms/bphys/HCAL_final/', "*.root"))
process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring("file:" + f for f in files))

#process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring('file:/cms/bphys/HCAL/HcalCalIterativePhiSym_job00493.root'))

