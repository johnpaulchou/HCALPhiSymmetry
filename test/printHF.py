import FWCore.ParameterSet.Config as cms

process = cms.Process("HFGEOM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")


# Load DB-backed calo geometry readers
process.load("Geometry.CaloEventSetup.CaloGeometryDBReader_cfi")

# For this analyzer we only need HCAL, which includes HF
process.CaloGeometryBuilder.SelectedCalos = cms.vstring("HCAL")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:run3_data", "")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

process.printHFGeometry = cms.EDAnalyzer("PrintHFGeometry")
process.p = cms.Path(process.printHFGeometry)
