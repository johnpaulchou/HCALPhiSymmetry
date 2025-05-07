import CRABClient
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'MUON0_2024I_singletest_2'
config.General.workArea = ''


config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'd101XphiSymRECO.py'
config.JobType.inputFiles = ['merged_new.json']
config.JobType.maxMemoryMB = 3000

config.section_("Data")
#config.Data.inputDataset = '/EGamma0/Run2024C-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
#config.Data.inputDataset = '/SingleNeutrino_E-10-Egun/Run3Winter24Reco-NoPU_BPixV5_NZS_133X_mcRun3_2023_realistic_postBPix_v5-v2/GEN-SIM-RECO'
#config.Data.inputDataset ='/SingleMuon/Run2022C-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
#config.Data.inputDataset ='/EGamma0/Run2024G-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
#config.Data.inputDataset ='/EGamma1/Run2024G-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
#config.Data.inputDataset ='/EGamma0/Run2024I-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
#config.Data.inputDataset ='/EGamma1/Run2024I-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'

#config.Data.inputDataset = '/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24DRPremix-140X_mcRun3_2024_realistic_v26-v2/AODSIM'
#config.Data.inputDataset ='/Muon0/Run2024G-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
#config.Data.inputDataset ='/Muon1/Run2024G-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
config.Data.inputDataset ='/Muon0/Run2024I-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
#config.Data.inputDataset ='/Muon1/Run2024I-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = '180'
config.Data.totalUnits = '1000'
#config.Data.lumiMask = 'Cert_Collisions2024_378981_379866_Golden_limfrom379347.json'
#config.Data.lumiMask = 'Cert_Collisions2022_355100_362760_Golden.json'
config.Data.lumiMask = 'Cert_Collisions2024_378981_386951_Golden.json'
config.Data.outLFNDirBase= '/store/user/sdonnell/met_veto'


config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.whitelist = ['T2_CH_CERN','T2_US_Caltech','T2_UK_London_IC','T2_US_Purdue','T2_US_Nebraska']
#config.Site.whitelist = ['T2_CH_CERN','T2_US_Caltech','T2_UK_London_IC']



