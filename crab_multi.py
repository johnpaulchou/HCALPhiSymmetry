import CRABClient
from WMCore.Configuration import Configuration
import os

full_run_range = (383811, 385801) # Example: from run 321000 to run 321500
num_parts = 4
run_ranges = [
    (383811, 384231),
    (384238, 384935),
    (384950, 385235),
    (385383, 385801)
]


for start_run, end_run in run_ranges:
        request_Name = f'MUON1_Run2024G_{start_run}_{end_run}_METFilters'
        config = Configuration()
        config.section_("General")
        config.General.requestName = request_Name
        config.General.workArea = ''


        config.section_("JobType")
        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = 'd101XphiSymRECO.py'
        config.JobType.inputFiles = ''
        config.JobType.maxMemoryMB = 3000
        
        config.section_("Data")
        #config.Data.inputDataset = '/EGamma0/Run2024C-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
        #config.Data.inputDataset = '/SingleNeutrino_E-10-Egun/Run3Winter24Reco-NoPU_BPixV5_NZS_133X_mcRun3_2023_realistic_postBPix_v5-v2/GEN-SIM-RECO'
        #config.Data.inputDataset ='/SingleMuon/Run2022C-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
        #config.Data.inputDataset ='/EGamma0/Run2024G-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
        #config.Data.inputDataset ='/EGamma1/Run2024G-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
        #config.Data.inputDataset ='/EGamma0/Run2024I-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
        #config.Data.inputDataset ='/EGamma1/Run2024I-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
        #config.Data.inputDataset ='/Muon0/Run2024G-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
        config.Data.inputDataset ='/Muon1/Run2024G-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
        #config.Data.inputDataset ='/Muon0/Run2024I-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
        #config.Data.inputDataset ='/Muon1/Run2024I-HcalCalIterativePhiSym-PromptReco-v1/ALCARECO'
        config.Data.inputDBS = 'global'
        config.Data.splitting = 'Automatic'
        #config.Data.unitsPerJob = '180'
        config.Data.totalUnits = '-1'
        #config.Data.lumiMask = 'Cert_Collisions2024_378981_379866_Golden_limfrom379347.json'
        #config.Data.lumiMask = 'Cert_Collisions2022_355100_362760_Golden.json'
        config.Data.lumiMask = 'Cert_Collisions2024_378981_386951_Golden.json'
        config.Data.outLFNDirBase= '/store/user/sdonnell/new_iter_w_MET'
        
        config.Data.runRange = f'{start_run}-{end_run}'
        #config.Data.ignoreLocality = True
        config.Data.partialDataset = True
        config.section_("Site")
        config.Site.storageSite = 'T3_US_FNALLPC'
        #config.Site.maxRunningJobs = '10'
        config.Site.whitelist = ['T2_CH_CERN','T2_US_Caltech','T2_UK_London_IC','T2_US_Purdue','T2_US_Nebraska']
        #config.Site.whitelist = ['T2_CH_CERN','T2_US_Caltech','T2_UK_London_IC']


        config_file = f'config_{start_run}_{end_run}.py'
        with open(config_file, 'w') as f:
                f.write(str(config))

        print(f"Submitting CRAB job for run range {start_run}-{end_run}...")
        os.system(f'crab submit -c {config_file}')
