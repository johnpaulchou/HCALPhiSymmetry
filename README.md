HCAL Phi Symmetry Calibrations

This includes all of the CMSSW and scripts for the phi symmetry (iterative method) calibrations.

This package uses an EDAnalyzer that should be built within CMSSW. The recommended release is CMSSW_15_0_15_patch4 (or greater).
The directory structure needs to be CMSSW_release/src/PhiSym/HCALPhiSymmetry, so you will want to do something like:
mkrel CMSSW_15_0_15_patch4
cd CMSSW_15_0_15_patch4
cd src
cmsenv
mkdir PhiSym
cd PhiSym
git clone https://github.com/johnpaulchou/HCALPhiSymmetry.git
scram b -j4





source /cvmfs/sft.cern.ch/lcg/views/LCG_107_cuda/x86_64-el8-gcc11-opt/setup.sh