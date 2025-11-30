HCAL Phi Symmetry Calibrations

This includes all of the CMSSW and scripts for the phi symmetry (iterative method) calibrations.

This package uses an EDAnalyzer that should be built within CMSSW. The recommended release is CMSSW_15_1_0. Other releases should probably work, but caveat emptor.
The directory structure needs to be CMSSW_release/src/PhiSym/HCALPhiSymmetry, so you will want to do something like:
```
cmsrel CMSSW_15_0_15_patch4
cd CMSSW_15_0_15_patch4/src
cmsenv
mkdir PhiSym
cd PhiSym
git clone https://github.com/johnpaulchou/HCALPhiSymmetry.git
scram b -j4
cd HCALPhiSymmetry
cmsRun test/phisymtree.py
```
The above commands will read an EDM file and generate a TTree that contains HB, HE, and HF hits above a certain energy threshold and the 3-vector of the selected triggers. The code that does this is found in src/phiSymTree.cc, and the file that controls the parameters is test/phisymtree.py.

The next step is to use the information in the TTree to generate histograms of the hit energy distributions. These histograms are weighted by the 3-vector (eta-phi) direction.
