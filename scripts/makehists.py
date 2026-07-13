#!/usr/bin/env python3

import ROOT
import common
import argparse
import fnmatch
import array
from collections import defaultdict
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
import numpy as np
import json

### creates a nested dictionary to manage the histograms
def nested_dict():
    return defaultdict(nested_dict)
### end nested_dict()




### gets a histogram from a nested dictionary based on the subdetnumber, ieta, iphi, depth, and mod
### if it doesn't exist, create it first
def get_hist_from_dict(hdict, subdetnum, ieta, iphi, depth, mod):
    node = hdict[subdetnum][ieta][iphi][depth]

    # If histogram doesn't exist, create it first
    if mod not in node:
        name=common.getHistName(subdetnum, ieta, iphi, depth, mod)
        binning=common.binning(subdetnum, ieta)
        node[mod] = ROOT.TH1D(name,name,binning[0],binning[1],binning[2])

    return node[mod]
### end get_hist_from_dict()



### write all of the histograms of the nested dictionary, with a recursive algorithm
def write_all_histograms(d, outfile):
    for k, v in d.items():
        if isinstance(v, ROOT.TH1):
            v.Write()
        else:
            write_all_histograms(v, outfile)
### end write_all_histograms()



### finds the leading trigger index that matches the string, if provided
def find_lead_trigger_index(nTrigs, trigPt, trigNames, hltpath):

    # loop over the triggers
    leadPt=-999.
    leadPtIndex=-1
    for itrig in range(nTrigs):

        # skip if there not a match, otherwise, find the leading one
        if hltpath is not None and not fnmatch.fnmatch(trigNames[itrig], hltpath): continue
        if trigPt[itrig]>leadPt:
            leadPt=trigPt[itrig]
            leadPtIndex=itrig

    return leadPtIndex
### end find_lead_trigger_index()
   




### main function
def main():
    
    # setup parser
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputfns", nargs='+', help="name(s) of .root file(s) to be processed")
    parser.add_argument("-o", "--outputfn", default="hists.root", help="name of the .root file that stores the histograms")
    parser.add_argument("--no-reweight", action='store_false', dest='reweight', help='Disable trigger reweighting')
    parser.add_argument("--force-hlt-path", dest='hltpath', help='Name of the HLT path to reweigh against. If not set, the leading trigger object is used from any path. If you want to do matching, make sure to put \' \' around the pattern, e.g. \'electron*_v*\'.')
    parser.add_argument("--no-pvcorr", action='store_false', dest='pvcorr', help='Disable the HF primary vertex correction')
    parser.add_argument("--hf-geometry-json", default="hf_geometry.json", dest="hfgeofile", help=".json file produced by the \"cmsRun printHF.py\" command")
    args=parser.parse_args()

    # load the HF geometry JSON file
    if args.pvcorr:
        with open(args.hfgeofile) as f:
            hfdata = json.load(f)

        # create an empty dictionary from (ieta,iphi,depth)-->(eta,phi,area,energy)
        emptyhfdict = {}
        for cell in hfdata["hf_cells"]:
            ieta = cell["ieta"]
            iphi = cell["iphi"]
            depth = cell["depth"]
            x = cell["x_cm"]
            y = cell["y_cm"]
            area = cell["etaSpan"]*cell["phiSpan"]
            emptyhfdict[ieta, iphi, depth]=[x,y,area,0.0]
            
    # create histogram dictionary
    hdict = nested_dict()

    # create TChain out of files
    chain = ROOT.TChain("phisym/phisymtree")
    for file in args.inputfns:
        chain.Add(file)
    nentries = chain.GetEntries()
    print("Total entries in TChain:", nentries)


    run    = array.array("L", [0])
    lumi   = array.array("L", [0])
    event  = array.array("L", [0])
    nHits      = array.array("i", [0])
    MAXNHITS   = 20000
    hitEnergy  = array.array("f", [0]*MAXNHITS)
    ieta       = array.array("i", [0]*MAXNHITS)
    iphi       = array.array("i", [0]*MAXNHITS)
    depth      = array.array("i", [0]*MAXNHITS)
    subdet     = array.array("i", [0]*MAXNHITS)
    nTrigs     = array.array("i", [0])
    MAXNTRIGS  = 100
    trigPt     = array.array("f", [0]*MAXNTRIGS)
    trigEta    = array.array("f", [0]*MAXNTRIGS)
    trigPhi    = array.array("f", [0]*MAXNTRIGS)
    trigNames = ROOT.std.vector("string")()
    pvx = array.array("f", [0])
    pvy = array.array("f", [0])
    pvz = array.array("f", [0])

    chain.SetBranchAddress("run",   run)
    chain.SetBranchAddress("lumi",  lumi)
    chain.SetBranchAddress("event", event)    
    chain.SetBranchAddress("nHits", nHits)
    chain.SetBranchAddress("hitEnergy", hitEnergy)
    chain.SetBranchAddress("ieta",      ieta)
    chain.SetBranchAddress("iphi",      iphi)
    chain.SetBranchAddress("depth",     depth)
    chain.SetBranchAddress("subdet",    subdet)
    chain.SetBranchAddress("nTrigs",  nTrigs)
    chain.SetBranchAddress("trigPt",  trigPt)
    chain.SetBranchAddress("trigEta", trigEta)
    chain.SetBranchAddress("trigPhi", trigPhi)
    chain.SetBranchAddress("trigNames", trigNames)
    chain.SetBranchAddress("pvx", pvx)
    chain.SetBranchAddress("pvy", pvy)
    chain.SetBranchAddress("pvz", pvz)

    # create an output file
    rootfileout = ROOT.TFile(args.outputfn, "RECREATE")

    # if we are going to do reweighting, we need to create an eta-phi heat map of the trigger object, first
    if args.reweight:

        # create TH2D histogram for the trigger eta-phi
        histtitle="leading trigger object location"
        if args.hltpath is not None: histtitle += " for path="+args.hltpath
        hTrigEtaPhi = ROOT.TH2D('hTrigEtaPhi', histtitle, 15, -2.6, 2.6, 50, -3.1416, 3.1416)
        hWeightEtaPhi = ROOT.TH2D('hWeightEtaPhi', histtitle, 15, -2.6, 2.6, 50, -3.1416, 3.1416)

        # loop over events in the chain
        for ievt in range(nentries):
            if ievt%10000==0: print("1st loop: "+str(ievt))
            chain.GetEntry(ievt)

            # get the index of the firing trigger, and fill the histogram, if found
            trigIndex = find_lead_trigger_index(nTrigs[0], trigPt, trigNames, args.hltpath)
            if trigIndex<0: continue
            hTrigEtaPhi.Fill(trigEta[trigIndex], trigPhi[trigIndex])

        # produce the weight histogram from the trigger map by normalizing the 2d histogram by eta slice
        for i in range(1,hTrigEtaPhi.GetNbinsX()+1):
            sum=0.
            for j in range(1,hTrigEtaPhi.GetNbinsY()+1):
                sum += 1./hTrigEtaPhi.GetBinContent(i,j)
            for j in range(1,hTrigEtaPhi.GetNbinsY()+1):
                val = 1./hTrigEtaPhi.GetBinContent(i,j)
                hWeightEtaPhi.SetBinContent(i,j,val/sum*hTrigEtaPhi.GetNbinsY())


    # loop over events in the chain
    for ievt in range(nentries):
        if ievt%10000==0: print("2nd loop: "+str(ievt))
        chain.GetEntry(ievt)

        # if we're re-weighting, find the index of the trigger, and get the weight for the hits
        hitweight=1.0
        if args.reweight:
            trigIndex = find_lead_trigger_index(nTrigs[0], trigPt, trigNames, args.hltpath)
            if trigIndex<0: continue # skip the event, if the trigger can't be found

            # skip events for which the trigger shows up in the underflow/overflow
            xbin=hWeightEtaPhi.GetXaxis().FindBin(trigEta[trigIndex])
            ybin=hWeightEtaPhi.GetYaxis().FindBin(trigPhi[trigIndex])
            if xbin<=0 or xbin>=hWeightEtaPhi.GetNbinsX()+1 or ybin<=0 or ybin>=hWeightEtaPhi.GetNbinsY()+1: continue
            hitweight = hWeightEtaPhi.GetBinContent(xbin,ybin)


        # perform the primary vertex correction for the HF hits
        if args.pvcorr:

            # check that there is a valid PV (skip the event if not)
            if pvx[0]<-900. or pvy[0]<-900.: continue

            # create a new hf dictionary that is a copy of the empty one
            hfdict = {key: value.copy() for key, value in emptyhfdict.items()}

            # loop over just the HF hits and fill the dictionary
            for i in range(nHits[0]):
                if subdet[i]!=4: continue
                key = (ieta[i], iphi[i], depth[i])
                hfdict[key][3] = hitEnergy[i]/hfdict[key][2]

            # create 4 interpolators: depth 1/2 times ieta +/- side
            interp_inputs = {
                (1, +1): [[], []],
                (1, -1): [[], []],
                (2, +1): [[], []],
                (2, -1): [[], []],
            }

            # fill the interpolator inputs
            for key, vals in hfdict.items():
                cell_ieta, cell_iphi, cell_depth = key
                x, y, area, energydensity = vals
                zside = +1 if cell_ieta > 0 else -1

                interp_key = (cell_depth, zside)
                if interp_key not in interp_inputs: continue

                interp_inputs[interp_key][0].append((x, y))
                interp_inputs[interp_key][1].append(energydensity)

            # create the interpolators
            hf_linear_interps = {}
            hf_nearest_interps = {}
            for interp_key, (points, values) in interp_inputs.items():
                hf_linear_interps[interp_key] = LinearNDInterpolator(points, values)
                hf_nearest_interps[interp_key] = NearestNDInterpolator(points, values)

            # fill the histograms
            for key, vals in hfdict.items():
                cell_ieta, cell_iphi, cell_depth = key
                x, y, area, energydensity = vals
                zside = +1 if cell_ieta > 0 else -1

                interp_key = (cell_depth, zside)
                if interp_key not in interp_inputs: continue

                newdensity = hf_linear_interps[interp_key](x + pvx[0], y + pvy[0])
                if np.isnan(newdensity):
                    newdensity = hf_nearest_interps[interp_key](x + pvx[0], y + pvy[0])
                newenergy = float(newdensity) * area
                oldenergy = energydensity * area

                # get the hist
                mod = event[0] % common.modulus
                hist = get_hist_from_dict(hdict, 4, cell_ieta, cell_iphi, cell_depth, mod)
                hist.Fill(newenergy, hitweight)
#                print("ieta="+str(cell_ieta)+"; iphi="+str(cell_iphi)+"; depth="+str(cell_depth)+"; old energy="+str(oldenergy)+"; new energy="+str(newenergy))
            
        # loop over the hits
        for i in range(nHits[0]):
            
            # get the histogram and fill it (skip the HF if we are correcting for the PV location)
            if not args.pvcorr or subdet[i]!=4:
                mod = event[0] % common.modulus
                hist = get_hist_from_dict(hdict, subdet[i], ieta[i], iphi[i], depth[i], mod)
                hist.Fill(hitEnergy[i], hitweight)
        
    # write all of the histograms to a file
    rootfileout.cd()
    if args.reweight:
        hTrigEtaPhi.Write()
        hWeightEtaPhi.Write()
    write_all_histograms(hdict, rootfileout)
    rootfileout.Close()


    # test the file for the channels
    common.testGoodChannel(args.outputfn)
            
### end main()




if __name__ == "__main__":
    main()
