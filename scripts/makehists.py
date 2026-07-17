#!/usr/bin/env python3

import ROOT
import common
import argparse
import fnmatch
import array
from collections import defaultdict
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay, cKDTree
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



### caches the histogram lookup after the first access to each channel/modulus
def get_cached_hist(hdict, hist_cache, subdetnum, ieta, iphi, depth, mod):
    key = (subdetnum, ieta, iphi, depth, mod)
    hist = hist_cache.get(key)
    if hist is None:
        hist = get_hist_from_dict(hdict, subdetnum, ieta, iphi, depth, mod)
        hist_cache[key] = hist
    return hist
### end get_cached_hist()



### enables only the branches needed for the current event loop
def set_active_branches(chain, branch_names):
    chain.SetBranchStatus("*", 0)
    for branch_name in branch_names:
        chain.SetBranchStatus(branch_name, 1)
### end set_active_branches()



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
    parser.add_argument("--numevents", type=int, default=-1, help="number of events to process")
    args=parser.parse_args()

    # load the HF geometry JSON file
    if args.pvcorr:
        with open(args.hfgeofile) as f:
            hfdata = json.load(f)

        # Group cells by depth and z side. The coordinates, areas, Delaunay
        # triangulations, and nearest-neighbor trees are fixed for all events.
        hf_groups = {
            (1, +1): {"cells": []},
            (1, -1): {"cells": []},
            (2, +1): {"cells": []},
            (2, -1): {"cells": []},
        }
        for cell in hfdata["hf_cells"]:
            cell_key = (cell["ieta"], cell["iphi"], cell["depth"])
            zside = +1 if cell["ieta"] > 0 else -1
            group_key = (cell["depth"], zside)
            if group_key not in hf_groups:
                continue

            hf_groups[group_key]["cells"].append((
                cell_key,
                cell["x_cm"],
                cell["y_cm"],
                cell["z_cm"],
                cell["avgArea_cm2"]
            ))

        # Map each HF channel directly to its group and array index. Resetting
        # the density arrays to zero each event accounts for unrecorded hits.
        hf_hit_lookup = {}
        for group_key, group in hf_groups.items():
            cells = group.pop("cells")
            group["keys"] = [cell[0] for cell in cells]
            group["points"] = np.asarray([(cell[1], cell[2]) for cell in cells], dtype=np.float64)
            z_values = np.asarray([cell[3] for cell in cells], dtype=np.float64)
            if not np.allclose(z_values, z_values[0]):
                raise ValueError("HF group {} contains multiple z planes".format(group_key))
            group["z_cm"] = float(z_values[0])
            group["areas"] = np.asarray([cell[4] for cell in cells], dtype=np.float64)
            group["densities"] = np.zeros(len(cells), dtype=np.float64)
            group["triangulation"] = Delaunay(group["points"])
            group["nearest_tree"] = cKDTree(group["points"])

            for index, cell_key in enumerate(group["keys"]):
                hf_hit_lookup[cell_key] = (group_key, index)
            
    # create histogram dictionary
    hdict = nested_dict()
    hist_cache = {}

    # create TChain out of files
    chain = ROOT.TChain("phisym/phisymtree")
    for file in args.inputfns:
        chain.Add(file)
    total_entries = chain.GetEntries()
    nentries = total_entries if args.numevents<0 else min(args.numevents, total_entries)
    print("Total entries to loop over: ", nentries)

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
    
    main_branches = [
        "event", "nHits", "hitEnergy", "ieta", "iphi", "depth", "subdet"
    ]
    if args.reweight:
        main_branches.extend(["nTrigs", "trigPt", "trigEta", "trigPhi", "trigNames"])
    if args.pvcorr:
        main_branches.extend(["pvx", "pvy", "pvz"])
    set_active_branches(chain, main_branches)

    chain.SetBranchAddress("event", event)    
    chain.SetBranchAddress("nHits", nHits)
    chain.SetBranchAddress("hitEnergy", hitEnergy)
    chain.SetBranchAddress("ieta",      ieta)
    chain.SetBranchAddress("iphi",      iphi)
    chain.SetBranchAddress("depth",     depth)
    chain.SetBranchAddress("subdet",    subdet)
    if args.reweight:
        chain.SetBranchAddress("nTrigs",  nTrigs)
        chain.SetBranchAddress("trigPt",  trigPt)
        chain.SetBranchAddress("trigEta", trigEta)
        chain.SetBranchAddress("trigPhi", trigPhi)
        chain.SetBranchAddress("trigNames", trigNames)
    if args.pvcorr:
        chain.SetBranchAddress("pvx", pvx)
        chain.SetBranchAddress("pvy", pvy)
        chain.SetBranchAddress("pvz", pvz)

    # create an output file
    rootfileout = ROOT.TFile(args.outputfn, "RECREATE")

    # if we are going to do reweighting, we need to create an eta-phi heat map of the trigger object, first
    if args.reweight:

        # The first pass only uses trigger branches.
        trigger_branches = ["nTrigs", "trigPt", "trigEta", "trigPhi", "trigNames"]
        set_active_branches(chain, trigger_branches)

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

        # Restore all branches needed by the main event loop.
        set_active_branches(chain, main_branches)


    # loop over events in the chain
    for ievt in range(nentries):
        if ievt%10000==0: print("2nd loop: "+str(ievt))
        chain.GetEntry(ievt)
        mod = event[0] % common.modulus

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
            if pvx[0]<-900. or pvy[0]<-900. or pvz[0]<-900.: continue

            # Start every event with zero density in every HF cell. Recorded hits below overwrite the corresponding entries
            for group in hf_groups.values(): group["densities"].fill(0.0)

            # Loop over the recorded HF hits and fill the density arrays
            for i in range(nHits[0]):
                if subdet[i]!=4: continue
                cell_key = (ieta[i], iphi[i], depth[i])
                location = hf_hit_lookup.get(cell_key)
                if location is None:
                    raise KeyError("No HF geometry for cell {} in event {}".format(cell_key, event[0]))
                group_key, index = location
                group = hf_groups[group_key]
                group["densities"][index] = hitEnergy[i]/group["areas"][index]

            # Evaluate each of the four interpolation groups in one call. The values and PV-shifted query points are independent for each event
            pv_transverse = np.asarray((pvx[0], pvy[0]), dtype=np.float64)
            for group in hf_groups.values():
                z_plane = group["z_cm"]
                scale = 1.0 - pvz[0]/z_plane
                query_points = pv_transverse + scale*group["points"]
                linear_interp = LinearNDInterpolator(group["triangulation"],group["densities"],fill_value=np.nan)
                new_densities = np.asarray(linear_interp(query_points), dtype=np.float64).reshape(-1)

                # Use the nearest source cell for points outside the hull
                outside = np.isnan(new_densities)
                if np.any(outside):
                    nearest_indices = group["nearest_tree"].query(query_points[outside], k=1)[1]
                    new_densities[outside] = group["densities"][nearest_indices]

                new_energies = new_densities*group["areas"]*scale**2
                for cell_key, newenergy in zip(group["keys"], new_energies):
                    cell_ieta, cell_iphi, cell_depth = cell_key
                    hist = get_cached_hist(hdict, hist_cache, 4, cell_ieta, cell_iphi, cell_depth, mod)
                    hist.Fill(float(newenergy), hitweight)
            
        # loop over the hits
        for i in range(nHits[0]):
            
            # get the histogram and fill it (skip the HF if we are correcting for the PV location)
            if not args.pvcorr or subdet[i]!=4:
                hist = get_cached_hist(hdict, hist_cache, subdet[i], ieta[i], iphi[i], depth[i], mod)
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
