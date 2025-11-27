#!/usr/bin/env python3

import ROOT
import common
import argparse
import fnmatch
from collections import defaultdict



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
    parser.add_argument("--no-reweight", action='store_false', dest='reweight', help='Disable reweighting (i.e. the default is to _perform_ reweighting)')
    parser.add_argument("--force-hlt-path", dest='hltpath', help='Name of the HLT path to reweigh against. If not set, the leading trigger object is used from any path. If you want to do matching, make sure to put \' \' around the pattern, e.g. \'electron*_v*\'.')
    args=parser.parse_args()


    # create histogram dictionary
    hdict = nested_dict()

    # create TChain out of files
    chain = ROOT.TChain("phisym/phisymtree")
    for file in args.inputfns:
        chain.Add(file)
    nentries = chain.GetEntries()
    print("Total entries in TChain:", nentries)

    # create an output file
    rootfileout = ROOT.TFile(args.outputfn, "RECREATE")

    # if we are going to do reweighting, we need to create an eta-phi heat map of the trigger object, first
    if args.reweight:

        # create TH2D histogram for the trigger eta-phi
        histtitle="leading trigger object location"
        if args.hltpath is not None: histtitle += " for path="+args.hltpath
        hTrigEtaPhi = ROOT.TH2D('hTrigEtaPhi', histtitle, 5, -2.6, 2.6, 5, -3.1415, 3.1415)

        # loop over events in the chain
        for ievt in range(nentries):
            chain.GetEntry(ievt)

            # get the index of the firing trigger, and fill the histogram, if found
            trigIndex = find_lead_trigger_index(chain.nTrigs, chain.trigPt, chain.trigNames, args.hltpath)
            if trigIndex<0: continue
            hTrigEtaPhi.Fill(chain.trigEta[trigIndex], chain.trigPhi[trigIndex])


    # normalize the 2d histogram
    normalize=hTrigEtaPhi.GetEntries()/hTrigEtaPhi.GetNbinsX()/hTrigEtaPhi.GetNbinsY()
    hTrigEtaPhi.Scale(1./normalize)
    
    # loop over events in the chain
    for ievt in range(nentries):
        chain.GetEntry(ievt)

        # if we're re-weighting, find the index of the trigger, and get the weight for the hits
        hitweight=1.0
        if args.reweight:
            trigIndex = find_lead_trigger_index(chain.nTrigs, chain.trigPt, chain.trigNames, args.hltpath)
            if trigIndex<0: continue # skip the event, if the trigger can't be found
            xbin=hTrigEtaPhi.GetXaxis().FindBin(chain.trigEta[trigIndex])
            ybin=hTrigEtaPhi.GetYaxis().FindBin(chain.trigPhi[trigIndex])
            value = hTrigEtaPhi.GetBinContent(xbin,ybin)
            if value==0: hitweight==0
            else:        hitweight=1./value
                             
        # loop over the hits
        for i in range(chain.nHits):
            
            # get the histogram and fill it
            mod = chain.event % common.modulus
            hist = get_hist_from_dict(hdict, chain.subdet[i], chain.ieta[i], chain.iphi[i], chain.depth[i], mod)
            hist.Fill(chain.hitEnergy[i], hitweight)

    # write all of the histograms to a file
    rootfileout.cd()
    if args.reweight: hTrigEtaPhi.Write()
    write_all_histograms(hdict, rootfileout)
    rootfileout.Close()


    # test the file for the channels
    common.testGoodChannel(args.outputfn)
            
### end main()




if __name__ == "__main__":
    main()
