#!/usr/bin/env python3

import ROOT
import common
import argparse
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
        node[mod] = ROOT.TH1D(name,name,100,0,20)

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


### main function
def main():
    
    # setup parser
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputfns", nargs='+', help="name(s) of .root file(s) to be processed")
    parser.add_argument("--outputfn", default="hists.root", help="name of the .root file that stores the histograms")
    parser.add_argument("--errfn", default="makehists.err", help="name of the file to write errors")
    args=parser.parse_args()

    
    # create histogram dictionary
    hdict = nested_dict()

    # create TChain out of files
    chain = ROOT.TChain("phisym/phisymtree")
    for file in args.inputfns:
        chain.Add(file)
    print("Total entries in TChain:", chain.GetEntries())

    # create an output file
    rootfileout = ROOT.TFile(args.outputfn, "RECREATE")

    # loop over events in the chain
    nentries = chain.GetEntries()
    for ievt in range(nentries):
        chain.GetEntry(ievt)

        # loop over the hits
        for i in range(chain.nHits):
            
            # get the histogram and fill it
            mod = chain.event % common.modulus
            hist = get_hist_from_dict(hdict, chain.subdet[i], chain.ieta[i], chain.iphi[i], chain.depth[i], mod)
            hist.Fill(chain.hitEnergy[i])

    # write all of the histograms to a file
    rootfileout.cd()
    write_all_histograms(hdict, rootfileout)
    rootfileout.Close()


    # test the file for the channels
    common.testGoodChannel(args.outputfn, args.errfn)
            
### end main()




if __name__ == "__main__":
    main()
