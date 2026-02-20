#!/usr/bin/env python3

import ROOT  # type: ignore
import argparse

def main():

    # setup parser
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputfilename", help="name of file containing the corrections")
    parser.add_argument("--rootfilename", help="name of .root file to write out histograms to", default="hists.root")
    parser.add_argument("--metanote", help="meta note to save in the .root file")
    args=parser.parse_args()

    
    hists = []
    for depth in range(1, 8):
        histname = "h_corr_depth"+str(depth)
        histtitle = "Correction factor for depth "+str(depth)
        hist = ROOT.TH2D(histname, histtitle, 83, -41.5, 41.5, 72, 0.5, 72.5)
        hists.append(hist)

    # don't assume the file is sorted in any particular way
    with open(args.inputfilename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            
            subdet, ieta, iphi, depth, corr, err = line.split()
            
            ieta = int(ieta)
            iphi = int(iphi)
            depth = int(depth)
            corr = float(corr)
            err = float(err)

            hist=hists[depth-1]
            
            # Find bin numbers explicitly
            bin_x = hist.GetXaxis().FindBin(ieta)
            bin_y = hist.GetYaxis().FindBin(iphi)

            hist.SetBinContent(bin_x, bin_y, corr)
            hist.SetBinError(bin_x, bin_y, err)

    # -----------------------------
    # Save output
    # -----------------------------
    outfile = ROOT.TFile(args.rootfilename, "RECREATE")
    outfile.cd()
    for hist in hists:
        hist.SetStats(0)
        hist.SetMinimum(0.5)
        hist.SetMaximum(1.5)
        hist.GetXaxis().SetTitle("i_{#eta}")
        hist.GetYaxis().SetTitle("i_{#phi}")
        hist.Write()
    if args.metanote:
        metanote = ROOT.TNamed("metanote", args.metanote)
        metanote.Write()
    outfile.Close()    

    
if __name__ == "__main__":
    main()
