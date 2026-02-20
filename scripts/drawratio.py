#!/usr/bin/env python3

import ROOT
import argparse

def main():

    # setup parser
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("num_filename", help="filename with the numerator histograms")
    parser.add_argument("den_filename", help="filename with the denominator histograms")
    args=parser.parse_args()

    out = ROOT.TFile("ratio.root", "RECREATE")

    f1 = ROOT.TFile.Open(args.num_filename)
    f2 = ROOT.TFile.Open(args.den_filename)

    if not f1 or f1.IsZombie(): sys.exit(f"Error opening {args.num_filename}")
    if not f2 or f2.IsZombie(): sys.exit(f"Error opening {args.den_filename}")

    # loop over all keys present in f1
    for key in f1.GetListOfKeys():
        obj = key.ReadObj()

        # Only consider TH2D histograms
        if not obj.InheritsFrom("TH2"): continue

        # get the object from f2 by the name of the first object
        name = obj.GetName()
        h1 = obj
        h2 = f2.Get(name)
        if not h2:
            print(f"Warning: {name} not found in second file, skipping")
            continue

        h_ratio = h1.Clone(f"{name}_ratio")
        h_ratio.SetTitle(f"{name} ratio;{h1.GetXaxis().GetTitle()};{h1.GetYaxis().GetTitle()}")
        h_ratio.Reset()

        h_ratio.Divide(h1, h2, 1.0, 1.0, "")

        # -----------------------------
        # Optional: protect against zero denominator
        # -----------------------------
        for ix in range(1, h_ratio.GetNbinsX() + 1):
            for iy in range(1, h_ratio.GetNbinsY() + 1):
                if h2.GetBinContent(ix, iy) == 0:
                    h_ratio.SetBinContent(ix, iy, 0.0)
                    h_ratio.SetBinError(ix, iy, 0.0)

        out.cd()
        h_ratio.SetMinimum(0.85)
        h_ratio.SetMaximum(1.15)
        h_ratio.Write()

    # end loop over keys
    
    out.Close()
    f1.Close()
    f2.Close()



if __name__ == "__main__":
    main()
