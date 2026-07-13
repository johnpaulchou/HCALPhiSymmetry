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

    # create some other histograms out of the ratio plots
    hHBMD1Res = ROOT.TH1D("hHBMD1Res","HB- Depth=1 residuals",100,0.5,1.5)
    hHBMD2Res = ROOT.TH1D("hHBMD2Res","HB- Depth=2 residuals",100,0.5,1.5)
    hHBMD3Res = ROOT.TH1D("hHBMD3Res","HB- Depth=3 residuals",100,0.5,1.5)
    hHBMD4Res = ROOT.TH1D("hHBMD4Res","HB- Depth=4 residuals",100,0.5,1.5)
    hHBPD1Res = ROOT.TH1D("hHBPD1Res","HB+ Depth=1 residuals",100,0.5,1.5)
    hHBPD2Res = ROOT.TH1D("hHBPD2Res","HB+ Depth=2 residuals",100,0.5,1.5)
    hHBPD3Res = ROOT.TH1D("hHBPD3Res","HB+ Depth=3 residuals",100,0.5,1.5)
    hHBPD4Res = ROOT.TH1D("hHBPD4Res","HB+ Depth=4 residuals",100,0.5,1.5)
    hHEMD1Res = ROOT.TH1D("hHEMD1Res","HE- Depth=1 residuals",100,0.5,1.5)
    hHEMD2Res = ROOT.TH1D("hHEMD2Res","HE- Depth=2 residuals",100,0.5,1.5)
    hHEMD3Res = ROOT.TH1D("hHEMD3Res","HE- Depth=3 residuals",100,0.5,1.5)
    hHEMD4Res = ROOT.TH1D("hHEMD4Res","HE- Depth=4 residuals",100,0.5,1.5)
    hHEMD5Res = ROOT.TH1D("hHEMD5Res","HE- Depth=5 residuals",100,0.5,1.5)
    hHEMD6Res = ROOT.TH1D("hHEMD6Res","HE- Depth=6 residuals",100,0.5,1.5)
    hHEMD7Res = ROOT.TH1D("hHEMD7Res","HE- Depth=7 residuals",100,0.5,1.5)
    hHEPD1Res = ROOT.TH1D("hHEPD1Res","HE+ Depth=1 residuals",100,0.5,1.5)
    hHEPD2Res = ROOT.TH1D("hHEPD2Res","HE+ Depth=2 residuals",100,0.5,1.5)
    hHEPD3Res = ROOT.TH1D("hHEPD3Res","HE+ Depth=3 residuals",100,0.5,1.5)
    hHEPD4Res = ROOT.TH1D("hHEPD4Res","HE+ Depth=4 residuals",100,0.5,1.5)
    hHEPD5Res = ROOT.TH1D("hHEPD5Res","HE+ Depth=5 residuals",100,0.5,1.5)
    hHEPD6Res = ROOT.TH1D("hHEPD6Res","HE+ Depth=6 residuals",100,0.5,1.5)
    hHEPD7Res = ROOT.TH1D("hHEPD7Res","HE+ Depth=7 residuals",100,0.5,1.5)
    hHFMD1Res = ROOT.TH1D("hHFMD1Res","HF- Depth=1 residuals",100,0.5,1.5)
    hHFMD2Res = ROOT.TH1D("hHFMD2Res","HF- Depth=2 residuals",100,0.5,1.5)
    hHFPD1Res = ROOT.TH1D("hHFPD1Res","HF+ Depth=1 residuals",100,0.5,1.5)
    hHFPD2Res = ROOT.TH1D("hHFPD2Res","HF+ Depth=2 residuals",100,0.5,1.5)
    
    hHBMD1Pull = ROOT.TH1D("hHBMD1Pull","HB- Depth=1 pull",100,-15,15)
    hHBMD2Pull = ROOT.TH1D("hHBMD2Pull","HB- Depth=2 pull",100,-15,15)
    hHBMD3Pull = ROOT.TH1D("hHBMD3Pull","HB- Depth=3 pull",100,-15,15)
    hHBMD4Pull = ROOT.TH1D("hHBMD4Pull","HB- Depth=4 pull",100,-15,15)
    hHBPD1Pull = ROOT.TH1D("hHBPD1Pull","HB+ Depth=1 pull",100,-15,15)
    hHBPD2Pull = ROOT.TH1D("hHBPD2Pull","HB+ Depth=2 pull",100,-15,15)
    hHBPD3Pull = ROOT.TH1D("hHBPD3Pull","HB+ Depth=3 pull",100,-15,15)
    hHBPD4Pull = ROOT.TH1D("hHBPD4Pull","HB+ Depth=4 pull",100,-15,15)
    hHEMD1Pull = ROOT.TH1D("hHEMD1Pull","HE- Depth=1 pull",100,-15,15)
    hHEMD2Pull = ROOT.TH1D("hHEMD2Pull","HE- Depth=2 pull",100,-15,15)
    hHEMD3Pull = ROOT.TH1D("hHEMD3Pull","HE- Depth=3 pull",100,-15,15)
    hHEMD4Pull = ROOT.TH1D("hHEMD4Pull","HE- Depth=4 pull",100,-15,15)
    hHEMD5Pull = ROOT.TH1D("hHEMD5Pull","HE- Depth=5 pull",100,-15,15)
    hHEMD6Pull = ROOT.TH1D("hHEMD6Pull","HE- Depth=6 pull",100,-15,15)
    hHEMD7Pull = ROOT.TH1D("hHEMD7Pull","HE- Depth=7 pull",100,-15,15)
    hHEPD1Pull = ROOT.TH1D("hHEPD1Pull","HE+ Depth=1 pull",100,-15,15)
    hHEPD2Pull = ROOT.TH1D("hHEPD2Pull","HE+ Depth=2 pull",100,-15,15)
    hHEPD3Pull = ROOT.TH1D("hHEPD3Pull","HE+ Depth=3 pull",100,-15,15)
    hHEPD4Pull = ROOT.TH1D("hHEPD4Pull","HE+ Depth=4 pull",100,-15,15)
    hHEPD5Pull = ROOT.TH1D("hHEPD5Pull","HE+ Depth=5 pull",100,-15,15)
    hHEPD6Pull = ROOT.TH1D("hHEPD6Pull","HE+ Depth=6 pull",100,-15,15)
    hHEPD7Pull = ROOT.TH1D("hHEPD7Pull","HE+ Depth=7 pull",100,-15,15)
    hHFMD1Pull = ROOT.TH1D("hHFMD1Pull","HF- Depth=1 pull",100,-15,15)
    hHFMD2Pull = ROOT.TH1D("hHFMD2Pull","HF- Depth=2 pull",100,-15,15)
    hHFPD1Pull = ROOT.TH1D("hHFPD1Pull","HF+ Depth=1 pull",100,-15,15)
    hHFPD2Pull = ROOT.TH1D("hHFPD2Pull","HF+ Depth=2 pull",100,-15,15)

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


        for ix in range(1, h_ratio.GetNbinsX() + 1):
            for iy in range(1, h_ratio.GetNbinsY() + 1):
                content=h_ratio.GetBinContent(ix,iy)
                err=h_ratio.GetBinError(ix,iy)
                if err>0: err = (err**2+0.01**2)**.5
                if err>0: pull=(1-content)/err
                else:     pull=-999
                ieta = -42+ix
                iphi = iy
                if content!=0 and (content<.7 or content>1.5):
                    print("Outlier = ("+str(ieta)+", "+str(iphi)+", "+name+")="+str(content))
                if "depth1" in name and ieta>=-16 and ieta<=-1: hHBMD1Res.Fill(content); hHBMD1Pull.Fill(pull)
                if "depth1" in name and ieta<=16 and ieta>=1:   hHBPD1Res.Fill(content); hHBPD1Pull.Fill(pull)
                if "depth2" in name and ieta>=-16 and ieta<=-1: hHBMD2Res.Fill(content); hHBMD2Pull.Fill(pull)
                if "depth2" in name and ieta<=16 and ieta>=1:   hHBPD2Res.Fill(content); hHBPD2Pull.Fill(pull)
                if "depth3" in name and ieta>=-16 and ieta<=-1: hHBMD3Res.Fill(content); hHBMD3Pull.Fill(pull)
                if "depth3" in name and ieta<=16 and ieta>=1:   hHBPD3Res.Fill(content); hHBPD3Pull.Fill(pull)
                if "depth4" in name and ieta>=-15 and ieta<=-1: hHBMD4Res.Fill(content); hHBMD4Pull.Fill(pull)
                if "depth4" in name and ieta<=15 and ieta>=1:   hHBPD4Res.Fill(content); hHBPD4Pull.Fill(pull)
                if "depth1" in name and ieta>=-28 and ieta<=-17: hHEMD1Res.Fill(content); hHEMD1Pull.Fill(pull)
                if "depth1" in name and ieta<=28 and ieta>=17:   hHEPD1Res.Fill(content); hHEPD1Pull.Fill(pull)
                if "depth1" in name and ieta==-29 and iphi%2==0: hHEMD1Res.Fill(content); hHEMD1Pull.Fill(pull)
                if "depth1" in name and ieta==29 and iphi%2==0:  hHEPD1Res.Fill(content); hHEPD1Pull.Fill(pull)
                if "depth2" in name and ieta>=-28 and ieta<=-17: hHEMD2Res.Fill(content); hHEMD2Pull.Fill(pull)
                if "depth2" in name and ieta<=28 and ieta>=17:   hHEPD2Res.Fill(content); hHEPD2Pull.Fill(pull)
                if "depth2" in name and ieta==-29 and iphi%2==0: hHEMD2Res.Fill(content); hHEMD2Pull.Fill(pull)
                if "depth2" in name and ieta==29 and iphi%2==0:  hHEPD2Res.Fill(content); hHEPD2Pull.Fill(pull)
                if "depth3" in name and ieta>=-29 and ieta<=-17: hHEMD3Res.Fill(content); hHEMD3Pull.Fill(pull)
                if "depth3" in name and ieta<=29 and ieta>=17:   hHEPD3Res.Fill(content); hHEPD3Pull.Fill(pull)
                if "depth4" in name and ieta>=-29 and ieta<=-16: hHEMD4Res.Fill(content); hHEMD4Pull.Fill(pull)
                if "depth4" in name and ieta<=29 and ieta>=16:   hHEPD4Res.Fill(content); hHEPD4Pull.Fill(pull)
                if "depth5" in name and ieta>=-29 and ieta<=-17: hHEMD5Res.Fill(content); hHEMD5Pull.Fill(pull)
                if "depth5" in name and ieta<=29 and ieta>=17:   hHEPD5Res.Fill(content); hHEPD5Pull.Fill(pull)
                if "depth6" in name and ieta>=-29 and ieta<=-17: hHEMD6Res.Fill(content); hHEMD6Pull.Fill(pull)
                if "depth6" in name and ieta<=29 and ieta>=17:   hHEPD6Res.Fill(content); hHEPD6Pull.Fill(pull)
                if "depth7" in name and ieta>=-29 and ieta<=-17: hHEMD7Res.Fill(content); hHEMD7Pull.Fill(pull)
                if "depth7" in name and ieta<=29 and ieta>=17:   hHEPD7Res.Fill(content); hHEPD7Pull.Fill(pull)
                if "depth1" in name and ieta>=-41 and ieta<=-30: hHFMD1Res.Fill(content); hHFMD1Pull.Fill(pull)
                if "depth1" in name and ieta<=41 and ieta>=30:   hHFPD1Res.Fill(content); hHFPD1Pull.Fill(pull)
                if "depth1" in name and ieta==-29 and iphi%2==1: hHFMD1Res.Fill(content); hHFMD1Pull.Fill(pull)
                if "depth1" in name and ieta==29 and iphi%2==1:  hHFPD1Res.Fill(content); hHFPD1Pull.Fill(pull)
                if "depth2" in name and ieta>=-41 and ieta<=-30: hHFMD2Res.Fill(content); hHFMD2Pull.Fill(pull)
                if "depth2" in name and ieta<=41 and ieta>=30:   hHFPD2Res.Fill(content); hHFPD2Pull.Fill(pull)
                if "depth2" in name and ieta==-29 and iphi%2==1: hHFMD2Res.Fill(content); hHFMD2Pull.Fill(pull)
                if "depth2" in name and ieta==29 and iphi%2==1:  hHFMD2Res.Fill(content); hHFPD2Pull.Fill(pull)
                            

                    
        out.cd()
        h_ratio.SetMinimum(0.85)
        h_ratio.SetMaximum(1.15)
        h_ratio.Write()
    # end loop over keys


    hHBMD1Res.Write()
    hHBMD2Res.Write()
    hHBMD3Res.Write()
    hHBMD4Res.Write()
    hHBPD1Res.Write()
    hHBPD2Res.Write()
    hHBPD3Res.Write()
    hHBPD4Res.Write()
    hHEMD1Res.Write()
    hHEMD2Res.Write()
    hHEMD3Res.Write()
    hHEMD4Res.Write()
    hHEMD5Res.Write()
    hHEMD6Res.Write()
    hHEMD7Res.Write()
    hHEPD1Res.Write()
    hHEPD2Res.Write()
    hHEPD3Res.Write()
    hHEPD4Res.Write()
    hHEPD5Res.Write()
    hHEPD6Res.Write()
    hHEPD7Res.Write()
    hHFMD1Res.Write()
    hHFMD2Res.Write()
    hHFPD1Res.Write()
    hHFPD2Res.Write()

    hHBMD1Pull.Write()
    hHBMD2Pull.Write()
    hHBMD3Pull.Write()
    hHBMD4Pull.Write()
    hHBPD1Pull.Write()
    hHBPD2Pull.Write()
    hHBPD3Pull.Write()
    hHBPD4Pull.Write()
    hHEMD1Pull.Write()
    hHEMD2Pull.Write()
    hHEMD3Pull.Write()
    hHEMD4Pull.Write()
    hHEMD5Pull.Write()
    hHEMD6Pull.Write()
    hHEMD7Pull.Write()
    hHEPD1Pull.Write()
    hHEPD2Pull.Write()
    hHEPD3Pull.Write()
    hHEPD4Pull.Write()
    hHEPD5Pull.Write()
    hHEPD6Pull.Write()
    hHEPD7Pull.Write()
    hHFMD1Pull.Write()
    hHFMD2Pull.Write()
    hHFPD1Pull.Write()
    hHFPD2Pull.Write()

    out.Close()
    f1.Close()
    f2.Close()



if __name__ == "__main__":
    main()
