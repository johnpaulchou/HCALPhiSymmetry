#!/usr/bin/env python3

import ROOT
import itertools
import array as ar
import argparse

### getHist
def getHist(filename, subdet, ieta, iphi, depth):
    rootfile = ROOT.TFile(filename, "READ")
    if not rootfile or rootfile.IsZombie():
        print("Error: Unable to open file "+filename+".")
        exit(1)

    # select directory
    if subdet=="HB":   histname = "phaseHF/eHBspec/E_"
    elif subdet=="HE": histname = "phaseHF/eHEspec/E_"
    else:              histname = "phaseHF/espec/E_"

    # account for +/- in the name
    if ieta<0: histname = histname + "-" + str(abs(ieta)) + "_" + str(iphi) + "_" + str(depth)
    else:      histname = histname + "+" + str(abs(ieta)) + "_" + str(iphi) + "_" + str(depth)

    hist = rootfile.Get(histname)
    if not hist or not isinstance(hist, ROOT.TH1) or hist is None:
        rootfile.Close()
        return None
    else:
        hist.SetDirectory(0)
        rootfile.Close()
        return hist    
### end getHist()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("subdetector", help="select the subdetector (HB, HE, or HF)", choices = ["HB","HE","HF"])
    parser.add_argument("-i","--iterations", help="set the number of iterations",type=int, default=10)
    args=parser.parse_args()
    
    # detector geometry and magic numbers
    subdets = [ "HB", "HE", "HF" ]
    ndepths = [ 4, 7, 4]
    minabsietas = [ 1, 16, 29]
    maxabsietas = [ 16, 29, 41]
    miniphis = [ 1, 1, 1]
    maxiphis = [ 72, 72, 71]
    minthresholds = [ 4., 4., 10.]
    maxthresholds = [ 100., 150., 150.]
    if   args.subdetector==subdets[0]: subdetindex=0
    elif args.subdetector==subdets[1]: subdetindex=1
    elif args.subdetector==subdets[2]: subdetindex=2
    else: exit(1)
    subdet = subdets[subdetindex]
    niterations = args.iterations
    inputfilename = "EGamma0_Run2024C.root"
    outputtxtfilename = "corrs_"+subdet+".txt"
    outputhistfilename = "hists_"+subdet+".root"

    # open files
    outputtxtfile = open(outputtxtfilename, "w")
    outputhistfile = ROOT.TFile(outputhistfilename, "RECREATE")
    
    # create some temporary histograms
    hSmoothed = ROOT.TH1D("hSmoothed","Smoothed",10000,0,250)
    
    # loop over depths
    for depth in range(1, ndepths[subdetindex]+1):
        
        # loop over positive and negative ieta
        ietas1 = range(minabsietas[subdetindex],maxabsietas[subdetindex]+1)
        ietas2 = range(-maxabsietas[subdetindex],-minabsietas[subdetindex]+1)
        ietas = itertools.chain(ietas1, ietas2)
        for ieta in ietas:

            # loop over iphi the first time
            foundHist=False
            meanE=0.
            nmeanE=0
            hists = [] # store histograms so that you don't have to read it twice
            for iphi in range(miniphis[subdetindex], maxiphis[subdetindex]+1):
                
                # get the histogram and compute stuff
                h=getHist(inputfilename, subdet, ieta, iphi, depth)
                hists.append(h)
                if h is None: continue
                foundHist=True
                h.SetAxisRange(minthresholds[subdetindex],maxthresholds[subdetindex])
                rLP = h.Integral()*h.GetMean() # average of the E^2
                meanE=meanE+rLP
                nmeanE=nmeanE+1
            # end first iphi loop

            # skip this if no iphi slice was found
            if not foundHist: continue
            
            # compute the mean energy^2 across all iphi ranges
            meanE=meanE/nmeanE
            if meanE<=0:
                print("all histograms empty for subdet="+subdet+", ieta="+str(ieta)+", depth="+str(depth))
                continue
            
            # loop over iphi a second time
            for iphi in range(miniphis[subdetindex], maxiphis[subdetindex]+1):

                # get the histogram again
                h=hists[iphi-miniphis[subdetindex]]
                if h is None: continue
                h.SetAxisRange(minthresholds[subdetindex],maxthresholds[subdetindex])
                
                # create a spline
                xcorr = []
                ycorr = []
                for k in range(1, h.GetNbinsX()+1):
                    xcorr.append(h.GetBinCenter(k))
                    ycorr.append(h.GetBinContent(k))
                spline = ROOT.TSpline5("spline"+subdet+str(ieta)+str(iphi)+str(depth),ar.array('d',xcorr),ar.array('d',ycorr),len(xcorr),"",10,20)

                outputhistfile.cd()
                h.Write("hOriginal"+subdet+str(ieta)+str(iphi)+str(depth))
                
                # loop over iterations
                corr = 1.0
                for iter in range(niterations):

                    # fill a new histogram with the corrected energies based off the spline
                    hSmoothed.Reset()
                    for k in range(1,hSmoothed.GetNbinsX()+1):
                        x=hSmoothed.GetBinCenter(k)
                        hSmoothed.Fill(x*corr,spline.Eval(x)/10.0)
                    hSmoothed.SetAxisRange(minthresholds[subdetindex],maxthresholds[subdetindex])
                    rLP = hSmoothed.Integral()*hSmoothed.GetMean()
                    deltaLP = (rLP-meanE)/meanE

                    # maximum change is 70%
                    if abs(deltaLP)>0.7: deltaLP=0.7*deltaLP/abs(deltaLP)

                    # compute the uncertainty on the correction
                    if rLP>0:
                        drLP = ((hSmoothed.GetMeanError()/hSmoothed.GetMean())**2+
                                1./hSmoothed.Integral()+
                                (deltaLP/(1.+iter**.5))**2 )**.5
                    else: drLP=1e-6

                    # make the new correction (if we need to keep going)
                    if abs(deltaLP)>0.003 or iter<2:
                        corr = corr*(1-deltaLP/(1.+iter**.5))
                        correrr = corr*drLP
                    else:
                        correrr = corr*drLP
                        break
                    # end loop over iterations

                # print result
                outputtxtfile.write(subdet+" "+str(ieta)+" "+str(iphi)+" "+str(depth)+" "+str(corr)+" "+str(correrr)+"\n")
                    
                # end loop over iphi
            # end loop over ieta
        # end loop over depth

    outputtxtfile.close()
    outputhistfile.Close()
### end main function
