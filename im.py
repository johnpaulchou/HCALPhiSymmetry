#!/usr/bin/env python3

import ROOT
import itertools
import array as ar
import argparse

### getHist
def getHist(filename, subdet, ieta, iphi, depth):

    # don't bother if it's not a good channel to begin with
    if not goodChannel(subdet, ieta, iphi, depth): return None

    # try to open the file
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

    # get the hist and return it
    hist = rootfile.Get(histname)
    if not hist or not isinstance(hist, ROOT.TH1) or hist is None:
        rootfile.Close()
        return None
    else:
        hist.SetDirectory(0)
        rootfile.Close()
        return hist    
    ### end getHist()

### goodChannel - determines if the channel should exist (or not)
def goodChannel(subdet, ieta, iphi, depth):

    # do some sanity checks, first
    if iphi<=0 or iphi>=73: return False
    if ieta==0 or abs(ieta)>=42: return False
    if depth<=0 or depth>7: return False

    # HB checks
    if subdet=="HB" and (depth==1 or depth==2 or depth==3) and abs(ieta)<=16: return True
    if subdet=="HB" and depth==4 and abs(ieta)<=15: return True

    # HE checks
    if subdet=="HE":
        if abs(ieta)==16 and depth==4: return True
        if abs(ieta)==17 and (depth==2 or depth==3): return True
        if abs(ieta)==18 and depth>=1 and depth<=5: return True
        if abs(ieta)>=19 and abs(ieta)<=20 and depth>=1 and depth<=6: return True
        if abs(ieta)>=21 and abs(ieta)<=25 and depth>=1 and depth<=6 and iphi%2==1: return True
        if abs(ieta)>=26 and abs(ieta)<=28 and depth>=1 and depth<=7 and iphi%2==1: return True
        if abs(ieta)==29 and depth>=1 and depth<=3 and iphi%2==1: return True

    # HF checks
    if subdet=="HF":
        if iphi%2==0: return False
        if abs(ieta)>=29 and abs(ieta)<=39 and depth>=1 and depth<=4: return True
        if abs(ieta)>=40 and abs(ieta)<=41 and iphi%4==3: return True
    
    return False
    ### end goodChannel()

def testGoodChannel():

    inputfilename = "EGamma0_Run2024C.root"
    subdets = [ "HB", "HE", "HF" ]
    ndepths = [ 4, 7, 2]
    minabsietas = [ 1, 16, 29]
    maxabsietas = [ 16, 29, 41]
    miniphis = [ 1, 1, 1]
    maxiphis = [ 72, 72, 71]

    for subdetindex,subdet in enumerate(subdets):
        for depth in range(1, ndepths[subdetindex]+1):
            ietas1 = range(minabsietas[subdetindex],maxabsietas[subdetindex]+1)
            ietas2 = range(-maxabsietas[subdetindex],-minabsietas[subdetindex]+1)
            ietas = itertools.chain(ietas1, ietas2)
            for ieta in ietas:
                for iphi in range(miniphis[subdetindex], maxiphis[subdetindex]+1):
                    val=goodChannel(subdet, ieta, iphi, depth)
                    h=getHist(inputfilename, subdet, ieta, iphi, depth)
                    if h is None and val==True:
                        print("No hist found for "+subdet+" "+str(ieta)+" "+str(iphi)+" "+str(depth))
                    elif h is not None and val==False:
                        print("Hist found for "+subdet+" "+str(ieta)+" "+str(iphi)+" "+str(depth))




if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("inputfilename", help="name of .root file containing the histograms to be processed")
    parser.add_argument("subdetector", help="select the subdetector (HB, HE, HF, or all)", choices = ["HB","HE","HF","all"])
    parser.add_argument("-i","--iterations", help="set the number of iterations",type=int, default=10)
    args=parser.parse_args()
    
    # detector geometry and magic numbers
    subdets = [ "HB", "HE", "HF" ]
    subdetnums = [ 1, 2, 4]
    ndepths = [ 4, 7, 4]
    minabsietas = [ 1, 16, 29]
    maxabsietas = [ 16, 29, 41]
    miniphi = 1
    maxiphi = 72
    minthresholds = [ 4., 4., 10.]
    maxthresholds = [ 100., 150., 150.]

    niterations = args.iterations
    inputfilename = args.inputfilename
    outputtxtfilename = "corrs_"+args.subdetector+".txt"
    outputhistfilename = "hists_"+args.subdetector+".root"

    # open files
    outputtxtfile = open(outputtxtfilename, "w")
    outputhistfile = ROOT.TFile(outputhistfilename, "RECREATE")
    
    # create some temporary histograms
    hSmoothed = ROOT.TH1D("hSmoothed","Smoothed",10000,0,250)

    # loop over subdetectors
    for subdetindex,subdet in enumerate(subdets):
        if subdet!=args.subdetector and args.subdetector!="all": continue
    
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
                for iphi in range(miniphi, maxiphi+1):
                
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
                for iphi in range(miniphi, maxiphi+1):

                    # get the histogram again
                    h=hists[iphi-miniphi]
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
                    outputtxtfile.write(str(subdetnums[subdetindex])+" "+str(ieta)+" "+str(iphi)+" "+str(depth)+" "+str(corr)+" "+str(correrr)+"\n")
                    
                    # end loop over iphi
                # end loop over ieta
            # end loop over depth
        # end loop over subdets

    outputtxtfile.close()
    outputhistfile.Close()
### end main function
