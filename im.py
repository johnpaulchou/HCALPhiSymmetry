#!/usr/bin/env python3

import ROOT
import itertools
import array as ar
import argparse
import sys
import scipy.optimize as spo
import scipy.integrate as integrate
import numpy as np
from scipy.stats import bootstrap
sys.stdout = open("output.txt", "w")
# detector geometry and magic numbers
subdets = [ "HB", "HE", "HF" ]
subdetnums = [ 1, 2, 4]
ndepths = [ 4, 7, 2]
minabsietas = [ 1, 16, 29]
maxabsietas = [ 16, 29, 41]
miniphi = 1
maxiphi = 72
minthresholds = [ 4., 4., 10.]
maxthresholds = [ 100., 150., 150.]
modulus = 20
doEFlow = True # if true, use the mean*integral, not the mean

### getHist() - opens a file and gets the histogram with the given detector ID
def getHist(rootfile, subdet, ieta, iphi, depth, mod, checkGoodChannel=True):

    # don't bother if it's not a good channel to begin with
    if checkGoodChannel and not goodChannel(subdet, ieta, iphi, depth, mod): return None

    # select directory
    if subdet=="HB":   histname = "phaseHF/eHBspec/E_"
    elif subdet=="HE": histname = "phaseHF/eHEspec/E_"
    else:              histname = "phaseHF/espec/E_"

    # account for +/- in the name
    if ieta<0: histname = histname + "-" + str(abs(ieta)) + "_" + str(iphi) + "_" + str(depth) + "_" + str(mod)
    else:      histname = histname + "+" + str(abs(ieta)) + "_" + str(iphi) + "_" + str(depth) + "_" + str(mod)

    # get the hist and return it
    rootfile.cd()
    hist = rootfile.Get(histname)
    if not hist or not isinstance(hist, ROOT.TH1) or hist is None:
        return None
    else:
        hist.SetDirectory(0)
        # REBIN THE HF histograms here
        if subdet=="HF": hist.Rebin(5)
        return hist
    ### end getHist()

### goodChannel() - determines if the channel should exist (or not)
def goodChannel(subdet, ieta, iphi, depth, mod):

    # do some sanity checks, first
    if iphi<=0 or iphi>=73: return False
    if ieta==0 or abs(ieta)>=42: return False
    if depth<=0 or depth>7: return False
    if mod>=modulus or mod<0: return False

    # HB checks
    if subdet=="HB" and depth>=1 and depth<=3 and abs(ieta)<=16: return True
    if subdet=="HB" and depth==4 and abs(ieta)<=15: return True

    # HE checks
    if subdet=="HE":
        if abs(ieta)==16 and depth==4: return True
        if abs(ieta)==17 and (depth==2 or depth==3): return True
        if abs(ieta)==18 and depth>=2 and depth<=5: return True
        if abs(ieta)>=19 and abs(ieta)<=20 and depth>=1 and depth<=6: return True
        if abs(ieta)>=21 and abs(ieta)<=25 and depth>=1 and depth<=6 and iphi%2==1: return True
        if abs(ieta)>=26 and abs(ieta)<=28 and depth>=1 and depth<=7 and iphi%2==1: return True
        if abs(ieta)==29 and depth>=1 and depth<=3 and iphi%2==1: return True

    # HF checks
    if subdet=="HF":
        if depth<1 or depth >2: return False
        if iphi%2==0: return False
        if abs(ieta)>=29 and abs(ieta)<=39: return True
        if abs(ieta)>=40 and abs(ieta)<=41 and iphi%4==3: return True

    return False
    ### end goodChannel()

# testGoodChannel() - tests whether the goodChannel is consistent with the file
def testGoodChannel(inputfilename):
    inputhistfile = ROOT.TFile(inputfilename, "READ")
    if not inputhistfile or inputhistfile.IsZombie():
        print("Error: Unable to open file "+inputfilename+".")
        exit(1)

    for subdetindex,subdet in enumerate(subdets):
        for depth in range(1, ndepths[subdetindex]+1):
            ietas1 = range(minabsietas[subdetindex],maxabsietas[subdetindex]+1)
            ietas2 = range(-maxabsietas[subdetindex],-minabsietas[subdetindex]+1)
            ietas = itertools.chain(ietas1, ietas2)
            for ieta in ietas:
                for iphi in range(miniphi, maxiphi+1):
                    for mod in range(modulus):
                        val=goodChannel(subdet, ieta, iphi, depth, mod)
                        h=getHist(inputhistfile, subdet, ieta, iphi, depth, mod, False)
                        if h is None and val==True:
                            outputerrfile.write("No hist found for "+subdet+" "+str(ieta)+" "+str(iphi)+" "+str(depth)+" "+str(mod))
                        elif h is not None and val==False:
                            outputerrfile.write("Hist found for "+subdet+" "+str(ieta)+" "+str(iphi)+" "+str(depth)+" "+str(mod))
    ### end testGoodChannel()

# defines the function to minimize
class minimizeFunc:
    def __init__(self, spline, mean, lo, hi):
        self.spline = spline
        self.mean = mean
        self.lo = lo
        self.hi = hi

    def splineF(self, x):
        return self.spline.Eval(x)

    def splineMeanF(self, x, corr):
        return x*corr*self.spline.Eval(x)

    def minimizer(self,corr):
        # integrate from the new bounds
        splineIntegral = integrate.quad(self.splineF,self.lo/corr[0],self.hi/corr[0])
        splineMean = integrate.quad(self.splineMeanF,self.lo/corr[0],self.hi/corr[0],args=(corr[0]))

        # compute the mean and return the difference with the target mean squared
        if doEFlow: return (splineMean[0] - self.mean)**2
        else: return (splineMean[0]/splineIntegral[0] - self.mean)**2

    def minimize(self):
        result = spo.minimize(self.minimizer, 1.0, bounds=[(0.10,3.)],tol=1e-4,method="Nelder-Mead")
        if result.success: return result.x[0]
        else: return -1.

### define the main function
def main():
    # setup parser
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputfilename", help="name of .root file containing the histograms to be processed")
    parser.add_argument("--subdetector", help="select the subdetector (HB, HE, HF, or all)", default="all", choices=["HB","HE","HF","all"])
    parser.add_argument("--side",help="side of the detector (+, -, or all)", default="all", choices=["+","-","all"])
    parser.add_argument("--ieta",help="ieta value (0 for all ietas)", type=int, default=0)
    parser.add_argument("--depth",help="depth value (0 for all depths)", type=int, default=0)
    args=parser.parse_args()

    # filenames
    inputfilename = args.inputfilename
    if args.ieta==0:
        if args.side=="+":   filestr = args.subdetector+"P"
        elif args.side=="-": filestr = args.subdetector+"M"
        else:                filestr = args.subdetector+"PM"
    else:
        filestr = args.subdetector+"_e"+str(args.ieta)
    if args.depth!=0:
        filestr = filestr + "_d"+str(args.depth)
    outputtxtfilename = "corrs_"+filestr+".txt"
    outputerrfilename = "err_"+filestr+".txt"
    outputhistfilename = "hists_"+filestr+".root"

    # open files
    outputtxtfile = open(outputtxtfilename, "w")
    outputerrfile = open(outputerrfilename, "w")
    outputhistfile = ROOT.TFile(outputhistfilename, "RECREATE")
    inputhistfile = ROOT.TFile(inputfilename, "READ")
    if not inputhistfile or inputhistfile.IsZombie():
        print("Error: Unable to open file "+inputfilename+".")
        exit(1)

    # loop over subdetectors
    for subdetindex,subdet in enumerate(subdets):
        if subdet!=args.subdetector and args.subdetector!="all": continue

        # loop over depths
        for depth in range(1, ndepths[subdetindex]+1):
            if depth!=args.depth and args.depth!=0: continue

            # loop over positive and negative ieta
            ietas1 = range(minabsietas[subdetindex],maxabsietas[subdetindex]+1)
            ietas2 = range(-maxabsietas[subdetindex],-minabsietas[subdetindex]+1)
            ietas = itertools.chain(ietas1, ietas2)
            for ieta in ietas:
                if args.side!="all" and ((ieta<0 and args.side!="-") or (ieta>0 and args.side!="+")): continue
                if ieta!=args.ieta and args.ieta!=0: continue

                # store the corrections here
                corrs = [[-1 for mod in range(modulus)] for iphi in range(maxiphi+1-miniphi)]

                # loop over moduli
                for mod in range(modulus):

                    # loop over iphi the first time
                    foundHist=False
                    meanE=0.
                    nmeanE=0
                    hists = [None]*(maxiphi-miniphi+1) # store histograms so that you don't have to read it twice
                    for iphi in range(miniphi, maxiphi+1):

                        # get the histogram and compute stuff
                        h=getHist(inputhistfile, subdet, ieta, iphi, depth, mod)

                        # if the histogram has issues, change it to "None" to be skipped later, and print a warning
                        if h is not None:
                            h.SetAxisRange(minthresholds[subdetindex],maxthresholds[subdetindex])
                            if h.Integral()<=0 or h.GetRMS()<=0:
                                print("Histogram Warning: "+subdet+", ieta="+str(ieta)+", iphi="+str(iphi)+", depth="+str(depth)+", mod="+str(mod)+" is empty or has some issues",file=sys.stderr)
                                h=None

                        # store it in the list, even if it is bad
                        hists[iphi-miniphi]=h

                        # skip bad histograms
                        if h is None: continue
                        foundHist=True
                        h.SetAxisRange(minthresholds[subdetindex],maxthresholds[subdetindex])
                        if doEFlow: meanE=meanE+h.GetMean()*h.Integral("width") # average energy*integral (normalized to bin size)
                        else:       meanE=meanE+h.GetMean() # average energy
                        nmeanE=nmeanE+1
                        # end first iphi loop

                    # skip this if no iphi slice was found
                    if not foundHist: continue

                    # compute the average across all iphi ranges
                    meanE=meanE/nmeanE
                    if meanE<=0:
                        print("all histograms empty for subdet="+subdet+", ieta="+str(ieta)+", depth="+str(depth)+", mod="+str(mod),file=sys.stderr)
                        continue

                    # loop over iphi a second time
                    for iphi in range(miniphi, maxiphi+1):

                        # get the histogram again
                        h=hists[iphi-miniphi]
                        if h is None: continue
                        h.SetAxisRange(minthresholds[subdetindex],maxthresholds[subdetindex])

                        # create a spline
                        xcorr = [0.0]*h.GetNbinsX()
                        ycorr = [0.0]*h.GetNbinsX()
                        for k in range(1, h.GetNbinsX()+1):
                            xcorr[k-1]=h.GetBinCenter(k)
                            ycorr[k-1]=h.GetBinContent(k)
                        spline = ROOT.TSpline5("spline"+subdet+str(ieta)+str(iphi)+str(depth),ar.array('d',xcorr),ar.array('d',ycorr),len(xcorr),"",10,20)

                        outputhistfile.cd()
                        if ieta<0:
                            h.Write("hOriginal"+subdet+"M"+str(abs(ieta))+"_"+str(iphi)+"_"+str(depth)+"_"+str(mod))
                            spline.Write("spline"+subdet+"M"+str(abs(ieta))+"_"+str(iphi)+"_"+str(depth)+"_"+str(mod))
                        else:
                            h.Write("hOriginal"+subdet+"P"+str(abs(ieta))+"_"+str(iphi)+"_"+str(depth)+"_"+str(mod))
                            spline.Write("spline"+subdet+"P"+str(abs(ieta))+"_"+str(iphi)+"_"+str(depth)+"_"+str(mod))

                        # compute the correction and store it
                        mini = minimizeFunc(spline, meanE, minthresholds[subdetindex], maxthresholds[subdetindex])
                        corr=mini.minimize()
                        if corr<0:
                            print("Convergence Failure Warning for: "+subdet+", ieta="+str(ieta)+", iphi="+str(iphi)+", depth="+str(depth)+", mod="+str(mod),file=sys.stderr)
                        corrs[iphi-miniphi][mod]=corr

                        # end loop over iphi
                    # end loop over moduli

                # last loop over iphi
                for iphi in range(miniphi, maxiphi+1):

                    # compute average and stddev of the correction over moduli
                    data = np.array(corrs[iphi-miniphi])
                    data = data[data >= 0] # eliminate negative data values
                    #print(iphi)
                    #print(data)
                    if len(data)<modulus:
                        outputerrfile.write("Could not find all moduli for " +str(ieta)+" "+str(iphi)+" "+str(depth)+ " "+str(subdetnums[subdetindex]))
                        #corrStr = -3
                        #corrErrStr = -0.00001
                        #outputerrfile.write(str(subdetnums[subdetindex])+" "+str(ieta)+" "+str(iphi)+" "+str(depth)+"\n")
                        continue
                    avgcorr=np.mean(data)
#                    stddevcorr=np.std(data, ddof=1, mean=avgcorr)/len(data)**.5
                    d=(data,)
                    stddevcorr=bootstrap(d, np.mean, confidence_level=0.68,n_resamples=999).standard_error

                    # write the results to the file
                    corrStr = f"{avgcorr:.5f}"
                    corrErrStr = f"{stddevcorr:.5f}"
                    outputtxtfile.write(str(subdetnums[subdetindex])+" "+str(ieta)+" "+str(iphi)+" "+str(depth)+" "+corrStr+" "+corrErrStr+"\n")

                    # end loop over iphi
                # end loop over ieta
            # end loop over depth
        # end loop over subdets

    outputtxtfile.close()
    outputerrfile.close()
    outputhistfile.cd()
    inputhistfile.Close()
    outputhistfile.Close()
### end main function


if __name__ == "__main__":
    main()
