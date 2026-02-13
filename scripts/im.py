#!/usr/bin/env python3

import ROOT
import itertools
import array as ar
import argparse
import sys
import scipy.optimize as spo
import scipy.integrate as integrate
import numpy as np
import common
from scipy.stats import bootstrap
from ctypes import *
#sys.stdout = open("output.txt", "w")

doEFlow = True # if true, use the mean*integral, not the integral


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

    def splineIntegrator(self,x0,x,y,b,c,d,e,f):
        return f/6*(x-x0)**6+e/5*(x-x0)**5+d/4*(x-x0)**4+c/3*(x-x0)**3+b/2*(x-x0)**2+y*(x-x0)
    
    def splineMeanIntegrator(self,corr,x0,x,y,b,c,d,e,f):
        return (f/7*corr*(x-x0)**7+f/6*corr*x0*(x-x0)**6+e/6*corr*(x-x0)**6+e/5*x0*corr*(x-x0)**5
                +d/5*corr*(x-x0)**5+d/4*x0*corr*(x-x0)**4+c/4*corr*(x-x0)**4+c/3*x0*corr*(x-x0)**3
                +b/3*corr*(x-x0)**3+b/2*x0*corr*(x-x0)**2+y/2*corr*(x**2-x0**2))

    def minimizer(self,corr):
        # integrate from the new bounds
        #numerical integration (legacy)
        #splineIntegral = integrate.quad(self.splineF,self.lo/corr[0],self.hi/corr[0],limit=1000)
        #splineMean = integrate.quad(self.splineMeanF,self.lo/corr[0],self.hi/corr[0],args=(corr[0]),limit=1000)
        
        #analytical integration (new)
        splineMean2= 0.0 
        splineIntegral2 = 0.0
        
        for i in range(self.spline.GetNp()-1):
            x0val=c_double(0.0)
            y=c_double(0.0)
            b=c_double(0.0)
            c=c_double(0.0)
            d=c_double(0.0)
            e=c_double(0.0)
            f=c_double(0.0)
            
            self.spline.GetCoeff(i,x0val,y,b,c,d,e,f)
            if(x0val.value>(self.hi/corr[0])):break
            x1val=c_double(0.0)
            y1=c_double(0.0)
            self.spline.GetKnot(i+1,x1val,y1)
            if(x1val.value<(self.lo/corr[0])):continue

            x1 = min(x1val.value,self.hi/corr[0])
            
            if(x0val.value<self.lo/corr[0]):
                splineIntegral2+=(self.splineIntegrator(x0val.value,x1,y.value,b.value,c.value,d.value,e.value,f.value)
                                  -self.splineIntegrator(x0val.value,self.lo/corr[0],y.value,b.value,c.value,d.value,e.value,f.value))
                splineMean2+=(self.splineMeanIntegrator(corr[0],x0val.value,x1,y.value,b.value,c.value,d.value,e.value,f.value)
                              -self.splineMeanIntegrator(corr[0],x0val.value,self.lo/corr[0],y.value,b.value,c.value,d.value,e.value,f.value))
            else: 
                splineIntegral2+=self.splineIntegrator(x0val.value,x1,y.value,b.value,c.value,d.value,e.value,f.value)
                splineMean2+=self.splineMeanIntegrator(corr[0],x0val.value,x1,y.value,b.value,c.value,d.value,e.value,f.value)
        
        # compute the mean and return the difference with the target mean squared
        if doEFlow: return (splineMean2 - self.mean)**2
        else: return (splineMean2/splineIntegral2 - self.mean)**2
        #else: return (splineIntegral[0] - self.mean)**2 #alt method

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

    # open files
    outputtxtfile = open(outputtxtfilename, "w")
    outputerrfile = open(outputerrfilename, "w")
    inputhistfile = ROOT.TFile(inputfilename, "READ")
    if not inputhistfile or inputhistfile.IsZombie():
        print("Error: Unable to open file "+inputfilename+".")
        exit(1)

    # loop over subdetectors
    for subdetindex,subdet in enumerate(common.subdetnums):
        if common.subdetName(subdet)!=args.subdetector and args.subdetector!="all": continue

        # loop over depths
        for depth in range(1, common.ndepths[subdetindex]+1):
            if depth!=args.depth and args.depth!=0: continue

            # loop over positive and negative ieta
            ietas1 = range(common.minabsietas[subdetindex],common.maxabsietas[subdetindex]+1)
            ietas2 = range(-common.maxabsietas[subdetindex],-common.minabsietas[subdetindex]+1)
            ietas = itertools.chain(ietas1, ietas2)
            for ieta in ietas:
                if args.side!="all" and ((ieta<0 and args.side!="-") or (ieta>0 and args.side!="+")): continue
                if ieta!=args.ieta and args.ieta!=0: continue

                # store the corrections here
                corrs = [[-1 for mod in range(common.modulus)] for iphi in range(common.maxiphi+1-common.miniphi)]

                # loop over moduli
                for mod in range(common.modulus):

                    # loop over iphi the first time
                    foundHist=False
                    meanE=0.
                    nmeanE=0
                    hists = [None]*(common.maxiphi-common.miniphi+1) # store histograms so that you don't have to read it twice
                    for iphi in range(common.miniphi, common.maxiphi+1):

                        # get the histogram and compute stuff
                        h=common.getHist(inputhistfile, subdet, ieta, iphi, depth, mod)

                        # if the histogram has issues, change it to "None" to be skipped later, and print a warning
                        if h is not None:
                            (minthresh,maxthresh)=common.thresholds(subdet,ieta)
                            h.SetAxisRange(minthresh,maxthresh)
                            if h.Integral()<=0 or h.GetRMS()<=0:
                                outputerrfile.write("Histogram Warning: "+str(subdet)+", ieta="+str(ieta)+", iphi="+str(iphi)+", depth="+str(depth)+", mod="+str(mod)+" is empty or has some issues\n")
                                h=None

                        # store it in the list, even if it is bad
                        hists[iphi-common.miniphi]=h

                        # skip bad histograms
                        if h is None: continue
                        foundHist=True
                        (minthresh,maxthresh)=common.thresholds(subdet,ieta)
                        h.SetAxisRange(minthresh,maxthresh)
                        if doEFlow: meanE=meanE+h.GetMean()*h.Integral("width") # average energy*integral (normalized to bin size)
                        else:       meanE=meanE+h.GetMean()
                        #else: meanE=meanE+h.Integral("width") #alt method
                        nmeanE=nmeanE+1
                        # end first iphi loop

                    # skip this if no iphi slice was found
                    if not foundHist: continue

                    # compute the average across all iphi ranges
                    meanE=meanE/nmeanE
                    if meanE<=0:
                        outputerrfile.write("all histograms empty for subdet="+subdet+", ieta="+str(ieta)+", depth="+str(depth)+", mod="+str(mod)+"\n")
                        continue

                    # loop over iphi a second time
                    for iphi in range(common.miniphi, common.maxiphi+1):

                        # get the histogram again
                        h=hists[iphi-common.miniphi]
                        if h is None: continue
                        (minthresh,maxthresh)=common.thresholds(subdet,ieta)
                        h.SetAxisRange(minthresh,maxthresh)

                        # create a spline
                        xcorr = [0.0]*h.GetNbinsX()
                        ycorr = [0.0]*h.GetNbinsX()
                        for k in range(1, h.GetNbinsX()+1):
                            xcorr[k-1]=h.GetBinCenter(k)
                            ycorr[k-1]=h.GetBinContent(k)
                        spline = ROOT.TSpline5("spline"+str(subdet)+str(ieta)+str(iphi)+str(depth),ar.array('d',xcorr),ar.array('d',ycorr),len(xcorr),"",10,20)

                        # compute the correction and store it
                        mini = minimizeFunc(spline, meanE, minthresh, maxthresh)
                        corr=mini.minimize()
                        if corr<0:
                            outputerrfile.write("Convergence Failure Warning for: "+subdet+", ieta="+str(ieta)+", iphi="+str(iphi)+", depth="+str(depth)+", mod="+str(mod)+"\n")
                        corrs[iphi-common.miniphi][mod]=corr

                        # end loop over iphi
                    # end loop over moduli

                # last loop over iphi
                for iphi in range(common.miniphi, common.maxiphi+1):

                    # compute average and stddev of the correction over moduli
                    data = np.array(corrs[iphi-common.miniphi])
                    data = data[data >= 0] # eliminate negative data values
                    if len(data)<common.modulus: continue # skip if we didn't have enough data
                    avgcorr=np.mean(data)
                    d=(data,)
                    stddevcorr=bootstrap(d, np.mean, confidence_level=0.68,n_resamples=999).standard_error

                    # write the results to the file
                    corrStr = f"{avgcorr:.5f}"
                    corrErrStr = f"{stddevcorr:.5f}"
                    outputtxtfile.write(str(common.subdetnums[subdetindex])+" "+str(ieta)+" "+str(iphi)+" "+str(depth)+" "+corrStr+" "+corrErrStr+"\n")

                    # end loop over iphi
                # end loop over ieta
            # end loop over depth
        # end loop over subdets

    outputtxtfile.close()
    outputerrfile.close()
    inputhistfile.Close()
### end main function


if __name__ == "__main__":
    main()
