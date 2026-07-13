#!/usr/bin/env python3

import ROOT
import itertools
import argparse
import sys
import numpy as np
import common
from scipy.stats import bootstrap


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

    # create the target dictionary here
    targets = {}
        
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

                # compute the targets
                target=common.computeTargetEflow(inputhistfile, subdet, ieta, depth)
                if target<0: continue
                targets[subdet, ieta, depth]=target

                # loop over iphi
                for iphi in range(common.miniphi, common.maxiphi+1):

                    # loop over moduli and compute the correction
                    corrs = []
                    for mod in range(common.modulus):
                        corr=common.computeCorrection(inputhistfile, subdet, ieta, iphi, depth, mod, target)
                        if corr<0 or corr>=3.:
                            if common.goodChannel(subdet, ieta, iphi, depth, mod):
                                outputerrfile.write("Convergence Failure Warning for: "+str(subdet)+", ieta="+str(ieta)+", iphi="+str(iphi)+", depth="+str(depth)+", mod="+str(mod)+"\n")
                            continue
                        corrs.append(corr)

                    # compute average and stddev of the correction over moduli
                    data = np.array(corrs)
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
    inputhistfile.Close()
### end main function


if __name__ == "__main__":
    main()
