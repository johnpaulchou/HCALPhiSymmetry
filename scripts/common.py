import ROOT
import itertools
import sys

# this is code and parameters common to all of the scripts

# detector geometry and magic numbers
subdetnums = [ 1, 2, 4]
ndepths = [ 4, 7, 2]
minabsietas = [ 1, 16, 29]
maxabsietas = [ 16, 29, 41]
miniphi = 1
maxiphi = 72
modulus = 2


### returns a name for a subdetector based on the subdet number
def subdetName(subdetnum):
    if   subdetnum==subdetnums[0]: return "HB"
    elif subdetnum==subdetnums[1]: return "HE"
    elif subdetnum==subdetnums[2]: return "HF"
### end subdetName()


### returns the (minimum, maximum) energy used in the iterative method
def thresholds(subdetnum, ieta):
    if subdetnum==subdetnums[0]:   return (4., 100.)
    elif subdetnum==subdetnums[1]: return (4., 150.)
    elif subdetnum==subdetnums[2]:
        if abs(ieta)<=35:          return (10., 150.)
        else:                      return (15., 300.)
### end thresholds()


### returns the (nbins, minx, maxx) used in the histogram binning
def binning(subdetnum, ieta):

    # the lower bin needs to match the threshold in phisymtree.py
    if subdetnum==subdetnums[0]:   return (796, 1., 200.) 
    elif subdetnum==subdetnums[1]: return (796, 1., 200.)
    elif subdetnum==subdetnums[2]:
        if abs(ieta)<=35:          return (495, 5., 500.)
        else:                      return (198, 5., 500.)
### end binning()




### generate the name of the histogram for a given channel
def getHistName(subdetnum, ieta, iphi, depth, mod):

    if ieta<0: histname = subdetName(subdetnum) + "M_" + str(abs(ieta)) + "_" + str(iphi) + "_" + str(depth) + "_" + str(mod)
    else:      histname = subdetName(subdetnum) + "P_" + str(abs(ieta)) + "_" + str(iphi) + "_" + str(depth) + "_" + str(mod)
    return histname
### end getHistName()


        
### getHist() - opens a file and gets the histogram with the given detector ID
def getHist(rootfile, subdetnum, ieta, iphi, depth, mod, checkGoodChannel=True):

    # don't bother if it's not a good channel to begin with
    if checkGoodChannel and not goodChannel(subdetnum, ieta, iphi, depth, mod): return None

    histname = getHistName(subdetnum, ieta, iphi, depth, mod)
    
    # get the hist and return it
    rootfile.cd()
    hist = rootfile.Get(histname)
    if not hist or not isinstance(hist, ROOT.TH1) or hist is None:
        return None
    else:
        hist.SetDirectory(0)
        return hist
    ### end getHist()

### goodChannel() - determines if the channel should exist (or not)
def goodChannel(subdetnum, ieta, iphi, depth, mod):

    # do some sanity checks, first
    if iphi<=0 or iphi>=73: return False
    if ieta==0 or abs(ieta)>=42: return False
    if depth<=0 or depth>7: return False
    if mod>=modulus or mod<0: return False

    # HB checks
    if subdetnum==1 and depth>=1 and depth<=3 and abs(ieta)<=16: return True
    if subdetnum==1 and depth==4 and abs(ieta)<=15: return True

    # HE checks
    if subdetnum==2:
        if abs(ieta)==16 and depth==4: return True
        if abs(ieta)==17 and (depth==2 or depth==3): return True
        if abs(ieta)==18 and depth>=2 and depth<=5: return True
        if abs(ieta)>=19 and abs(ieta)<=20 and depth>=1 and depth<=6: return True
        if abs(ieta)>=21 and abs(ieta)<=25 and depth>=1 and depth<=6 and iphi%2==1: return True
        if abs(ieta)>=26 and abs(ieta)<=28 and depth>=1 and depth<=7 and iphi%2==1: return True
        if abs(ieta)==29 and depth>=1 and depth<=3 and iphi%2==1: return True

    # HF checks
    if subdetnum==4:
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

    for subdetindex,subdetnum in enumerate(subdetnums):
        for depth in range(1, ndepths[subdetindex]+1):
            ietas1 = range(minabsietas[subdetindex],maxabsietas[subdetindex]+1)
            ietas2 = range(-maxabsietas[subdetindex],-minabsietas[subdetindex]+1)
            ietas = itertools.chain(ietas1, ietas2)
            for ieta in ietas:
                for iphi in range(miniphi, maxiphi+1):
                    for mod in range(modulus):
                        val=goodChannel(subdetnum, ieta, iphi, depth, mod)
                        h=getHist(inputhistfile, subdetnum, ieta, iphi, depth, mod, False)
                        if h is None and val==True:
                            sys.stderr.write("No hist found for "+subdetName(subdetnum)+" ieta="+str(ieta)+" iphi="+str(iphi)+" depth="+str(depth)+" mod="+str(mod)+"\n")
                        elif h is not None and val==False:
                            sys.stderr.write("Hist found for "+subdetName(subdetnum)+" ieta="+str(ieta)+" iphi="+str(iphi)+" depth="+str(depth)+" mod="+str(mod)+"\n")
### end testGoodChannel()
