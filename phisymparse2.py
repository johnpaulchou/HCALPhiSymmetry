import ROOT
import argparse
import math
import pandas as pd
#import numpy as np
#import matplotlib.pyplot as plt
#from scipy.stats import norm
from itertools import chain
ROOT.gROOT.SetBatch(True)
import os
import sys
import re

def str_to_float_list(s):
    return [abs(float(th)) for th in s.split(',')]

def strip_float(fl): #for image file names
    return str(fl).replace('.','-').rstrip('0')

def get_df(filename): #for calibration method outputs
    df=pd.read_csv(filename,sep='\s+',header=None)
    df=df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    df.columns = ['Subdetector', 'iEta', 'iPhi', 'Depth', 'Correction', 'Error']
    filt= df[(df['Depth'] == 1) & (df['Subdetector'] == 1)& (df['iEta']<=-1 )& (df['iEta']>=-16 ) &((df['iPhi']!=16)&(df['iEta']!=-16))
                      &((df['iPhi']!=44)&(df['iEta']!=-8)) &((df['iPhi']!=10)&(df['iEta']!=-2))]
    filt = filt.reset_index(drop=True)
    return df

def get_df2(filename): #for outliers
    df=pd.read_csv(filename,sep='\s+',header=None)
    df=df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    df.columns = ['Subdetector', 'iEta', 'iPhi', 'Depth', 'Correction1', 'Error1','Correction2','Error2','Pull','Residual']
    filt= df[(df['Depth'] == 1) & (df['Subdetector'] == 1)& (df['iEta']<=-1 )& (df['iEta']>=-16 ) &((df['iPhi']!=16)&(df['iEta']!=-16))
                      &((df['iPhi']!=44)&(df['iEta']!=-8)) &((df['iPhi']!=10)&(df['iEta']!=-2))]
    filt = filt.reset_index(drop=True)
    return df

def get_filename(filename):
    return os.path.basename(filename)[:filename.rindex('.')] #cut off the directory

def get_file_extension(filename):
    return filename.split('.')[-1]

#def get_corr_err_uncert(filename,subdetector,ieta,iphi,depth):
    #df=get_df(filename)
    #corr=df.loc[(df['Subdetector']==subdetector) & (df['iEta']==ieta) & (df['iPhi']==iphi) & (df['Depth']==depth,'Correction')]
    #corr=df.query('Subdetector==@subdetector and iEta==@ieta and iPhi==@iphi and Depth==@depth')['Correction'].iloc[0]
    #corre=df.query('Subdetector==@subdetector and iEta==@ieta and iPhi==@iphi and Depth==@depth')['Error'].iloc[0]
    #return (corr,corre,corre/corr)

def get_pull(df1,df2,subdetector,ieta,iphi,depth):
    corr1=df1.query('Subdetector==@subdetector and iEta==@ieta and iPhi==@iphi and Depth==@depth')['Correction'].iloc[0]
    corre1=df1.query('Subdetector==@subdetector and iEta==@ieta and iPhi==@iphi and Depth==@depth')['Error'].iloc[0]
    corr2=df2.query('Subdetector==@subdetector and iEta==@ieta and iPhi==@iphi and Depth==@depth')['Correction'].iloc[0]
    corre2=df2.query('Subdetector==@subdetector and iEta==@ieta and iPhi==@iphi and Depth==@depth')['Error'].iloc[0]
    #(corr1,corre1,uncert1)=get_corr_err_uncert(file1,subdetector,ieta,iphi,depth)
    #(corr2,corre2,uncert2)=get_corr_err_uncert(file2,subdetector,ieta,iphi,depth)
    return (corr1-corr2)/math.sqrt(math.pow(corre1,2.0)+math.pow(corre2,2.0))

def get_resid(df1,df2,subdetector,ieta,iphi,depth):
    corr1=df1.query('Subdetector==@subdetector and iEta==@ieta and iPhi==@iphi and Depth==@depth')['Correction'].iloc[0]
    corr2=df2.query('Subdetector==@subdetector and iEta==@ieta and iPhi==@iphi and Depth==@depth')['Correction'].iloc[0]
    return 1.0-corr2/corr1

def corr_plots(df,filename,depth):
    corrs=ROOT.TH2D('corrs_f_'+filename+'_d'+str(depth),'Corrections for Depth '+str(depth),83,-41,41,73,0,72)
    for row in df.itertuples():
        corr=row.Correction
        if (row.Depth!=depth):continue
        if(row.Subdetector==4 and (row.iEta==29 or row.iEta==-29)): corrs.SetBinContent(row.iEta+42,row.iPhi+1,corr)
        else:corrs.SetBinContent(row.iEta+42,row.iPhi,corr)
    outroot.WriteObject(corrs,corrs.GetName())

def correlation_plot(df1,df2,rem_ones):
    graphhb = ROOT.TGraph()
    if rem_ones: graphhb.SetName('correlation_f1_'+filename1+"_f2_"+filename2+'_hb_ones_removed')
    else: graphhb.SetName('correlation_f1_'+filename1+"_f2_"+filename2+'_hb')
    graphhb.SetTitle("Correlation for Barrel")
    graphhb.SetMarkerStyle(20)
    graphhb.SetMarkerSize(0.5)
    graphhb.SetMarkerColor(ROOT.kBlue)
    graphhb.GetXaxis().SetTitle(filename1)  # Set X-axis title
    graphhb.GetYaxis().SetTitle(filename2)
    graphhb.GetXaxis().SetLimits(0, 12)
    graphhe = ROOT.TGraph()
    if rem_ones: graphhe.SetName('correlation_f1_'+filename1+"_f2_"+filename2+'_he_ones_removed')
    else: graphhe.SetName('correlation_f1_'+filename1+"_f2_"+filename2+'_he')
    graphhe.SetTitle("Correlation for Endcap")
    graphhe.SetMarkerStyle(20)
    graphhe.SetMarkerSize(0.5)
    graphhe.SetMarkerColor(ROOT.kBlue)
    graphhe.GetXaxis().SetTitle(filename1)  # Set X-axis title
    graphhe.GetYaxis().SetTitle(filename2)
    graphhf = ROOT.TGraph()
    if rem_ones: graphhf.SetName('correlation_f1_'+filename1+"_f2_"+filename2+'_hf_ones_removed')
    else: graphhf.SetName('correlation_f1_'+filename1+"_f2_"+filename2+'_hf')
    graphhf.SetTitle("Correlation for Forward")
    graphhf.SetMarkerStyle(20)
    graphhf.SetMarkerSize(0.5)
    graphhf.SetMarkerColor(ROOT.kBlue)
    graphhf.GetXaxis().SetTitle(filename1)  # Set X-axis title
    graphhf.GetYaxis().SetTitle(filename2)
    #x_minb, x_maxb = float("inf"), float("-inf")
    #x_mine, x_maxe = float("inf"), float("-inf")
    #x_minf, x_maxf = float("inf"), float("-inf")
    for row in df1.itertuples():
        (subdet,ieta,iphi,depth)=(row.Subdetector,row.iEta,row.iPhi,row.Depth)
        corr1=row.Correction
        if(((df2['Subdetector']==subdet) & (df2['iEta']==ieta) & (df2['iPhi']==iphi) & (df2['Depth']==depth)).any()):
            corr2=df2.query('Subdetector==@subdet and iEta==@ieta and iPhi==@iphi and Depth==@depth')['Correction'].iloc[0]
            if(corr2!=0):
                if (rem_ones==True and (corr2==1 or corr1==1)):
                    continue
                if (subdet==4) :
                    graphhf.AddPoint(corr1, corr2)
                    #x_minf = min(x_minf, corr1)
                    #x_maxf = max(x_maxf, corr1)
                if (subdet==2) :
                    graphhe.AddPoint(corr1, corr2)
                    #x_mine = min(x_mine, corr1)
                    #x_maxe = max(x_maxe, corr1)
                if (subdet==1) :
                    graphhb.AddPoint(corr1, corr2)
                    #x_minb = min(x_minb, corr1)
                    #x_maxb = max(x_maxb, corr1)
        #otherwise the channel does not exist in file 2. This is handled later in the pull method.
    outroot.WriteObject(graphhb,graphhb.GetName())
    outroot.WriteObject(graphhe,graphhe.GetName())
    outroot.WriteObject(graphhf,graphhf.GetName())

def pull_plots(df1,df2,outpulls,outres,outliers,depth): #pulls and residuals actually
    if(not (df1['Depth']==depth).any() or not (df2['Depth']==depth).any()):return 0 #if the depth doesn't exist at all, don't bother
    #plot corrections:
    #plot pulls
    pulls_th2d=ROOT.TH2D('pull_th2d_d'+str(depth),'Pulls for Depth '+str(depth),83,-41,41,73,0,72)
    resid_th2d=ROOT.TH2D('resid_th2d_d'+str(depth),'Residuals for Depth '+str(depth),83,-41,41,73,0,72)

    residuals_hb=[]
    residuals_he=[]
    residuals_hf=[]
    pulls_hb=[]
    pulls_he=[]
    pulls_hf=[]
    for eta in range(0,pulls_th2d.GetNbinsX()+1):
        for phi in range(0,pulls_th2d.GetNbinsY()+1): 
            pulls_th2d.SetBinContent(eta,phi,-sys.float_info.max)
            resid_th2d.SetBinContent(eta,phi,-sys.float_info.max)

    for row1 in df1.itertuples():
        (subdet,ieta,iphi)=(row1.Subdetector,row1.iEta,row1.iPhi)
        #corr1=df1.query('Subdetector==@subdet and iEta==@ieta and iPhi==@iphi and Depth==@depth')['Correction'].iloc[0]
        corre1=row1.Error
        #corre1=df1.query('Subdetector==@subdet and iEta==@ieta and iPhi==@iphi and Depth==@depth')['Error'].iloc[0]
        #corr2=df2.query('Subdetector==@subdet and iEta==@ieta and iPhi==@iphi and Depth==@depth')['Correction'].iloc[0]
        if(((df2['Subdetector']==subdet) & (df2['iEta']==ieta) & (df2['iPhi']==iphi) & (df2['Depth']==depth)).any()):
            corre2=df2.query('Subdetector==@subdet and iEta==@ieta and iPhi==@iphi and Depth==@depth')['Error'].iloc[0]
            #channel1=(row1.Subdetector,row1.iEta,row1.iPhi,row1.Depth)
            if((math.pow(corre1,2.0)+math.pow(corre2,2.0))==0.0 or row1.Depth!=depth): continue
            corr1=row1.Correction
            corr2=df2.query('Subdetector==@subdet and iEta==@ieta and iPhi==@iphi and Depth==@depth')['Correction'].iloc[0]
            if(corr2!=0):
                residual=get_resid(df1,df2,subdet,ieta,iphi,depth)
                match subdet:
                    case 1: residuals_hb.append(residual)
                    case 2: residuals_he.append(residual)
                    case 4: residuals_hf.append(residual)
                    case _: print("joever")
                if(subdet==4 and (ieta==29 or ieta==-29)):resid_th2d.SetBinContent(ieta+42,iphi+1,residual)
                else: resid_th2d.SetBinContent(ieta+42,iphi,residual)
                pull=get_pull(df1,df2,subdet,ieta,iphi,depth) #Sean indented this line to->

                outpulls.write(str(subdet)+" "+str(ieta)+" "+str(iphi)+' '+str(depth)+' '+str(pull)+'\n')#
                #print(pull)
                outres.write(str(subdet)+" "+str(ieta)+" "+str(iphi)+' '+str(depth)+' '+str(residual)+'\n')#
                match subdet:#
                    case 1: pulls_hb.append(pull)#
                    case 2: pulls_he.append(pull)#
                    case 4: pulls_hf.append(pull)#
                    case _: print('joever')#
                if(subdet==4 and (ieta==29 or ieta==-29)): pulls_th2d.SetBinContent(ieta+42,iphi+1,pull)#
                else:pulls_th2d.SetBinContent(ieta+42,iphi,pull)
                if(abs(pull)>4 and abs(residual)>0.12): outliers.write(str(subdet)+" "+str(ieta)+" "+str(iphi)+" "+str(depth)+" "+str(corr1)+" "+str(corre1)+" "+str(corr2)+" "+str(corre2)+" "+str(pull)+" "+str(residual)+'\n') #when we have a large pull AND large residual, write to outlier file.
            else:#
                print(f'Channel from {filename1} not in {filename2}: {subdetectors[subdet]} ieta {ieta} iphi {iphi} depth {depth}\n') #this line, end of intent
            #We dont need the error unless a channel in 1 file is not present in the other, before indent this just output whenever any line wasnt present in file2, I think. And it was calculating pulls where corr2 was 0...

    pulls_1d=[pulls_hb,pulls_he,pulls_hf]
    for pull_list in pulls_1d:
        if pull_list==[]:continue
        pull_subdet=[]
        if pull_list==pulls_hb:
            pull_subdet='HB'
            pulls = ROOT.TH1D("pulls_1d_"+pull_subdet+"_d"+str(depth), "Pulls for "+pull_subdet + ' Depth '+str(depth), 50, -5, 5)
        elif pull_list==pulls_he:
            pull_subdet='HE'
            pulls = ROOT.TH1D("pulls_1d_"+pull_subdet+"_d"+str(depth), "Pulls for "+pull_subdet + ' Depth '+str(depth), 50, -5, 5)
        else: #elif pull_list==pulls_hf:
            pull_subdet='HF'
            pulls = ROOT.TH1D("pulls_1d_"+pull_subdet+"_d"+str(depth), "Pulls for "+pull_subdet + ' Depth '+str(depth), 50, -10, 10)
        for value in pull_list:
            pulls.Fill(value)
        outroot.WriteObject(pulls,pulls.GetName())

    pulls_th2d.SetMinimum(min(chain(pulls_hb,pulls_he,pulls_hf)))
    outroot.WriteObject(pulls_th2d,pulls_th2d.GetName())

    residuals_1d=[residuals_hb,residuals_he,residuals_hf]

    for res_list in residuals_1d:
        if res_list==[]:continue
        res_subdet='HB'
        if res_list==residuals_he: res_subdet='HE'
        elif res_list==residuals_hf: res_subdet='HF'
        residuals_log = ROOT.TH1D("res_1d_"+res_subdet+"_d"+str(depth)+"_log", "Log(residual) for "+res_subdet + ' Depth '+ str(depth), 50, min(res_list), max(res_list))
        residuals = ROOT.TH1D("res_1d_"+res_subdet+"_d"+str(depth), "Residuals for "+res_subdet + ' Depth '+ str(depth), 50, -0.15, 0.15)
        for value in res_list:
            residuals.Fill(value)
            residuals_log.Fill(value)
        outroot.WriteObject(residuals_log,residuals_log.GetName())
        outroot.WriteObject(residuals,residuals.GetName())

    resid_th2d.SetMinimum(min(chain(residuals_hb,residuals_he,residuals_hf)))
    outroot.WriteObject(resid_th2d,resid_th2d.GetName())

def uncert_plots(df,filename,outfile,depth):
    uncerts=[]
    unc_hist=ROOT.TH2D(f"uncert_f_{filename}_d{depth}",'Relative Uncertainty at Depth '+str(depth),83,-41,41,72,0,72)
    for row in df.itertuples():
        corr=row.Correction
        corre=row.Error
        if corr==0 or row.Depth!=depth:continue
        if(row.Subdetector==4 and (row.iEta==29 or row.iEta==-29)):unc_hist.SetBinContent(row.iEta+42,row.iPhi+1,corre/corr)
        else: unc_hist.SetBinContent(row.iEta+42,row.iPhi,corre/corr)
        uncerts.append(corre/corr)
        outfile.write(str(row.Subdetector)+' '+str(row.iEta)+' '+str(row.iPhi)+' '+str(depth)+' '+str(corr)+' '+str(corre/corr)+'\n')
    outroot.WriteObject(unc_hist,unc_hist.GetName())
    ROOT.gDirectory.Remove(unc_hist)
    
def average_plots(df,inroot_file):
    inroot = ROOT.TFile.Open(inroot_file)
    ROOT.gROOT.SetBatch(True)
    folder_map = {"1": "eHBspec", "2": "eHEspec", "4": "espec"}
    min_x=4
    max_x=150

    for row in df.itertuples():
        (subdet,ieta_formatted,iphi,depth)=(str(row.Subdetector),str(row.iEta),str(row.iPhi),str(row.Depth))
        if subdet not in folder_map:
            raise Exception(f"Folder key {subdet} not recognized.")
        if subdet=='4':max_x=500
        ieta = ieta_formatted if ieta_formatted.startswith('-') else f"+{ieta_formatted}"
        hist_folder=folder_map[subdet]
        ref_hist_path=f'phaseHF/{hist_folder}/E_{ieta}_{iphi}_{depth}_1'

        sum_hist = None
        sum_hist_j = None

        # Loop and sum histograms (over i and j) within x-range 4 to 100
        for i in range(72):
            for j in range(20):
                hist_path=f'phaseHF/{hist_folder}/E_{ieta}_{i+1}_{depth}_{j}'
                hist=inroot.Get(hist_path)

                if not hist:continue

                bin_min=hist.GetXaxis().FindBin(min_x)
                bin_max=hist.GetXaxis().FindBin(max_x)

                if sum_hist is None:
                    sum_hist=hist.Clone()
                    sum_hist.Reset()

                for bin_idx in range(bin_min,bin_max+1):
                    sum_hist.AddBinContent(bin_idx,hist.GetBinContent(bin_idx))
        
        # Loop and sum histograms only over j (fixed i=35) within x-range 4 to 100
        for j in range(20):
            hist_path_j=f'phaseHF/{hist_folder}/E_{ieta}_{iphi}_{depth}_{j}'
            hist_j=inroot.Get(hist_path_j)

            if not hist_j:continue

            bin_min=hist_j.GetXaxis().FindBin(min_x)
            bin_max=hist_j.GetXaxis().FindBin(max_x)

            if sum_hist_j is None:
                sum_hist_j=hist_j.Clone()
                sum_hist_j.Reset()
            
            for bin_idx in range(bin_min,bin_max+1):
                sum_hist_j.AddBinContent(bin_idx,hist_j.GetBinContent(bin_idx))

        hist_ref_original=inroot.Get(ref_hist_path)

        bin_min_ref=hist_ref_original.GetXaxis().FindBin(min_x)
        bin_max_ref=hist_ref_original.GetXaxis().FindBin(max_x)
        hist_name_match = re.search(r'(E_[^/]+)$', ref_hist_path)
        hist_name = hist_name_match.group(1) if hist_name_match else ref_hist_path

        hist_ref=hist_ref_original.Clone()
        hist_ref.Reset()

        for bin_idx in range(bin_min_ref, bin_max_ref + 1):
            hist_ref.SetBinContent(bin_idx, hist_ref_original.GetBinContent(bin_idx))

        integral_ref_sum = sum_hist.Integral(bin_min_ref, bin_max_ref)
        integral_ref_sum_j = sum_hist_j.Integral(bin_min_ref, bin_max_ref)

        sum_hist.Scale(1.0 / integral_ref_sum)
        sum_hist_j.Scale(1.0 / integral_ref_sum_j)
        sum_hist.SetLineColor(ROOT.kRed)
        sum_hist.SetLineWidth(2)
        sum_hist.SetTitle(f"Averaged Histogram vs. {ref_hist_path}")
        sum_hist_j.SetLineColor(ROOT.kGreen + 2)
        sum_hist_j.SetLineWidth(2)
        # Reference histogram
        integral_ref = hist_ref.Integral(bin_min_ref, bin_max_ref)
        if integral_ref != 0:
            hist_ref.Scale(1.0 / integral_ref)
        hist_ref.SetLineColor(ROOT.kBlue)
        hist_ref.SetLineWidth(2)
        if hist_folder=="espec" and abs(int(ieta))>=40 :
            #print("rebinning")
            hist_ref.Rebin(30)
            sum_hist_j.Rebin(30)
            sum_hist.Rebin(30)
        elif hist_folder=="espec":
            #print("rebinning")
            hist_ref.Rebin(5)
            sum_hist.Rebin(5)
            sum_hist_j.Rebin(5)

        if hist_folder == "eHBspec":
            sum_hist.GetXaxis().SetRangeUser(min_x, max_x)
            #print(sum_hist.GetMean())
            hist_ref.GetXaxis().SetRangeUser(min_x, max_x)
            sum_hist_j.GetXaxis().SetRangeUser(min_x, max_x)
        if hist_folder == "eHEspec":
            sum_hist.GetXaxis().SetRangeUser(min_x, max_x)
            #print(sum_hist.GetMean())
            hist_ref.GetXaxis().SetRangeUser(min_x, max_x)
            sum_hist_j.GetXaxis().SetRangeUser(min_x, max_x)
        if hist_folder == "espec":
            sum_hist.GetXaxis().SetRangeUser(min_x, max_x)
            #print(sum_hist.GetMean())
            hist_ref.GetXaxis().SetRangeUser(min_x, max_x)
            sum_hist_j.GetXaxis().SetRangeUser(min_x, max_x)
        
        sum_hist.SetTitle( hist_folder+ " " +hist_name+ " (Mod(20)=1) versus averaged Ring+N")
        save_name = ref_hist_path.replace('/', '_').replace(':', '_')
        sum_hist.SetName(f"averaged_vs_{hist_folder}_{get_filename(inroot_file)}_{save_name}")
        hist_ref.SetName(f"ref_averaged_vs_{hist_folder}_{get_filename(inroot_file)}_{save_name}")
        sum_hist_j.SetName(f"j_averaged_vs_{hist_folder}_{get_filename(inroot_file)}_{save_name}")
        outroot.WriteObject(sum_hist,sum_hist.GetName())
        outroot.WriteObject(hist_ref,hist_ref.GetName())
        outroot.WriteObject(sum_hist_j,sum_hist_j.GetName())

def draw_hist(hist): #for correlation plot it takes TGraph, not hist, BTW.
    histname=hist.GetName()

    if 'd1' in histname: depth=1
    elif 'd2' in histname: depth=2
    elif 'd3' in histname: depth=3
    elif 'd4' in histname: depth=4
    elif 'd5' in histname: depth=5
    elif 'd6' in histname: depth=6
    else: depth=7

    if '1d' in histname:
        ROOT.gStyle.SetOptStat(1111)
        ROOT.gStyle.SetOptFit(1111)
        gaussian = ROOT.TF1("gaussian", "gaus", hist.GetMinimum(), hist.GetMaximum())
        if 'res' in histname and 'log' not in histname: gaussian=ROOT.TF1('gaussian','gaus',-0.15,0.15)
        hist.Fit(gaussian)
        canvas = ROOT.TCanvas("canvas", "Histogram with Gaussian Fit", 1800, 1200) #maybe move canvas outside of conditionals?
        canvas.SetRightMargin(0.15)
        if 'log' in histname: canvas.SetLogy()
        hist.Draw()
        gaussian.Draw("same")
        if 'res' in histname: hist.GetXaxis().SetTitle('Residuals')
        elif 'pull' in histname: hist.GetXaxis().SetTitle('Pulls')
        hist.GetYaxis().SetTitle('Entries')
        if 'HB' in histname: subdet='HB'
        elif 'HE' in histname: subdet='HE'
        else: subdet='HF'
        if 'pull' in histname: canvas.SaveAs(outfolder+'freq_pulls_f1_'+filename1+'_f2_'+filename2+'_'+subdet+'_d'+str(depth)+".png")
        elif 'res' in histname:
            if 'log' in histname: canvas.SaveAs(outfolder+'freq_resid_f1_'+filename1+'_f2_'+filename2+'_'+subdet+'_d'+str(depth)+"_log.png")
            else: canvas.SaveAs(outfolder+'freq_resid_f1_'+filename1+'_f2_'+filename2+'_'+subdet+'_d'+str(depth)+".png")
        canvas.Close()
    elif 'th2d' in histname:
        ROOT.gStyle.SetOptStat(0)
        canvas = ROOT.TCanvas("canvas", "Histogram with Gaussian Fit", 1800, 1200)
        canvas.SetRightMargin(0.15)
        hist.Draw('COLZ')
        if 'pull' in histname:
            canvas.SaveAs(outfolder+'th2d_pulls_f1_'+filename1+'_f2_'+filename2+'_d'+str(depth)+'.png')
            canvas.Close()
            for pthreshold in pull_thresholds:
                hist.SetTitle('\\text{Pulls for Depth }'+str(depth)+'\\text{ within }\\pm'+str(pthreshold))
                hist.SetMaximum(pthreshold)
                hist.SetMinimum(-pthreshold)
                zoomed_pulls_canvas = ROOT.TCanvas("canvas3", "Canvas Title",1800, 1200)
                zoomed_pulls_canvas.SetRightMargin(0.15)
                hist.Draw('COLZ')
                zoomed_pulls_canvas.SaveAs(outfolder+'th2d_pulls_'+strip_float(pthreshold)+'_f1_'+filename1+'_f2_'+filename2+'_d'+str(depth)+'.png')
                zoomed_pulls_canvas.Close()
        elif 'res' in histname:
            canvas.SaveAs(outfolder+'th2d_resid_f1_'+filename1+'_f2_'+filename2+'_d'+str(depth)+'.png')
            canvas.Close()
            for rthreshold in resid_thresholds:
                hist.SetTitle('\\text{Residuals for Depth }'+str(depth)+'\\text{ within }\\pm'+str(rthreshold))
                hist.SetMaximum(rthreshold)
                hist.SetMinimum(-rthreshold)
                zoomed_resid_canvas = ROOT.TCanvas("canvas5", "Canvas Title", 1800, 1200)
                zoomed_resid_canvas.SetRightMargin(0.15)
                hist.Draw('COLZ')
                zoomed_resid_canvas.SaveAs(outfolder+'th2d_resid_'+strip_float(rthreshold)+'_f1_'+filename1+'_f2_'+filename2+'_d'+str(depth)+'.png')
                zoomed_resid_canvas.Close()
        else: #just to be safe
            canvas.Close()
    elif 'uncert' in histname:
        ROOT.gStyle.SetOptStat(0)
        canvas = ROOT.TCanvas("canvas", "Canvas Title", 1800, 1200)
        canvas.SetRightMargin(0.15)
        hist.Draw("COLZ")
        if filename1 in histname: imagename=outfolder+"uncert_f_"+filename1+"_d"+str(depth)+".png"
        else: imagename=outfolder+"uncert_f_"+filename2+"_d"+str(depth)+".png"
        canvas.SaveAs(imagename)
        canvas.Close()
    elif 'correlation' in histname: #this uses TGraph 
        ROOT.gStyle.SetOptStat(0)
        canvas = ROOT.TCanvas("canvas", "Scatter Plot", 1800, 1200)
        xy_line = ROOT.TF1("identity", "x", hist.GetMinimum(), hist.GetMaximum())  # y = x
        xy_line.SetLineColor(ROOT.kBlack)  # Black color
        xy_line.SetLineWidth(2)  # Line thickness
        pearson_r = hist.GetCorrelationFactor()
        fit_function = ROOT.TF1("fit", "pol1", 0, 10)
        hist.Fit(fit_function, "Q")
        hist.Draw("AP")
        xy_line.Draw("SAME")
        legend = ROOT.TLegend(0.15, 0.75, 0.4, 0.85)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.03)
        dummy_graph = ROOT.TGraph()  # Invisible graph
        legend.AddEntry(dummy_graph, f"#rho = {pearson_r:.4f}", "")
        legend.AddEntry(xy_line, "y = x", "l")  # "l" = line
        legend.AddEntry(fit_function, "Linear Fit", "l")
        legend.Draw()
        imagename=outfolder+histname+'.png'
        canvas.SaveAs(imagename)
        canvas.Close()
    elif 'corrs' in histname:
        ROOT.gStyle.SetOptStat(0)
        canvas = ROOT.TCanvas("canvas", "Canvas Title", 1800, 1200)
        canvas.SetRightMargin(0.15)
        hist.Draw('COLZ')
        if filename1 in histname: imagename=outfolder+"corrs_f_"+filename1+"_d"+str(depth)+".png"
        else: imagename=outfolder+"corrs_f_"+filename2+"_d"+str(depth)+".png"
        canvas.SaveAs(imagename)
        canvas.Close()
    elif 'ref' in histname:
        min_x=4
        max_x=500 if 'espec' in histname else 150
        ROOT.gStyle.SetOptStat(0)
        sum_hist=outroot.Get(histname[4:])
        sum_hist_j=outroot.Get('j_'+histname[4:])
        canvas = ROOT.TCanvas("canvas", "Histogram Comparison", 1800, 1200)
        max_y = max(sum_hist.GetMaximum(), hist.GetMaximum()) * 1.1
        sum_hist.SetMaximum(max_y)
        hist.SetMaximum(max_y)
        sum_hist.Draw('HIST')
        hist.Draw('HIST SAME')
        
        bin_min_ref=hist.GetXaxis().FindBin(min_x)
        bin_max_ref=hist.GetXaxis().FindBin(max_x)
        integral_ref = hist.Integral(bin_min_ref, bin_max_ref)
        if integral_ref != 0:
            hist.Scale(1.0 / integral_ref)

        integral_ref_sum = sum_hist.Integral(bin_min_ref, bin_max_ref)
        integral_ref_sum_j = sum_hist_j.Integral(bin_min_ref, bin_max_ref)
        meanxintegral_hist_ref=hist.GetMean()*integral_ref

        meanxintegral_sum_hist=sum_hist.GetMean()*integral_ref_sum/72

        legend = ROOT.TLegend(0.55, 0.65, 0.85, 0.85)
        legend.AddEntry(hist, f"E (mean x = {meanxintegral_hist_ref:.2f})", "l")
        legend.AddEntry(sum_hist, f"Averaged Ring and N: (mean x = {meanxintegral_sum_hist:.2f})", "l")
        legend.Draw()

        canvas.SaveAs(outfolder+sum_hist.GetName()+'.png')

        max_y = max(sum_hist_j.GetMaximum(), sum_hist.GetMaximum()) * 1.1
        sum_hist_j.SetMaximum(max_y)
        sum_hist.SetMaximum(max_y)
        sum_hist.SetTitle(sum_hist.GetTitle().replace("(Mod(20)=1)","averaged N"))
        sum_hist.Draw('HIST')
        sum_hist_j.Draw('HIST SAME')
        meanxintegral_sum_hist_j=sum_hist_j.GetMean()*integral_ref_sum_j

        legend = ROOT.TLegend(0.55, 0.65, 0.85, 0.85)
        legend.AddEntry(sum_hist_j, f"Averaged N: (mean x = {meanxintegral_sum_hist_j:.2f})", "l")
        legend.AddEntry(sum_hist, f"Averaged Ring and N: (mean x = {meanxintegral_sum_hist:.2f})", "l")
        legend.Draw() 

        canvas.SaveAs(outfolder+sum_hist.GetName()+'_N.png')
        canvas.Close()


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("file1",help="first file")
    parser.add_argument("file2",help="second file")
    parser.add_argument('--outdir',type=str,help='output directory',default='')
    parser.add_argument('--pullthresholds',type=str_to_float_list,help='pull thresholds, comma separated, abs',default='10')
    parser.add_argument('--resthresholds',type=str_to_float_list,help='residual thresholds, comma separated, abs',default='0.15')
    args=parser.parse_args()

    global outfolder
    outfolder=str(args.outdir)
    if outfolder!='':
        if outfolder[-1]!='/': outfolder=outfolder+'/'
        if not(os.path.exists(outfolder) and os.path.isdir(outfolder)): os.makedirs(outfolder)

    global subdetectors
    subdetectors={1:'HB',2:'HE',4:'HF'}
    ROOT.gStyle.SetOptStat(0)
    file1=args.file1
    file2=args.file2

    global filename1,filename2 #get pure filenames to use globally for new output filenames
    filename1=get_filename(file1)
    filename2=get_filename(file2)

    global outroot #root file to store everything
    outroot=ROOT.TFile(outfolder+'rootfile_f1_'+filename1+'_f2_'+filename2+'.root','recreate')

    if get_file_extension(file1)=='root' and get_file_extension(file2)=='txt':
        df=get_df2(file2)
        average_plots(df,file1)

    elif get_file_extension(file2)=='root' and get_file_extension(file1)=='txt':
        df=get_df2(file1)
        average_plots(df,file2)

    elif get_file_extension(file1)=='txt' and get_file_extension(file2)=='txt':
        df1=get_df(file1) #get dfs
        df2=get_df(file2)
        
        global pull_thresholds, resid_thresholds
        pull_thresholds=args.pullthresholds
        resid_thresholds=args.resthresholds

        #text outputs. change as appropriate.
        pulls=open(outfolder+'pulls_f1_'+filename1+"_f2_"+filename2+'.txt','w')
        uncerts1=open(outfolder+'uncerts_'+filename1+'.txt','w')
        uncerts2=open(outfolder+'uncerts_'+filename2+'.txt','w')
        residuals=open(outfolder+'residuals_f1_'+filename1+"_f2_"+filename2+'.txt','w')
        outliers=open(outfolder+'outliers_f1_'+filename1+"_f2_"+filename2+'.txt','w')
        #corrs1=open('corrs_'+get_filename(file1)+'.txt','w')
        #corrs2=open('corrs_'+get_filename(file2)+'.txt','w')

        correlation_plot(df1,df2,False)
        correlation_plot(df1,df2,True)
        for depth in range(1,8):
            pull_plots(df1,df2,pulls,residuals,outliers,depth)
            uncert_plots(df1,filename1,uncerts1,depth)
            uncert_plots(df2,filename2,uncerts2,depth)
            corr_plots(df1,filename1,depth)
            corr_plots(df2,filename2,depth)

        pulls.close()
        uncerts1.close()
        uncerts2.close()
        residuals.close()
        outliers.close()
    
    else:
        raise Exception("Invalid file inputs. Inputs must be either both .txt files or 1 .root file and 1 .txt file.")

    for key in outroot.GetListOfKeys():
        draw_hist(key.ReadObj())

    outroot.Close()
if __name__ == "__main__":
    main()
