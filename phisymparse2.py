import ROOT
import argparse
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from itertools import chain
ROOT.gROOT.SetBatch(True)
import os

def get_df(filename):
    df=pd.read_csv(filename,sep='\s+',header=None)
    df=df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    df.columns = ['Subdetector', 'iEta', 'iPhi', 'Depth', 'Correction', 'Error']
    filt= df[(df['Depth'] == 1) & (df['Subdetector'] == 1)& (df['iEta']<=-1 )& (df['iEta']>=-16 ) &((df['iPhi']!=16)&(df['iEta']!=-16))
                      &((df['iPhi']!=44)&(df['iEta']!=-8)) &((df['iPhi']!=10)&(df['iEta']!=-2))]
    filt = filt.reset_index(drop=True)
    return df

def get_filename(filename):
    return os.path.basename(filename)[:-4] #cut off the directory

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

def corr_plots(file1,depth):
    df=get_df(file1)
    corrs=ROOT.TH2D('corrs','Corrections for Depth '+str(depth),83,-41,41,73,0,72)
    for row in df.itertuples():
        corr=row.Correction
        if (row.Depth!=depth):continue
        if(row.Subdetector==4 and (row.iEta==29 or row.iEta==-29)): corrs.SetBinContent(row.iEta+42,row.iPhi+1,corr)
        else:corrs.SetBinContent(row.iEta+42,row.iPhi,corr)
    canvas = ROOT.TCanvas("canvas", "Canvas Title",1800, 1200)
    canvas.SetRightMargin(0.15)
    corrs.Draw("COLZ")
    imagename=outfolder+get_filename(file1)+"_d"+str(depth)+".png"
    canvas.SaveAs(imagename)

def correlation_plot(file1,file2):
    df1=get_df(file1)
    df2=get_df(file2)
    graphhb = ROOT.TGraph()
    graphhb.SetTitle("Correlation for Barrel")
    graphhb.SetMarkerStyle(20)
    graphhb.SetMarkerSize(0.5)
    graphhb.SetMarkerColor(ROOT.kBlue)
    graphhb.GetXaxis().SetTitle(get_filename(file1))  # Set X-axis title
    graphhb.GetYaxis().SetTitle(get_filename(file2))
    graphhb.GetXaxis().SetLimits(0, 12)
    graphhe = ROOT.TGraph()
    graphhe.SetTitle("Correlation for Endcap")
    graphhe.SetMarkerStyle(20)
    graphhe.SetMarkerSize(0.5)
    graphhe.SetMarkerColor(ROOT.kBlue)
    graphhe.GetXaxis().SetTitle(get_filename(file1))  # Set X-axis title
    graphhe.GetYaxis().SetTitle(get_filename(file2))
    graphhf = ROOT.TGraph()
    graphhf.SetTitle("Correlation for Forward")
    graphhf.SetMarkerStyle(20)
    graphhf.SetMarkerSize(0.5)
    graphhf.SetMarkerColor(ROOT.kBlue)
    graphhf.GetXaxis().SetTitle(get_filename(file1))  # Set X-axis title
    graphhf.GetYaxis().SetTitle(get_filename(file2))
    x_minb, x_maxb = float("inf"), float("-inf")
    x_mine, x_maxe = float("inf"), float("-inf")
    x_minf, x_maxf = float("inf"), float("-inf")
    for row in df1.itertuples():
        (subdet,ieta,iphi,depth)=(row.Subdetector,row.iEta,row.iPhi,row.Depth)
        corr1=row.Correction
        corr2=df2.query('Subdetector==@subdet and iEta==@ieta and iPhi==@iphi and Depth==@depth')['Correction'].iloc[0]
        if(((df2['Subdetector']==subdet) & (df2['iEta']==ieta) & (df2['iPhi']==iphi) & (df2['Depth']==depth)).any()):
            if(corr2!=0):
                if (subdet==4) :
                    graphhf.AddPoint(corr1, corr2)
                    x_minf = min(x_minf, corr1)
                    x_maxf = max(x_maxf, corr1)
                if (subdet==2) :
                    graphhe.AddPoint(corr1, corr2)
                    x_mine = min(x_mine, corr1)
                    x_maxe = max(x_maxe, corr1)
                if (subdet==1) :
                    graphhb.AddPoint(corr1, corr2)
                    x_minb = min(x_minb, corr1)
                    x_maxb = max(x_maxb, corr1)
        #otherwise the channel does not exist in file 2. This is handled later in the pull method. 
    canvas = ROOT.TCanvas("canvas", "Scatter Plot", 1800, 1200)
    xy_lineB = ROOT.TF1("identityB", "x", x_minb, x_maxb)  # y = x
    xy_lineB.SetLineColor(ROOT.kBlack)  # Black color
    xy_lineB.SetLineWidth(2)  # Line thickness
    xy_lineB.SetLineStyle(3)  # Dotted line (ROOT style 3)
    pearson_r_hb = graphhb.GetCorrelationFactor()
    fit_functionb = ROOT.TF1("fit", "pol1", 0, 10)
    graphhb.Fit(fit_functionb, "Q")
    graphhb.Draw("AP")
    xy_lineB.Draw("SAME")
    legend = ROOT.TLegend(0.15, 0.75, 0.4, 0.85)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)
    dummy_graph = ROOT.TGraph()  # Invisible graph
    legend.AddEntry(dummy_graph, f"#rho = {pearson_r_hb:.4f}", "")
    legend.AddEntry(xy_lineB, "y = x", "l")  # "l" = line
    legend.AddEntry(fit_functionb, "Linear Fit", "l")
    legend.Draw()
    imagename=outfolder+"correlation_"+get_filename(file1)+get_filename(file2)+"_hb"+".png"
    canvas.SaveAs(imagename)
    xy_lineE = ROOT.TF1("identityE", "x", x_mine, x_maxe)
    xy_lineE.SetLineColor(ROOT.kBlack)  # Black color
    xy_lineE.SetLineWidth(2)  # Line thickness
    xy_lineE.SetLineStyle(3)  # Dotted line (ROOT style 3)
    pearson_r_he = graphhe.GetCorrelationFactor()
    fit_functione = ROOT.TF1("fit", "pol1", 0, 10)
    graphhe.Fit(fit_functione, "Q")
    graphhe.Draw("AP")
    xy_lineE.Draw("SAME")
    legend = ROOT.TLegend(0.15, 0.75, 0.4, 0.85)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)
    legend.AddEntry(dummy_graph, f"#rho = {pearson_r_he:.4f}", "")
    legend.AddEntry(xy_lineE, "y = x", "l")  # "l" = line
    legend.AddEntry(fit_functione, "Linear Fit", "l")
    legend.Draw()
    imagename=outfolder+"correlation_"+get_filename(file1)+get_filename(file2)+"_he"+".png"
    canvas.SaveAs(imagename)
    xy_lineF = ROOT.TF1("identityF", "x", x_minf, x_maxf)
    xy_lineF.SetLineColor(ROOT.kBlack)  # Black color
    xy_lineF.SetLineWidth(2)  # Line thickness
    xy_lineF.SetLineStyle(3)  # Dotted line (ROOT style 3)
    pearson_r_hf = graphhf.GetCorrelationFactor()
    fit_functionf = ROOT.TF1("fit", "pol1", 0, 10)
    graphhf.Fit(fit_functionf, "Q")
    graphhf.Draw("AP")
    xy_lineF.Draw("SAME")
    legend = ROOT.TLegend(0.15, 0.75, 0.4, 0.85)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)
    legend.AddEntry(dummy_graph, f"#rho = {pearson_r_hf:.4f}", "")
    legend.AddEntry(xy_lineF, "y = x", "l")  # "l" = line
    legend.AddEntry(fit_functionf, "Linear Fit", "l")
    legend.Draw()
    imagename=outfolder+"correlation_"+get_filename(file1)+get_filename(file2)+"_hf"+".png"
    canvas.SaveAs(imagename)
    
def pull_plots(file1,file2,outpulls,outres,outliers,depth):
    df1=get_df(file1)
    df2=get_df(file2)
    #plot corrections:
    #plot pulls
    pulls_th2d=ROOT.TH2D('pull_th2d','Pulls for Depth '+str(depth),83,-41,41,73,0,72)
    resid_th2d=ROOT.TH2D('resid_th2d','Residuals for Depth '+str(depth),83,-41,41,73,0,72)

    residuals_hb=[]
    residuals_he=[]
    residuals_hf=[]
    pulls_hb=[]
    pulls_he=[]
    pulls_hf=[]
    for eta in range(0,pulls_th2d.GetNbinsX()+1):
        for phi in range(0,pulls_th2d.GetNbinsY()+1):
            pulls_th2d.SetBinContent(eta,phi,-1000)
            resid_th2d.SetBinContent(eta,phi,-1000)

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
                if(abs(pull)>4 and abs(residual)>0.12): outliers.write(str(subdet)+" "+str(ieta)+" "+str(iphi)+' '+str(depth)+' '+str(pull)+' '+str(residual)+'\n') #when we have a large pull AND large residual, write to outlier file.
            else:#
                print(f'Channel from {file1} not in {file2}: {subdetectors[subdet]} ieta {ieta} iphi {iphi} depth {depth}\n') #this line, end of intent
            #We dont need the error unless a channel in 1 file is not present in the other, before indent this just output whenever any line wasnt present in file2, I think. And it was calculating pulls where corr2 was 0...

    pulls_1d=[pulls_hb,pulls_he,pulls_hf]
    for pull_list in pulls_1d:
        if pull_list==[]:continue
        pull_subdet=[]
        ROOT.gStyle.SetOptStat(1111)  # Display stats (1 for entries, 1 for mean, 1 for RMS, etc.)
        ROOT.gStyle.SetOptFit(1111)
        if pull_list==pulls_hb:
            pull_subdet='HB'
            pulls = ROOT.TH1D("pulls_1d_"+pull_subdet, "Pulls for "+pull_subdet + ' Depth '+str(depth), 50, -5, 5)
        if pull_list==pulls_he:
            pull_subdet='HE'
            pulls = ROOT.TH1D("pulls_1d_"+pull_subdet, "Pulls for "+pull_subdet + ' Depth '+str(depth), 50, -5, 5)
        elif pull_list==pulls_hf:
            pull_subdet='HF'
            pulls = ROOT.TH1D("pulls_1d_"+pull_subdet, "Pulls for "+pull_subdet + ' Depth '+str(depth), 50, -10, 10)
        for value in pull_list:
            pulls.Fill(value)
        gaussian = ROOT.TF1("gaussian", "gaus", min(pull_list), max(pull_list))
        pulls.Fit(gaussian)
        canvas = ROOT.TCanvas("canvas", "Histogram with Gaussian Fit", 1800, 1200)
        canvas.SetRightMargin(0.15)
        pulls.Draw()
        gaussian.Draw("same")
        pulls.GetXaxis().SetTitle("Pulls")
        pulls.GetYaxis().SetTitle("Entries")
        canvas.SaveAs(outfolder+'freq_pulls_f1_'+get_filename(file1)+'_f2_'+get_filename(file2)+'_'+pull_subdet+'_d'+str(depth)+".png")
        canvas.Close()

    ROOT.gStyle.SetOptStat(0)
    pulls_canvas = ROOT.TCanvas("canvas", "Canvas Title1", 1800, 1200)
    pulls_canvas.SetRightMargin(0.15)
    pulls_th2d.SetMinimum(min(chain(pulls_hb,pulls_he,pulls_hf)))
    pulls_th2d.Draw('COLZ')
    pulls_canvas.SaveAs(outfolder+'th2d_pulls_f1_'+get_filename(file1)+'_f2_'+get_filename(file2)+'_d'+str(depth)+'.png')
    pulls_canvas.Close()

    residuals_1d=[residuals_hb,residuals_he,residuals_hf]

    for res_list in residuals_1d:
        if res_list==[]:continue
        res_subdet='HB'
        if res_list==residuals_he: res_subdet='HE'
        elif res_list==residuals_hf: res_subdet='HF'
        residuals_log = ROOT.TH1D("res_1d_"+res_subdet, "Residuals for "+res_subdet + ' Depth '+ str(depth), 50, min(res_list), max(res_list))
        residuals = ROOT.TH1D("res_1d_"+res_subdet, "Residuals for "+res_subdet + ' Depth '+ str(depth), 50, -0.05, 0.05)
        for value in res_list:
            residuals.Fill(value)
            residuals_log.Fill(value)
        gaussian = ROOT.TF1("gaussian2", "gaus", min(res_list), max(res_list))
        residuals_log.Fit(gaussian)
        ROOT.gStyle.SetOptStat(1111)  # Display stats (1 for entries, 1 for mean, 1 for RMS, etc.)
        ROOT.gStyle.SetOptFit(1111)
        canvas = ROOT.TCanvas("canvas", "Histogram with Gaussian Fit", 1800, 1200)
        canvas.SetRightMargin(0.15)
        canvas.SetLogy()
        residuals_log.Draw()
        gaussian.Draw("same")
        residuals_log.GetXaxis().SetTitle("Residuals")
        residuals_log.GetYaxis().SetTitle("Entries")
        canvas.SaveAs(outfolder+'freq_resid_f1_'+get_filename(file1)+'_f2_'+get_filename(file2)+'_'+res_subdet+'_d'+str(depth)+"log"+".png")
        canvas.Close()

        gaussian = ROOT.TF1("gaussian2", "gaus", -0.05, 0.05)
        residuals.Fit(gaussian)
        ROOT.gStyle.SetOptStat(1111)  # Display stats (1 for entries, 1 for mean, 1 for RMS, etc.)
        ROOT.gStyle.SetOptFit(1111)
        canvas = ROOT.TCanvas("canvas", "Histogram with Gaussian Fit", 1800, 1200)
        canvas.SetRightMargin(0.15)
        residuals.Draw()
        gaussian.Draw("same")
        residuals.GetXaxis().SetTitle("Residuals")
        residuals.GetYaxis().SetTitle("Entries")
        canvas.SaveAs(outfolder+'freq_resid_f1_'+get_filename(file1)+'_f2_'+get_filename(file2)+'_'+res_subdet+'_d'+str(depth)+".png")
        canvas.Close()
    #for res_plot in residuals_1d:
        #plt.figure()
        #plt.hist(res_plot, bins=50, color='skyblue', edgecolor='black',density=True, alpha=0.6)
        #plt.xlabel('Residual')
        #plt.ylabel('Frequency')
        #mu, std = norm.fit(res_plot)
        #xmin, xmax = plt.xlim()
        #x = np.linspace(xmin, xmax, 100)
        #p = norm.pdf(x, mu, std)
        #plt.plot(x, p, 'r--', linewidth=2, label=f'Gaussian fit: $\mu$={mu:.2f}, $\sigma$={std:.2f}')

        #plt.savefig('freq_resid_f1_'+file1[:-4]+'_f2'+file2[:-4]+'_'+residuals_1d[res_plot]+'_d'+str(depth)+".png")
    ROOT.gStyle.SetOptStat(0)
    resid_canvas = ROOT.TCanvas("canvas2", "Canvas Title2", 1800, 1200)
    resid_canvas.SetRightMargin(0.15)
    resid_th2d.SetMinimum(min(chain(residuals_hb,residuals_he,residuals_hf)))
    resid_th2d.Draw('COLZ')
    resid_canvas.SaveAs(outfolder+'th2d_resid_f1_'+get_filename(file1)+'_f2_'+get_filename(file2)+'_d'+str(depth)+'.png')
    resid_canvas.Close()

    pulls_th2d.SetTitle('Pulls for Depth '+str(depth)+' within +-5')
    pulls_th2d.SetMaximum(5.0)
    pulls_th2d.SetMinimum(-5.0)
    zoomed_pulls_canvas = ROOT.TCanvas("canvas3", "Canvas Title",1800, 1200)
    zoomed_pulls_canvas.SetRightMargin(0.15)
    pulls_th2d.Draw('COLZ')
    zoomed_pulls_canvas.SaveAs(outfolder+'th2d_pulls_5_f1_'+get_filename(file1)+'_f2_'+get_filename(file2)+'_d'+str(depth)+'.png')
    zoomed_pulls_canvas.Close()

    pulls_th2d.SetTitle('Pulls for Depth '+str(depth)+' within +-10')
    pulls_th2d.SetMaximum(10.0)
    pulls_th2d.SetMinimum(-10.0)
    zoomed_pulls_canvas1 = ROOT.TCanvas("canvas4", "Canvas Title",1800, 1200)
    zoomed_pulls_canvas1.SetRightMargin(0.15)
    pulls_th2d.Draw('COLZ')
    zoomed_pulls_canvas1.SaveAs(outfolder+'th2d_pulls_10_f1_'+get_filename(file1)+'_f2_'+get_filename(file2)+'_d'+str(depth)+'.png')
    zoomed_pulls_canvas1.Close()

    resid_th2d.SetTitle('Residuals for Depth '+str(depth)+' within +-0.05')
    resid_th2d.SetMaximum(0.05)
    resid_th2d.SetMinimum(-0.05)
    zoomed_resid_canvas = ROOT.TCanvas("canvas5", "Canvas Title", 1800, 1200)
    zoomed_resid_canvas.SetRightMargin(0.15)
    resid_th2d.Draw('COLZ')
    zoomed_resid_canvas.SaveAs(outfolder+'th2d_resid_005_f1_'+get_filename(file1)+'_f2_'+get_filename(file2)+'_d'+str(depth)+'.png')
    zoomed_resid_canvas.Close()


def uncert_plots(file1,outfile,depth):
    df=get_df(file1)
    uncerts=[]
    unc_hist=ROOT.TH2D('dep1','Relative Uncertainty at Depth '+str(depth),83,-41,41,72,0,72)
    for row in df.itertuples():
        corr=row.Correction
        corre=row.Error
        if corr==0 or row.Depth!=depth:continue
        if(row.Subdetector==4 and (row.iEta==29 or row.iEta==-29)):unc_hist.SetBinContent(row.iEta+42,row.iPhi+1,corre/corr)
        else: unc_hist.SetBinContent(row.iEta+42,row.iPhi,corre/corr)
        uncerts.append(corre/corr)
        outfile.write(str(row.Subdetector)+' '+str(row.iEta)+' '+str(row.iPhi)+' '+str(depth)+' '+str(corr)+' '+str(corre/corr)+'\n')
    canvas = ROOT.TCanvas("canvas", "Canvas Title", 1800, 1200)
    canvas.SetRightMargin(0.15)
    unc_hist.Draw("COLZ")
    imagename=outfolder+"uncert_"+get_filename(file1)+"_d"+str(depth)+".png"
    canvas.SaveAs(imagename)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("file1",help="first file")
    parser.add_argument("file2",help="second file")
    parser.add_argument('--outdir',type=str,help='output directory',default='')
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
    #iter_file = "Run2023_Iter.txt"
    #iter_file='corrs_all_egamma.txt'
    #mom_file="mom.txt"
    #it1='corrs_all_EGAMMA0_2024I.txt'

    #text outputs. change as appropriate.
    pulls=open(outfolder+'pulls_f1_'+get_filename(file1)+"_f2_"+get_filename(file2)+'.txt','w')
    uncerts1=open(outfolder+'uncerts_'+get_filename(file1)+'.txt','w')
    uncerts2=open(outfolder+'uncerts_'+get_filename(file2)+'.txt','w')
    residuals=open(outfolder+'residuals_f1_'+get_filename(file1)+"_f2_"+get_filename(file2)+'.txt','w')
    outliers=open(outfolder+'outliers_f1_'+get_filename(file1)+"_f2_"+get_filename(file2)+'.txt','w')
    #corrs1=open('corrs_'+get_filename(file1)+'.txt','w')
    #corrs2=open('corrs_'+get_filename(file2)+'.txt','w')

    #uncert_plots(iter_file,1)
    #pull_plots(file1,file2,pulls,residuals,1)
    correlation_plot(file1,file2)
    for depth in range(1,8):
        pull_plots(file1,file2,pulls,residuals,outliers,depth)
        uncert_plots(file1,uncerts1,depth)
        uncert_plots(file2,uncerts2,depth)
        corr_plots(file1,depth)
        corr_plots(file2,depth)

    pulls.close()
    uncerts1.close()
    uncerts2.close()
    residuals.close()

if __name__ == "__main__":
    main()
