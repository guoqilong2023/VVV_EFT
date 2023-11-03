import ROOT
import Histograms
reload(Histograms)
from Histograms import *
import GetParameter_Range

import math
from prettytable import PrettyTable
import ROOT
import json

def F_Get_CDRatio(hC,hD):
    c = ROOT.TCanvas("canvas","",700,1000)
    g_R_est = ROOT.TGraphAsymmErrors()
    g_R_est.Divide(hC, hD, "pois w")
    g_R_est.Draw('pela')
    c.Modified()
    c.Update()
    c.SaveAs("plots/abcd/ratio.pdf")
    return g_R_est

def F_PredictA(hB,hC,hD,yields=False):
    if yields :
        hB.Rebin(hB.GetNbinsX())
        hC.Rebin(hC.GetNbinsX())
        hD.Rebin(hD.GetNbinsX())
    gR = F_Get_CDRatio(hC,hD)
    Nbins = hB.GetNbinsX()
    h_Est = hB.Clone("h_Est")
    for b in xrange(1, Nbins+1):
        n_b = hB.GetBinContent(b)
        ## graph points count from 0, thanks root, also -1 for empty 1st bin
        bg = b-1
        r = gR.GetY()[bg]
        r_eup = gR.GetErrorYhigh(bg)
        r_edw = gR.GetErrorYlow(bg)
        print " R ", r, r_eup, r_edw
        est_e = r*n_b
        if n_b > 0:
            unc_Eup = est_e*math.sqrt(1./n_b + (r_eup/r)**2)
            unc_Edw = est_e*math.sqrt(1./n_b + (r_edw/r)**2)
        else :
            unc_Eup = 0
            unc_Edw = 0
        h_Est.SetBinContent(b,est_e)
        h_Est.SetBinError(b,unc_Eup)
    return h_Est

def F_DMR(Histograms,PostFix):
    hDMR = Histograms["h_data"].Clone("h_DMR"+PostFix)
    hDMR.Add(Histograms["h_TTbar"],-1)
    hDMR.Add(Histograms["h_STop"],-1)
    hDMR.Add(Histograms["h_VV"],-1)
    hDMR.Add(Histograms["h_QCD"],-1)
    return hDMR

def F_DMR_2(Histograms,HistogramsTop,Region):
    hTop = F_Top(Histograms,HistogramsTop,Region)
    return Histograms["h_data"+Region]-hTop-Histograms["h_VV"+Region]-Histograms["h_QCD"+Region]

def F_DMR_3(Histograms,HistogramsTop,Region):
    hTop = F_Top_2(Histograms,HistogramsTop,Region)
    return Histograms["h_data"+Region]-hTop-Histograms["h_VV"+Region]-Histograms["h_QCD"+Region]

def F_Top(HistogramSR,HistogramsTop,Region):
    Hist = (HistogramsTop["h_TTbar"+Region+"Top"]+HistogramsTop["h_STop"+Region+"Top"]).Clone("h_Top"+Region)
    Hist.Scale((HistogramSR["h_TTbar"+Region]+HistogramSR["h_STop"+Region]).Integral()/Hist.Integral())
    return Hist

def F_Top_2(HistogramSR,HistogramsTop,Region):
    Hist = HistogramsTop["h_data"+Region+"Top"]-HistogramsTop["h_WJets"+Region+"Top"]-HistogramsTop["h_VV"+Region+"Top"]-HistogramsTop["h_QCD"+Region+"Top"]
    Hist.Scale((HistogramSR["h_TTbar"+Region]+HistogramSR["h_STop"+Region]).Integral()/Hist.Integral())
    return Hist

def F_EstTable(h_Est,h_Data):
    hR = h_Est.Clone("hR")
    hR.Divide(h_Data)
    tb = PrettyTable()
    tb.field_names = ["iBin","Estimated","[Data-rest]","Est/[Data-rest]"]
    Nbins = h_Data.GetNbinsX()
    for b in xrange(1, Nbins+1):
        Est = h_Est.GetBinContent(b)
        Data = h_Data.GetBinContent(b)
        ratio = hR.GetBinContent(b)
        e_Est = h_Est.GetBinError(b)
        e_Data = h_Data.GetBinError(b)
        e_ratio = hR.GetBinError(b)
        
        Est = "%.3g +- %.3g"%(Est,e_Est)
        Data = "%.3g +- %.3g"%(Data,e_Data)
        ratio = "%.3g +- %.3g"%(ratio,e_ratio)
        raws = ["Bin%s"%(b),Est,Data,ratio]
        tb.add_row(raws)
    print tb
    
from prettytable import PrettyTable
def F_ProcessBinYields(Histograms,nround = 1):
    tb = PrettyTable()
    process = ['Bin']
    for Hist in Histograms:
        process.append(Hist['Name'])
    tb.field_names = process
    Nbins = Histograms[0]['Hist'].GetNbinsX()
    for b in xrange(1, Nbins+1):
        L = Histograms[0]['Hist'].GetBinLowEdge(b)
        H = Histograms[0]['Hist'].GetBinLowEdge(b+1)
        bin_ = "%s,%s"%(str(int(L)),str(int(H)))
        raws = [bin_]
        for Hist in Histograms:
            v = round(Hist['Hist'].GetBinContent(b),nround)
            e = round(Hist['Hist'].GetBinError(b),nround)
            raws.append("%s +- %s"%(str(v),str(e)))
        tb.add_row(raws)
    print tb
    return tb

from array import array

def F_RebinHistograms(Histograms,Ngroup=None,varbins=None,postfix=''):
    Histograms_new = {}
    if Ngroup:
        for n in Histograms :
            Histogram = Histograms[n].Clone(Histograms[n].GetName()+postfix)
            Histogram.Rebin(Ngroup)
            Histograms_new[n] = Histogram
    if varbins:
        histo_bin = array('d',varbins)
        for n in Histograms :
            Name = Histograms[n].GetName()+postfix
            Histogram = Histograms[n].Rebin(len(histo_bin)-1,Name,histo_bin)
            Histograms_new[n] = Histogram
    return Histograms_new

def F_MergeHistograms(HistogramsList,postfix=""):
    process = [p for p in HistogramsList[0]]
    Histograms_new = {}
    for p in process :
        for index,Histograms in enumerate(HistogramsList):
            if index == 0 :
                Name = Histograms[p].GetName()+postfix
                Histograms_new[p] = Histograms[p].Clone(Name)
            else :
                Histograms_new[p].Add(Histograms[p])
    return Histograms_new

def F_TF_Mutiply_TH(TF,histogram):
    Name = histogram.GetName()+"_TF1"
    histout = histogram.Clone(Name)
    for i in range(1,histout.GetNbinsX()+1):
        v = histogram.GetBinContent(i)*TF.Eval(histout.GetBinCenter(i))
        histout.SetBinContent(i,v)
    return histout

def F_GroupHistogramDict(Histograms,postfix = ""):
    data = Histograms['h_data'].Clone(Histograms['h_data'].GetName()+postfix)
    Top = Histograms['h_TTbar'].Clone(Histograms['h_TTbar'].GetName()+postfix)
    Top.Add(Histograms['h_STop'])
    WJets = Histograms['h_WJets'].Clone(Histograms['h_WJets'].GetName()+postfix)
    Others = Histograms['h_VV'].Clone(Histograms['h_VV'].GetName()+postfix)
    Others.Add(Histograms['h_QCD'])
    new_Histograms = {
        "data" : data,
        "Top" : Top,
        "WJets" : WJets,
        "Others" : Others,
    }
    return new_Histograms

def F_GroupHistogram(Histograms):
    Hists = [
        {
            "Name" : "data",
            "Hist" : Histograms['h_data'],
        },
        {
            "Name" : "Top",
            "Hist" : Histograms['h_TTbar']+Histograms['h_STop'],   
        },
        {
            "Name" : "W+jets",
            "Hist" : Histograms['h_WJets'],   
        },
        {
            "Name" : "Others",
            "Hist" : Histograms['h_VV']+Histograms['h_QCD'],   
        }
    ]
    return Hists

def F_GroupHistogram_Nodata(Histograms):
    Hists = [
        {
            "Name" : "W+jets",
            "Hist" : Histograms['h_WJets'],   
        },
        {
            "Name" : "Top",
            "Hist" : Histograms['h_TTbar']+Histograms['h_STop'],   
        },
        {
            "Name" : "Others",
            "Hist" : Histograms['h_VV']+Histograms['h_QCD'],   
        }
    ]
    return Hists

def F_GroupHistogram_Ratio(Histograms):
    histogram_MC = Histograms['h_WJets'].Clone('h_MC') ; histogram_MC.Add(Histograms['h_TTbar']) ; histogram_MC.Add(Histograms['h_STop']) ; histogram_MC.Add(Histograms['h_VV']) ; histogram_MC.Add(Histograms['h_QCD'])  
    histogram_ratio = construct_ratio_histogram("ratio",Histograms['h_data'],histogram_MC,)
    Hists = [
        {
            "Name" : "data",
            "Hist" : Histograms['h_data'],
        },
        {
            "Name" : "Top",
            "Hist" : Histograms['h_TTbar']+Histograms['h_STop'],   
        },
        {
            "Name" : "W+jets",
            "Hist" : Histograms['h_WJets'],   
        },
        {
            "Name" : "Others",
            "Hist" : Histograms['h_VV']+Histograms['h_QCD'],   
        },
        {
            "Name" : "data/MC",
            "Hist" : histogram_ratio,   
        }
    ]
    return Hists

def F_MergeHistograms_3years(Histograms_dict,region,variable):
    Hist = []
    channels = ["Electron","Muon"]
    years = ["16","17","18"]
    for c in channels:
        for y in years:
            name = '%s%s%s%s'%(c,region,y,variable)
            if name in Histograms_dict:
                Hist.append(Histograms_dict[name])
    Histograms = F_MergeHistograms(Hist,postfix='Lep_161718_%s_%s'%(region,variable))
    return Histograms


def F_Pos_Histogram(Histogram):
    for i in range(1,Histogram.GetNbinsX()+1):
        if Histogram.GetBinContent(i) < 0:
            Histogram.SetBinContent(i,0)
    return Histogram

def F_Store_Histogram(Histograms,OutFile):
    OutFile = os.path.normpath(OutFile)     
    OutputDir = os.path.dirname(OutFile)
    if len(OutputDir) > 0:
        if not os.path.isdir(OutputDir) :
            os.makedirs(OutputDir)     
    outf = ROOT.TFile( OutFile , "recreate")
    for ih in Histograms : 
        Histograms[ih].Write()
    outf.Close()
    print "create %s"%(OutFile)

    
import csv

def F_WriteCsv(data,csvFiles = "tmp.csv"):
    with open(csvFiles, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(data)

def F_HistogramCsv(Histograms,nround = 1,csvFiles="tmp.csv"):
    csv_contents = []
    process = ['Bin']
    for Hist in Histograms:
        process.append(Hist['Name'])
        process.append(Hist['Name']+"err")
    csv_contents.append(process)
    Nbins = Histograms[0]['Hist'].GetNbinsX()
    for b in xrange(1, Nbins+1):
        L = Histograms[0]['Hist'].GetBinLowEdge(b)
        H = Histograms[0]['Hist'].GetBinLowEdge(b+1)
        bin_ = "%s,%s"%(str(int(L)),str(int(H)))
        raws = [str(b)]
        for Hist in Histograms:
            v = round(Hist['Hist'].GetBinContent(b),nround)
            e = round(Hist['Hist'].GetBinError(b),nround)
            raws.append(str(v))
            raws.append(str(e))
        csv_contents.append(raws)
    F_WriteCsv(csv_contents,csvFiles)
    return csvFiles

def F_GetROOFit_Params(Histogram,Func,npol,nI=1000,Random=False):
    HistLabel = "SR1.1.b"
    plotdir = "DVMTF/%s/"%(HistLabel.replace(".","_"))
    h_Ratio,hint,Fit,chi,Ndf,dict_ = F_FitHistogram(Histogram,Func,min=0,max=4000,maxY=1.3,PlotName=plotdir+Func,ColorFunc=2,HistColor=1,HistLabel=HistLabel,NPad=1,Name="Rb",DrawTF=True,confidence_level = 0.997)
    r2 = dict_['r']
    npol_ = int(Func.replace("pol",""))+1
    X = F_GetParams(r2,npol_,nI)
    vars_ = F_RooRealVars(X)
    Vars = []
    for i in range(npol+1):
        if i < len(X):
            Vars.append(vars_[i])
        else:
            Name = "A%s"%(str(i))
            tmp = ROOT.RooRealVar(Name, Name, 0, -10, 10)
            Vars.append(tmp)
    if Random:
        for var in Vars:
            var_ = var.getVal()
            rand = ROOT.TRandom3(1234)
            var_ = var_*rand.Gaus(1.,1)
            var.setVal(var_)
    return Vars

def F_PosDef(result):
    cov = ROOT.TMatrixDSym(result.covarianceMatrix()); 
    eigen = ROOT.TMatrixDSymEigen(cov);
    vectors = eigen.GetEigenVectors();
    values  = eigen.GetEigenValues();
    pos_def = True
    for i in values:
        if i < 0:
            pos_def = False
    return pos_def

def F_ROOFit_Polynomial(Histogram,x,dh,func,Random=False):
    npol = int(func.replace("pol",""))
    Fitpol = func
    vars_ = F_GetROOFit_Params(Histogram,Fitpol,npol,Random=Random)
    pol = F_RooPolynomial(x,vars_)
    result = pol.fitTo(dh, RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit"))
    result = pol.fitTo(dh, RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"), RooFit.Offset(True))
    if F_PosDef(result):
        return pol,result,vars_
    result = pol.fitTo(dh, RooFit.Extended(),RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"), RooFit.Offset(True))
    if F_PosDef(result):
        return pol,result,vars_
    return None,None,None

def F_ROOFit_Polynomial_ntimes(Histogram,x,dh,func,n=10):
    pol,result,vars_ = F_ROOFit_Polynomial(Histogram,x,dh,func)
    if pol:
        return pol,result,vars_
    for i in range(n):
        pol,result,vars_ = F_ROOFit_Polynomial(Histogram,x,dh,func,Random=True)
        if pol :
            return pol,result,vars_      
        
def F_pdf_YList(x,pdf,npoint):
    List = []
    x_min = x.getMin()
    x_max = x.getMax()
    step = float(x_max - x_min)/float(npoint)
    for i in range(npoint+1):
        x_ = x_min+step*i
        x.setVal(x_)
        y_ =  pdf.getVal()
        List.append((x_,y_))
    return List        

def F_TGraphAsymmErrors(Points):
    ErrorBand = ROOT.TGraphAsymmErrors(len(Points))
    XS = [i for i in Points]
    XS.sort()
    for i,X_ in enumerate(XS):
#         print X_,Points[X_][0],Points[X_][2],Points[X_][1]
        errdn = Points[X_][0]-Points[X_][1]
        errup = Points[X_][2]-Points[X_][0]
        ErrorBand.SetPoint(i,X_,Points[X_][0]);
        ErrorBand.SetPointError(i, 0.,0.,errdn,errup)
    return ErrorBand

def F_DecoPDF_ErrorBand(params,x,pdf,npoint,nset):
    Lists = []
    Points = {}
    List = F_pdf_YList(x,pdf,npoint)
    for i in range(npoint+1):
        X_ = List[i][0]
        Points[X_] = []
        Points[X_].append(List[i][1])
    rand = ROOT.TRandom3(1234)
    while len(Lists)<nset :
        for ipar in params:
            val = rand.Gaus(0.,1)
            ipar.setVal( val )
        List = F_pdf_YList(x,pdf,npoint)
        for ipar in params:
            ipar.setVal(0)
        Ys = [i[1] for i in List]
        Ys.sort()
        if Ys[-1] > 0:
            Lists.append(List)
    for i in range(npoint+1):
        Ys = [List[i][1] for List in Lists]
        Xs = [List[i][0] for List in Lists]
        Ys.sort()
        X_ = Xs[0]
        sigma = 0.68 ;
        sigma = 0.95 ;
        sigma = 0.997 ;
        Up = 1-(1-sigma)/2
        Dn = (1-sigma)/2
        Yup = Ys[int(Up*nset)] ; Ydn = Ys[int(Dn*nset)]
        Points[X_].append(Ydn) ; Points[X_].append(Yup)
    return Points


import numpy as np

def fit_results(**kwargs):
    # def fit_results(par_var,ratio,path,filename="fit_parameter.txt"):
    iop   = kwargs.get("iop","")
    ibin  = kwargs.get("ibin","")
    xs = kwargs.get("xs","")
    ys = kwargs.get("ys","")
    pro = kwargs.get("pro","aQGC")
    process = kwargs.get("process","WWW")
    plotpath = kwargs.get("plotpath","/eos/user/q/qiguo/www/VVV/Limit/Plots/parabola/")

    if not os.path.isdir( plotpath ):
        os.makedirs( plotpath ) 
    plotname = "%s/%s_%s_%s.png"%(plotpath,iop,ibin,process)

    ratio = [] ; par_var = []
    x = np.array(xs).astype(np.double)
    y = np.array(ys).astype(np.double)
    gr = ROOT.TGraph(len(x), x, y)
    low = min(xs)
    high = max(xs)
    fitFunc = ROOT.TF1("fit_result","[0]*(x**2) + [1]*x + 1",low,high) 
    fitFunc.SetLineColor(ROOT.kBlue) 
    r = gr.Fit("fit_result","ESR") 
    gr.SetLineWidth(2) 
    gr.SetLineColor(ROOT.kBlue) ;
    gr.SetMarkerStyle(4) 
    gr.SetTitle("Ratio: [%s]/[%s->SM]"%(pro,pro))
    chi2   = r.Chi2() ;
    par0   = fitFunc.GetParameter(0);
    par1   = fitFunc.GetParameter(1);
    err0   = fitFunc.GetParError(0) ;
    err1   = fitFunc.GetParError(1) ;
    c1= ROOT.TCanvas("c1","fitFunc",500,500) 
    c1.SetGridx(1) 
    c1.SetGridy(1) 
    gr.Draw("AP")
    c1.SaveAs(plotname)
    return par0,par1

import ROOT
import os,sys
import math
ROOT.gInterpreter.ProcessLine('#include "Variable.h"')

def UnderOverFlow1D(h):
    Bins=h.GetNbinsX();
    h.SetBinContent( 1,  h.GetBinContent(1)+h.GetBinContent(0) );  h.SetBinError(   1,  math.sqrt( h.GetBinError(1)*h.GetBinError(1) + h.GetBinError(0)*h.GetBinError(0)) );
    h.SetBinContent( Bins,  h.GetBinContent(Bins)+h.GetBinContent(Bins+1) );  h.SetBinError(   Bins,  math.sqrt( h.GetBinError(Bins)*h.GetBinError(Bins) + h.GetBinError(Bins+1)*h.GetBinError(Bins+1)) );
    return h;

from optparse   import OptionParser    
parser = OptionParser()
parser.add_option('--P'   , action="store",type="string",dest="Process"    ,default=None)
parser.add_option('--UN'   , action="store",type="string",dest="UN"    ,default=None)
(opts, args) = parser.parse_args()

input_processes = [opts.Process]
UNname = opts.UN

# ======================================= configurations =================================
# ======================================= configurations =================================
# ======================================= configurations =================================
# ======================================= configurations =================================

XSs = {
    "WWW_Dim6": "764.624328",
    "WWZ_Dim6": "763.260128",
    "WZZ_Dim6": "337.473507",
    "ZZZ_Dim6": "75.784416",
    'WWW_Dim8': '630.1894555568599', 
    'WWZ_Dim8': '265.93670603964983', 
    'WZZ_Dim8': '177.04690209086584', 
    'ZZZ_Dim8': '46.407877176107704',
    "WWW_Dim6_old": "764.624328",
    "WWZ_Dim6_old": "763.260128",
    "WZZ_Dim6_old": "337.473507",
    "ZZZ_Dim6_old": "75.784416",
    "WWW_Dim6_old2": "764.624328",
    "WWZ_Dim6_old2": "763.260128",
    "WZZ_Dim6_old2": "337.473507",
    "ZZZ_Dim6_old2": "75.784416",
}

ROOTFiles = {
    "WWW_Dim6": "WWW_1Jet_xqcut15_Dim6_cW_cHd_cHWB_cHW_4F-RunIISummer20UL18NanoAODv9_All.root",
    "WWZ_Dim6": "WWZ_1Jet_xqcut15_Dim6_WWZ_1Jet_xqcut15_12Operators_4F_NoFilter_ReSubmit-RunIISummer20UL18NanoAODv9_All.root",
    "WZZ_Dim6": "WZZ_1J.root",
    "ZZZ_Dim6": "ZZZ_1Jet_xqcut15_12Operators_4F_NoFilter-RunIISummer20UL18NanoAODv9_All.root",
    
    "WWW_Dim8": "WWW_All.root",
    "WWZ_Dim8": "WWZ_All.root",
    "WZZ_Dim8": "WZZ_All.root",
    "ZZZ_Dim8": "ZZZ_All.root",

    "WWW_Dim6_old": "Dim6_WWW_1Jet_xqcut15_NoFilter.root",
    "WWZ_Dim6_old": "Dim6_WWZ_1J.root",
    "WZZ_Dim6_old": "Dim6_WZZ_1J.root",
    "ZZZ_Dim6_old": "Dim6_ZZZ_1J.root",
    
}


MuSR = "((abs(Lep1fatJet2_LeptonPDGID) == 13 && Lep1fatJet2_Muon_pfRelIso04_all < 0.1 )|| (abs(Lep1fatJet2_LeptonPDGID) == 11))  && Lep1fatJet2_LeptonPt > 30 && Lep1fatJet2_FatJet_pt > 200 && Lep1fatJet2_FatJet_pt_2 > 200 && Lep1fatJet2_FatJet_msoftdrop > 40 && Lep1fatJet2_FatJet_msoftdrop_2 > 40 && Lep1fatJet2_Muon_pfRelIso04_all < 0.1 && Lep1fatJet2_LeptonPt > 55 && Lep1fatJet2_MET_pt > 50 && Lep1fatJet2_FatJet_msoftdrop > 65 && Lep1fatJet2_FatJet_msoftdrop < 110 && Lep1fatJet2_FatJet_msoftdrop_2 > 65 && Lep1fatJet2_FatJet_msoftdrop_2 < 110 && Common_nb_tight == 0 && Lep1fatJet2_FatJet_particleNetMD_WvsQCD > 0.64 && Lep1fatJet2_FatJet_particleNetMD_WvsQCD_2 > 0.64 && Lep1fatJet2_LeptonicWPt > 150 "

MuSR_jerup = "((abs(Lep1fatJet2_LeptonPDGID) == 13 && Lep1fatJet2_Muon_pfRelIso04_all < 0.1 )|| (abs(Lep1fatJet2_LeptonPDGID) == 11))  && Lep1fatJet2_LeptonPt > 30 && Lep1fatJet2_FatJet_pt_jerup > 200 && Lep1fatJet2_FatJet_pt_jerup_2 > 200 && Lep1fatJet2_FatJet_msoftdrop_jerup > 40 && Lep1fatJet2_FatJet_msoftdrop_jerup_2 > 40 && Lep1fatJet2_Muon_pfRelIso04_all < 0.1 && Lep1fatJet2_LeptonPt > 55 && Lep1fatJet2_MET_pt > 50 && Lep1fatJet2_FatJet_msoftdrop_jerup > 65 && Lep1fatJet2_FatJet_msoftdrop_jerup < 110 && Lep1fatJet2_FatJet_msoftdrop_jerup_2 > 65 && Lep1fatJet2_FatJet_msoftdrop_jerup_2 < 110 && Common_nb_tight == 0 && Lep1fatJet2_FatJet_particleNetMD_WvsQCD > 0.64 && Lep1fatJet2_FatJet_particleNetMD_WvsQCD_2 > 0.64 && Lep1fatJet2_LeptonicWPt > 150 "

MuSR_jerdn = "((abs(Lep1fatJet2_LeptonPDGID) == 13 && Lep1fatJet2_Muon_pfRelIso04_all < 0.1 )|| (abs(Lep1fatJet2_LeptonPDGID) == 11))  && Lep1fatJet2_LeptonPt > 30 && Lep1fatJet2_FatJet_pt_jerdn > 200 && Lep1fatJet2_FatJet_pt_jerdn_2 > 200 && Lep1fatJet2_FatJet_msoftdrop_jerdn > 40 && Lep1fatJet2_FatJet_msoftdrop_jerdn_2 > 40 && Lep1fatJet2_Muon_pfRelIso04_all < 0.1 && Lep1fatJet2_LeptonPt > 55 && Lep1fatJet2_MET_pt > 50 && Lep1fatJet2_FatJet_msoftdrop_jerdn > 65 && Lep1fatJet2_FatJet_msoftdrop_jerdn < 110 && Lep1fatJet2_FatJet_msoftdrop_jerdn_2 > 65 && Lep1fatJet2_FatJet_msoftdrop_jerdn_2 < 110 && Common_nb_tight == 0 && Lep1fatJet2_FatJet_particleNetMD_WvsQCD > 0.64 && Lep1fatJet2_FatJet_particleNetMD_WvsQCD_2 > 0.64 && Lep1fatJet2_LeptonicWPt > 150 "

MuSR_jesup = "((abs(Lep1fatJet2_LeptonPDGID) == 13 && Lep1fatJet2_Muon_pfRelIso04_all < 0.1 )|| (abs(Lep1fatJet2_LeptonPDGID) == 11))  && Lep1fatJet2_LeptonPt > 30 && Lep1fatJet2_FatJet_pt_jesup > 200 && Lep1fatJet2_FatJet_pt_jesup_2 > 200 && Lep1fatJet2_FatJet_msoftdrop_jesup > 40 && Lep1fatJet2_FatJet_msoftdrop_jesup_2 > 40 && Lep1fatJet2_Muon_pfRelIso04_all < 0.1 && Lep1fatJet2_LeptonPt > 55 && Lep1fatJet2_MET_pt > 50 && Lep1fatJet2_FatJet_msoftdrop_jesup > 65 && Lep1fatJet2_FatJet_msoftdrop_jesup < 110 && Lep1fatJet2_FatJet_msoftdrop_jesup_2 > 65 && Lep1fatJet2_FatJet_msoftdrop_jesup_2 < 110 && Common_nb_tight == 0 && Lep1fatJet2_FatJet_particleNetMD_WvsQCD > 0.64 && Lep1fatJet2_FatJet_particleNetMD_WvsQCD_2 > 0.64 && Lep1fatJet2_LeptonicWPt > 150 "

MuSR_jesdn = "((abs(Lep1fatJet2_LeptonPDGID) == 13 && Lep1fatJet2_Muon_pfRelIso04_all < 0.1 )|| (abs(Lep1fatJet2_LeptonPDGID) == 11))  && Lep1fatJet2_LeptonPt > 30 && Lep1fatJet2_FatJet_pt_jesdn > 200 && Lep1fatJet2_FatJet_pt_jesdn_2 > 200 && Lep1fatJet2_FatJet_msoftdrop_jesdn > 40 && Lep1fatJet2_FatJet_msoftdrop_jesdn_2 > 40 && Lep1fatJet2_Muon_pfRelIso04_all < 0.1 && Lep1fatJet2_LeptonPt > 55 && Lep1fatJet2_MET_pt > 50 && Lep1fatJet2_FatJet_msoftdrop_jesdn > 65 && Lep1fatJet2_FatJet_msoftdrop_jesdn < 110 && Lep1fatJet2_FatJet_msoftdrop_jesdn_2 > 65 && Lep1fatJet2_FatJet_msoftdrop_jesdn_2 < 110 && Common_nb_tight == 0 && Lep1fatJet2_FatJet_particleNetMD_WvsQCD > 0.64 && Lep1fatJet2_FatJet_particleNetMD_WvsQCD_2 > 0.64 && Lep1fatJet2_LeptonicWPt > 150 "

weight_center_expr = 'return (%s) * Common_event_lepSF * Common_event_mediumBtagSF * Common_event_prefireWeight * Common_event_puWeight;'

weight_BtagUp_expr = 'return (%s) * Common_event_lepSF * Common_event_mediumBtagSFup * Common_event_prefireWeight * Common_event_puWeight;'
weight_BtagDn_expr = 'return (%s) * Common_event_lepSF * Common_event_mediumBtagSFdn * Common_event_prefireWeight * Common_event_puWeight;'
weight_lepelUp_expr = 'return (%s) * Common_event_lepSFelup * Common_event_mediumBtagSF * Common_event_prefireWeight * Common_event_puWeight;'
weight_lepelDn_expr = 'return (%s) * Common_event_lepSFeldn * Common_event_mediumBtagSF * Common_event_prefireWeight * Common_event_puWeight;'
weight_lepmuUp_expr = 'return (%s) * Common_event_lepSFmuup * Common_event_mediumBtagSF * Common_event_prefireWeight * Common_event_puWeight;'
weight_lepmuDn_expr = 'return (%s) * Common_event_lepSFmudn * Common_event_mediumBtagSF * Common_event_prefireWeight * Common_event_puWeight;'

weight_PUUp_expr = 'return (%s) * Common_event_lepSF * Common_event_mediumBtagSF * Common_event_prefireWeight * Common_event_puWeightup;'
weight_PUDn_expr = 'return (%s) * Common_event_lepSF * Common_event_mediumBtagSF * Common_event_prefireWeight * Common_event_puWeightdn;'

# ======================================= parameters =================================

if UNname == "jer":
    Bins = [
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_center_expr,""),
        (("Plot",";x;y",35,500,4000),"MJJlv_jerup",MuSR_jerup,weight_center_expr,""),
        (("Plot",";x;y",35,500,4000),"MJJlv_jerdn",MuSR_jerdn,weight_center_expr,""),
    ]

if UNname == "jes":
    Bins = [
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_center_expr,""),
        (("Plot",";x;y",35,500,4000),"MJJlv_jesup",MuSR_jesup,weight_center_expr,""),
        (("Plot",";x;y",35,500,4000),"MJJlv_jesdn",MuSR_jesdn,weight_center_expr,""),
    ]

if UNname == "BtagSF":
    Bins = [
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_center_expr,""),
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_BtagUp_expr,"_BTagSFUp"),
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_BtagDn_expr,"_BTagSFDn"),
    ]

if UNname == "lepelSF":
    Bins = [
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_center_expr,""),
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_lepelUp_expr,"_lepelSFUp"),
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_lepelDn_expr,"_lepelSFDn"),
    ]

if UNname == "lepmuSF":
    Bins = [
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_center_expr,""),
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_lepmuUp_expr,"_lepmuSFUp"),
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_lepmuDn_expr,"_lepmuSFDn"),
    ]

if UNname == "prefire":
    Bins = [
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_center_expr,""),
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_center_expr,"_prefireUp"),
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_center_expr,"_prefireDn"),
    ]

if UNname == "trigger":
    Bins = [
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_center_expr,""),
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_center_expr,"_triggerUp"),
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_center_expr,"_triggerDn"),
    ]

if UNname == "pileup":
    Bins = [
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_center_expr,""),
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_PUUp_expr,"_pileupUp"),
        (("Plot",";x;y",35,500,4000),"MJJlv",MuSR,weight_PUDn_expr,"_pileupDn"),
    ]

GetParameter_Range_Inst = GetParameter_Range.GetParameter_Range()

processes = input_processes
# for i in ROOTFiles:
#     processes.append(ROOTFiles[i])
processes_Add = []

# dim6_OPs = ["SM","cW",]
# dim8_OPs = ["SM","FT0",]
dim6_OPs = ["SM","cW","cHbox","cHDD","cHW","cHB","cHWB","cHl3","cHq1","cHq3","cHu","cHd","cll1"]
dim8_OPs = ['FT1','FT0','FT3','FT2','FT5','FT4','FT7','FT6','FM6','FM7','FM4','FM5','FM2','FM3','FM0','FM1','FS0','FS1','FS2','SM','FT8','FT9']

dim6_BSM_OPs = [i for i in dim6_OPs if i != "SM"]
dim8_BSM_OPs = [i for i in dim8_OPs if i != "SM"]

Path = "/eos/user/q/qiguo/ROOTFILE/VVVEFT/Ntuple/UL18/V1/"
LUMI = 138.

# =================================================================
# =================================================================
# =================================================================
# =================================================================
# =================================================================

# import ipdb

countindex = 0
Histos = {}
for Bin,Var,Cut,weight_expr,postfix in Bins:
    for process in processes : 
        File = ROOTFiles[process]
        File = Path+File
        rdf = ROOT.RDataFrame("t",File)
        Parameter_Range,Weight_Index = GetParameter_Range_Inst.Parameter_Range(process)
        
        if "Dim6" in process:
            OPs = dim6_OPs
        if "Dim8" in process:
            OPs = dim8_OPs

        XS = float(XSs[process])
        filein = ROOT.TFile( File )
        Nevents = int(filein.Get("Wgt__h_nevents").GetBinContent(1))
        lumiWeight = (XS * LUMI) / float(abs(Nevents));
        filein.Close()
        
        for OP in OPs:
            if not OP in Parameter_Range:
                continue
            for grid in Parameter_Range[OP]:
                Variables = [
                    {"Lep1fatJet2_LeptonicWPt_jerup":'''
                    TLorentzVector NeutrinoP4,LeptonicWP4,lep;
                    lep.SetPtEtaPhiE(Lep1fatJet2_LeptonPt,Lep1fatJet2_LeptonEta,Lep1fatJet2_LeptonPhi,Lep1fatJet2_LeptonE);
                    NeutrinoP4 = getNeutrinoP4(Lep1fatJet2_MET_pt_jerup, Lep1fatJet2_MET_phi_jerup, lep, 1);
                    LeptonicWP4 = NeutrinoP4+lep;
                    return LeptonicWP4.Pt();
                    '''},
                    {"Lep1fatJet2_LeptonicWPt_jerdn":'''
                    TLorentzVector NeutrinoP4,LeptonicWP4,lep;
                    lep.SetPtEtaPhiE(Lep1fatJet2_LeptonPt,Lep1fatJet2_LeptonEta,Lep1fatJet2_LeptonPhi,Lep1fatJet2_LeptonE);
                    NeutrinoP4 = getNeutrinoP4(Lep1fatJet2_MET_pt_jerdn, Lep1fatJet2_MET_phi_jerdn, lep, 1);
                    LeptonicWP4 = NeutrinoP4+lep;
                    return LeptonicWP4.Pt();
                    '''},
                    {"Lep1fatJet2_LeptonicWPt_jesup":'''
                    TLorentzVector NeutrinoP4,LeptonicWP4,lep;
                    lep.SetPtEtaPhiE(Lep1fatJet2_LeptonPt,Lep1fatJet2_LeptonEta,Lep1fatJet2_LeptonPhi,Lep1fatJet2_LeptonE);
                    NeutrinoP4 = getNeutrinoP4(Lep1fatJet2_MET_pt_jesup, Lep1fatJet2_MET_phi_jesup, lep, 1);
                    LeptonicWP4 = NeutrinoP4+lep;
                    return LeptonicWP4.Pt();
                    '''},
                    {"Lep1fatJet2_LeptonicWPt_jesdn":'''
                    TLorentzVector NeutrinoP4,LeptonicWP4,lep;
                    lep.SetPtEtaPhiE(Lep1fatJet2_LeptonPt,Lep1fatJet2_LeptonEta,Lep1fatJet2_LeptonPhi,Lep1fatJet2_LeptonE);
                    NeutrinoP4 = getNeutrinoP4(Lep1fatJet2_MET_pt_jesdn, Lep1fatJet2_MET_phi_jesdn, lep, 1);
                    LeptonicWP4 = NeutrinoP4+lep;
                    return LeptonicWP4.Pt();
                    '''},
                    {"MJJlv":'''
                TLorentzVector AK8_1, AK8_2, leptonicW; 
                AK8_1.SetPtEtaPhiM( Lep1fatJet2_FatJet_pt, Lep1fatJet2_FatJet_eta , Lep1fatJet2_FatJet_phi, Lep1fatJet2_FatJet_msoftdrop ); 
                AK8_2.SetPtEtaPhiM( Lep1fatJet2_FatJet_pt_2, Lep1fatJet2_FatJet_eta_2 , Lep1fatJet2_FatJet_phi_2, Lep1fatJet2_FatJet_msoftdrop_2 ); 
                leptonicW.SetPtEtaPhiM( Lep1fatJet2_LeptonicWPt, Lep1fatJet2_LeptonicWEta , Lep1fatJet2_LeptonicWPhi, Lep1fatJet2_LeptonicWM ); 
                return ( AK8_1 + AK8_2 + leptonicW ).M();
                    '''},
                    {"MJJlv_jerup":'''
                TLorentzVector AK8_1, AK8_2, leptonicW; 
                AK8_1.SetPtEtaPhiM( Lep1fatJet2_FatJet_pt_jerup, Lep1fatJet2_FatJet_eta , Lep1fatJet2_FatJet_phi, Lep1fatJet2_FatJet_msoftdrop_jerup ); 
                AK8_2.SetPtEtaPhiM( Lep1fatJet2_FatJet_pt_jerup_2, Lep1fatJet2_FatJet_eta_2 , Lep1fatJet2_FatJet_phi_2, Lep1fatJet2_FatJet_msoftdrop_jerup_2 ); 
                TLorentzVector NeutrinoP4,lep;
                lep.SetPtEtaPhiE(Lep1fatJet2_LeptonPt,Lep1fatJet2_LeptonEta,Lep1fatJet2_LeptonPhi,Lep1fatJet2_LeptonE);
                NeutrinoP4 = getNeutrinoP4(Lep1fatJet2_MET_pt_jerup, Lep1fatJet2_MET_phi_jerup, lep, 1);
                leptonicW = NeutrinoP4+lep;
                return ( AK8_1 + AK8_2 + leptonicW ).M();
                    '''},
                    {"MJJlv_jerdn":'''
                TLorentzVector AK8_1, AK8_2, leptonicW; 
                AK8_1.SetPtEtaPhiM( Lep1fatJet2_FatJet_pt_jerdn, Lep1fatJet2_FatJet_eta , Lep1fatJet2_FatJet_phi, Lep1fatJet2_FatJet_msoftdrop_jerdn ); 
                AK8_2.SetPtEtaPhiM( Lep1fatJet2_FatJet_pt_jerdn_2, Lep1fatJet2_FatJet_eta_2 , Lep1fatJet2_FatJet_phi_2, Lep1fatJet2_FatJet_msoftdrop_jerdn_2 ); 
                TLorentzVector NeutrinoP4,lep;
                lep.SetPtEtaPhiE(Lep1fatJet2_LeptonPt,Lep1fatJet2_LeptonEta,Lep1fatJet2_LeptonPhi,Lep1fatJet2_LeptonE);
                NeutrinoP4 = getNeutrinoP4(Lep1fatJet2_MET_pt_jerdn, Lep1fatJet2_MET_phi_jerdn, lep, 1);
                leptonicW = NeutrinoP4+lep;
                return ( AK8_1 + AK8_2 + leptonicW ).M();
                    '''},
                    {"MJJlv_jesup":'''
                TLorentzVector AK8_1, AK8_2, leptonicW; 
                AK8_1.SetPtEtaPhiM( Lep1fatJet2_FatJet_pt_jesup, Lep1fatJet2_FatJet_eta , Lep1fatJet2_FatJet_phi, Lep1fatJet2_FatJet_msoftdrop_jesup ); 
                AK8_2.SetPtEtaPhiM( Lep1fatJet2_FatJet_pt_jesup_2, Lep1fatJet2_FatJet_eta_2 , Lep1fatJet2_FatJet_phi_2, Lep1fatJet2_FatJet_msoftdrop_jesup_2 ); 
                TLorentzVector NeutrinoP4,lep;
                lep.SetPtEtaPhiE(Lep1fatJet2_LeptonPt,Lep1fatJet2_LeptonEta,Lep1fatJet2_LeptonPhi,Lep1fatJet2_LeptonE);
                NeutrinoP4 = getNeutrinoP4(Lep1fatJet2_MET_pt_jesup, Lep1fatJet2_MET_phi_jesup, lep, 1);
                leptonicW = NeutrinoP4+lep;
                return ( AK8_1 + AK8_2 + leptonicW ).M();
                    '''},
                    {"MJJlv_jesdn":'''
                TLorentzVector AK8_1, AK8_2, leptonicW; 
                AK8_1.SetPtEtaPhiM( Lep1fatJet2_FatJet_pt_jesdn, Lep1fatJet2_FatJet_eta , Lep1fatJet2_FatJet_phi, Lep1fatJet2_FatJet_msoftdrop_jesdn ); 
                AK8_2.SetPtEtaPhiM( Lep1fatJet2_FatJet_pt_jesdn_2, Lep1fatJet2_FatJet_eta_2 , Lep1fatJet2_FatJet_phi_2, Lep1fatJet2_FatJet_msoftdrop_jesdn_2 ); 
                TLorentzVector NeutrinoP4,lep;
                lep.SetPtEtaPhiE(Lep1fatJet2_LeptonPt,Lep1fatJet2_LeptonEta,Lep1fatJet2_LeptonPhi,Lep1fatJet2_LeptonE);
                NeutrinoP4 = getNeutrinoP4(Lep1fatJet2_MET_pt_jesdn, Lep1fatJet2_MET_phi_jesdn, lep, 1);
                leptonicW = NeutrinoP4+lep;
                return ( AK8_1 + AK8_2 + leptonicW ).M();
                    '''},
                    {"Lep1fatJet2_FatJet_particleNetMD_WvsQCD":'''
                return (Lep1fatJet2_FatJet_particleNetMD_Xcc+Lep1fatJet2_FatJet_particleNetMD_Xqq)/(Lep1fatJet2_FatJet_particleNetMD_Xcc+Lep1fatJet2_FatJet_particleNetMD_Xqq+Lep1fatJet2_FatJet_particleNetMD_QCD);
                    ''',
                    },
                    {"Lep1fatJet2_FatJet_particleNetMD_WvsQCD_2":'''
                return (Lep1fatJet2_FatJet_particleNetMD_Xcc_2+Lep1fatJet2_FatJet_particleNetMD_Xqq_2)/(Lep1fatJet2_FatJet_particleNetMD_Xcc_2+Lep1fatJet2_FatJet_particleNetMD_Xqq_2+Lep1fatJet2_FatJet_particleNetMD_QCD_2);
                    '''
                    },
                ]

                weight = weight_expr%(lumiWeight)

                Variables.append({"weight":weight})

                reweightname = "weight_%s_%s"%(OP,str(grid).replace("-","m").replace(".","p"))
                gridid = 0
                for index,i in enumerate(Parameter_Range[OP]):
                    if i == grid:
                        gridid = Weight_Index[OP][0]+index
                        break

                reweight = '''
                double weight_ = 0;
                if(Common_LHEWeight_mg_reweighting.size()>0){
                    return weight*Common_LHEWeight_mg_reweighting[%s];
                }
                if(Common_LHEReweightingWeight.size()>0){
                    return weight*Common_LHEReweightingWeight[%s];
                }
                cout <<weight_<<endl;
                return weight_;
                '''%(str(gridid),str(gridid))
                Variables.append({reweightname:reweight})

                weightname = reweightname

                for ivar in Variables:
                    for var in ivar:
                        if var not in [str(i) for i in rdf.GetColumnNames()]:
                            rdf = rdf.Define(var,ivar[var])
                rdf = rdf.Filter(Cut)
                HistName = "h_%s_%s_%s_%s"%(process,OP,str(grid).replace("-","m").replace(".","p"),Var+postfix)
                Histogram = rdf.Histo1D(Bin,Var,weightname)
                # Histogram = UnderOverFlow1D(Histogram)
                Histos[HistName] = Histogram
                # Histos[HistName] = Histos[HistName]
                # Histos[HistName] = UnderOverFlow1D(Histos[HistName])
                countindex += 1
                if countindex%100 == 0:
                    print process,countindex,"done"

varbins=[500,1600,2600,3400,4000]
histo_bin = array('d',varbins)
for h in Histos:
    Histos[h] = Histos[h].GetValue().Clone(h)
    Histos[h] = Histos[h].Rebin(len(histo_bin)-1,h,histo_bin)
    Histos[h] = UnderOverFlow1D(Histos[h])
    print h,Histos[h].Integral(),[ Histos[h].GetBinContent(i) for i in range(1,Histos[h].GetNbinsX()+1)]


# =================================================================
# =================================================================
# =================================================================
# =================================================================
# =================================================================

import re

ParamHists = {}


for process in processes:
    pro = "dim"
    ibin = 4

    Parameters = {}
    
    for Bin,Var,Cut,weight_expr,postfix in Bins:
        if "Dim6" in process:
            BSM_OPs = dim6_BSM_OPs
        if "Dim8" in process:
            BSM_OPs = dim8_BSM_OPs
        Parameter_Range,Weight_Index = GetParameter_Range_Inst.Parameter_Range(process)
        for iop in BSM_OPs:
            if iop not in Parameter_Range:
                continue
            sm_name = "h_%s_SM_0_%s"%(process,Var+postfix)
            hsmname = "hsm_%s_%s"%(process,Var+postfix)
            hsm = Histos[sm_name].Clone(hsmname)
            # hsm = hsm.Rebin(len(histo_bin)-1,hsmname,histo_bin)

            hlname = "hl_%s_%s"%(process,Var+postfix)
            hl = hsm.Clone(hlname)

            hqname = "hq_%s_%s"%(process,Var+postfix)
            hq = hsm.Clone(hqname)

            for ibin in range(1,len(varbins)):
                values = {}
                for hn in Histos:
                    ss = r'h_%s_%s_(.*)_%s$'%(process,iop,Var+postfix)
                    if re.search(ss,hn):
                        param = float(re.search(ss,hn).group(1).replace("m","-").replace("p","."))
                        if abs(param) < 0.00001:
                            param = param*(10**12)
                        h_tmp = Histos[hn].Clone("h_tmp")
                        # h_tmp = h_tmp.Rebin(len(histo_bin)-1,"htmp",histo_bin)
                        if hsm.GetBinContent(ibin) == 0:
                            y_ = 0
                        else:
                            y_ = h_tmp.GetBinContent(ibin)/hsm.GetBinContent(ibin)
                        values[param] = y_
                params = [i for i in values]
                params.sort()
                xs = []
                ys = []
                for ip in params:
                    xs.append(ip)
                    ys.append(values[ip])
                print iop,Var,postfix,"bin",ibin, "==========================="
                print xs
                print ys
                par0,par1 = fit_results(iop=iop,ibin=str(ibin),xs=xs,ys=ys,pro=pro,process=process+"_"+Var+postfix)
                print par0,par1
                Parameters[(Var+postfix,process,ibin)] = (par0,par1,hsm.GetBinContent(ibin))

                hq.SetBinContent(ibin,par0*hsm.GetBinContent(ibin))
                hl.SetBinContent(ibin,(1+par1+par0)*hsm.GetBinContent(ibin))

            ParamHists[(Var+postfix,process,iop,"sm")] = hsm
            ParamHists[(Var+postfix,process,iop,"l")] = hl
            ParamHists[(Var+postfix,process,iop,"q")] = hq
            
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================

ROOTFileDir = "ROOTFiles"
if not os.path.isdir(ROOTFileDir):
    os.makedirs(ROOTFileDir)

Hists_EFT = {}
for process in processes:
    if "Dim6" in process:
        BSM_OPs = dim6_BSM_OPs
    if "Dim8" in process:
        BSM_OPs = dim8_BSM_OPs
    Parameter_Range,Weight_Index = GetParameter_Range_Inst.Parameter_Range(process)
    for iop in BSM_OPs:
        if iop not in Parameter_Range:
            continue
        Hists_To_Store = {}
        InterFile = "%s/%s_%s_%s.root"%(ROOTFileDir,process,iop,UNname)

        for Bin,Var,Cut,weight_expr,postfix in Bins:
            if ("up" in Var) or ("dn" in Var):
                namepostfix = Var.replace('MJJlv','')
            else:
                namepostfix = postfix
            Hists_EFT["%s_%s_slq_%s"%(process,iop,namepostfix)] = ParamHists[(Var+postfix, process, iop, 'l')].Clone("h_sm_lin_quad_%s_%s%s"%(process,iop,namepostfix))
            Hists_EFT["%s_%s_q_%s"%(process,iop,namepostfix)]   = ParamHists[(Var+postfix, process, iop, 'q')].Clone("h_quad_%s_%s%s"%(process,iop,namepostfix))
            Hists_EFT["%s_%s_s_%s"%(process,iop,namepostfix)]   = ParamHists[(Var+postfix, process, iop, 'sm')].Clone("h_sm_%s_%s%s"%(process,iop,namepostfix))

            Hists_To_Store["%s_%s_slq_%s"%(process,iop,namepostfix)] = ParamHists[(Var+postfix, process, iop, 'l')].Clone("h_sm_lin_quad_%s_%s%s"%(process,iop,namepostfix))
            Hists_To_Store["%s_%s_q_%s"%(process,iop,namepostfix)]   = ParamHists[(Var+postfix, process, iop, 'q')].Clone("h_quad_%s_%s%s"%(process,iop,namepostfix))
            Hists_To_Store["%s_%s_s_%s"%(process,iop,namepostfix)]   = ParamHists[(Var+postfix, process, iop, 'sm')].Clone("h_sm_%s_%s%s"%(process,iop,namepostfix))

        for h in Histos:
            if iop in h:
                Hists_To_Store[h] = Histos[h]
        F_Store_Histogram(Hists_To_Store,InterFile)




#     Hists_EFT["VVV_%s_slq"%(iop)] = Hists_EFT["%s_%s_slq"%(processes[0],iop)].Clone("h_sm_lin_quad_VVV_%s"%(iop))
#     Hists_EFT["VVV_%s_q"%(iop)] = Hists_EFT["%s_%s_q"%(processes[0],iop)].Clone("h_quad_VVV_%s"%(iop))
#     Hists_EFT["VVV_%s_s"%(iop)] = Hists_EFT["%s_%s_s"%(processes[0],iop)].Clone("h_sm_VVV_%s"%(iop))

#     for process in processes[1:]:
#         Hists_EFT["VVV_%s_slq"%(iop)].Add(Hists_EFT["%s_%s_slq"%(process,iop)])
#         Hists_EFT["VVV_%s_q"%(iop)].Add(Hists_EFT["%s_%s_q"%(process,iop)])
#         Hists_EFT["VVV_%s_s"%(iop)].Add(Hists_EFT["%s_%s_s"%(process,iop)])

for i in Hists_EFT:
    print i,"========",Hists_EFT[i].GetName(), Hists_EFT[i].Integral()