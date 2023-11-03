import ROOT
import os

from ROOT import gROOT, TPaveLabel, TPie, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack, TGraph, TGraphErrors,TChain,TArrow, TCanvas, TMatrixDSym, TMath, TText, TPad, TVectorD, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHist,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite

def F_FitHistogram(histogram,function,min,max,PlotName,**kwargs):
    ColorFunc = kwargs.get("ColorFunc",1)
    HistColor = kwargs.get("HistColor",1)
    HistLabel = kwargs.get("HistLabel",1)
    NPad = kwargs.get("NPad",1)
    YLabel = kwargs.get("YLabel","TF")
    DrawTF = kwargs.get("DrawTF",False)
    Name = kwargs.get("Name","")
    confidence_level= kwargs.get("confidence_level",0.95)
    maxY= kwargs.get("maxY",3)
    AddFitFunction= kwargs.get("AddFitFunction",None)
    AddHistogram= kwargs.get("AddHistogram",None)
    nbinsDraw = kwargs.get("nbinsDraw",None)
    
    histogram.SetLineColor(HistColor)
    histogram.SetLineWidth(3)
    histogram.SetLineStyle(7)
    histogram.SetMarkerColor(HistColor)
    histogram.SetMarkerStyle(8);
    histogram.SetMarkerSize(0);
    histogram.SetStats(0)


    polfunc = ROOT.TF1("polfunc"+Name, function, min, max)
    
    r = ROOT.TFitResultPtr() ; 
    r = histogram.Fit(polfunc,"S,N,0")
    r.NormalizeErrors();
    values = r.GetConfidenceIntervals(confidence_level, True)
    print values

#     print "Chi2/Ndf = %0.3f"%(r.Chi2()/r.Ndf())
    polfunc.SetLineColor(ColorFunc);
    polfunc.SetLineWidth(4)
    polfunc.SetLineStyle(1)

    ptstats = ROOT.TPaveStats(0.61,0.75,0.97,0.9,"brNDC");
    ptstats.SetBorderSize(1);
    ptstats.SetFillColor(0);
    #ptstats.SetTextAlign(12);
    ptstats.SetTextFont(42);
    ptstats.SetTextSize(0.025);
    chiText = "#chi^{2}/ndf=%0.3f/%s(%0.3f)"%(r.Chi2(),r.Ndf(),r.Chi2()/r.Ndf())
    chi = r.Chi2()/r.Ndf()
    ptstats.AddText(chiText)

    ptstats.AddText("p0       = %0.1f #pm %0.1f "%(r.Parameter(0),r.ParError(0)));
    ptstats.AddText("p1       = %0.1e #pm %0.1e "%(r.Parameter(1),r.ParError(1)));
    if function in ["pol2","pol3","pol4","pol5","pol6"]: ptstats.AddText("p2       = %0.1e #pm %0.1e "%(r.Parameter(2),r.ParError(2)));
    if function in [       "pol3","pol4","pol5","pol6"]: ptstats.AddText("p3       = %0.1e #pm %0.1e "%(r.Parameter(3),r.ParError(3)));
    if function in [              "pol4","pol5","pol6"]: ptstats.AddText("p4       = %0.1e #pm %0.1e "%(r.Parameter(4),r.ParError(4)));
    if function in [                     "pol5","pol6"]: ptstats.AddText("p5       = %0.1e #pm %0.1e "%(r.Parameter(5),r.ParError(5)));
    if function in [                            "pol6"]: ptstats.AddText("p6       = %0.1e #pm %0.1e "%(r.Parameter(6),r.ParError(6)));
    if function in [                            "pol7"]: ptstats.AddText("p7       = %0.1e #pm %0.1e "%(r.Parameter(7),r.ParError(7)));
    if function in [                            "pol8"]: ptstats.AddText("p8       = %0.1e #pm %0.1e "%(r.Parameter(9),r.ParError(8)));
    ptstats.SetOptStat(0);
    ptstats.SetOptFit(111);

    if not nbinsDraw:
        nbinsDraw = 100*histogram.GetNbinsX()
    hintUpDraw  = TH1D("hintUp"+Name,"hintUp"+Name+"more_bin",nbinsDraw,histogram.GetBinLowEdge(1),histogram.GetBinLowEdge(histogram.GetNbinsX()+1)); hintUpDraw.Sumw2();
    hintDnDraw  = TH1D("hintDn"+Name,"hintDn"+Name+"more_bin",nbinsDraw,histogram.GetBinLowEdge(1),histogram.GetBinLowEdge(histogram.GetNbinsX()+1)); hintDnDraw.Sumw2();
    
    hint = histogram.Clone("hint"+Name)
    hintUp = histogram.Clone("hintUp"+Name)
    hintDn = histogram.Clone("hintDn"+Name)
    for i in range(hint.GetNbinsX()):
        hint.SetBinContent(i+1,polfunc.Eval(hint.GetBinCenter(i+1)))
        hintUp.SetBinContent(i+1,polfunc.Eval(hint.GetBinCenter(i+1))+values[i])
        hintDn.SetBinContent(i+1,polfunc.Eval(hint.GetBinCenter(i+1))-values[i])
        hint.SetBinError(i+1,values[i])
        
    for i in range(hintUpDraw.GetNbinsX()):
        hintUpDraw.SetBinContent(i+1,hintUp.Interpolate(hintUpDraw.GetBinCenter(i+1)))
        hintDnDraw.SetBinContent(i+1,hintDn.Interpolate(hintDnDraw.GetBinCenter(i+1)))
        

    hint.SetMarkerSize(0)
    hint.SetLineWidth(0)
    hint.SetStats(0);
    hint.SetFillColor(polfunc.GetLineColor());
    hint.SetFillStyle(3003)
    hintUpDraw.SetLineColorAlpha(ROOT.kGreen, 0.45)
    hintDnDraw.SetLineColorAlpha(ROOT.kBlue, 0.45)
    
    
    NPad = kwargs.get('NPad',2)
    LabelHistogram = kwargs.get('LabelHistogram')
    
    if DrawTF:
        if NPad == 2:   
            canvas_controlplot = ROOT.TCanvas("TFFit_canvas", "canvas", 500,500);
            fPads1 = ROOT.TPad("TFFit_pad1", "", 0.0, 0.29, 1.00, 1.00);
            fPads2 = ROOT.TPad("TFFit_pad2", "", 0.0, 0.00, 1.00, 0.29);
            fPads1.SetBottomMargin(0.1);fPads1.SetLeftMargin( 0.20);fPads1.SetRightMargin( 0.03);
            fPads2.SetLeftMargin(  0.10 );fPads2.SetRightMargin(0.03);fPads2.SetBottomMargin(0.25);
            fPads1.Draw(); fPads2.Draw(); fPads1.cd();
        if NPad == 1 :   
            canvas_controlplot = ROOT.TCanvas("TFFit_canvas", "canvas", 600,500);
            canvas_controlplot.SetBottomMargin(0.1);canvas_controlplot.SetLeftMargin(0.1); canvas_controlplot.SetRightMargin(0.05);
        
        theLeg = ROOT.TLegend(0.3, 0.7, 0.3, 0.9, "", "NDC");theLeg.SetName("theLegend"); theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42);theLeg.SetTextSize(.05);
        theLeg.SetFillColor(0);theLeg.SetBorderSize(0);theLeg.SetLineColor(0);theLeg.SetLineWidth(0);theLeg.SetLineStyle(0);theLeg.SetTextFont(42);#theLeg.SetNColumns(2);
        theLeg.AddEntry(histogram,HistLabel,"P");
        theLeg.AddEntry(polfunc,function+" Fit","L");
        
        histogram.Draw("e")
#         histogram.GetYaxis().SetRangeUser(  0,1.3*ROOT.TMath.Max(0,F_GetMaxPlusErr(hint)) );
        histogram.GetYaxis().SetRangeUser(  0,maxY );
        histogram.GetYaxis().SetTitle(YLabel)
        histogram.GetYaxis().SetTitleOffset(0.45)
        histogram.GetYaxis().SetTitleSize(0.05)
        histogram.GetYaxis().SetLabelSize(0.05)
        histogram.GetXaxis().SetLabelSize(0.05)
        histogram.GetXaxis().SetTitleOffset(1.5)
        histogram.GetXaxis().SetTitleSize(0.05)

        histogram.Draw("e")
        
        if AddFitFunction :
            polfunc.SetLineWidth(1)
            polfunc.Draw("same")
#             hint.Draw("e3 same")
            ptstats.Draw("same")
            hintUpDraw.SetLineColor(4)
            hintDnDraw.SetLineColor(4)
            AddFitFunction.SetLineWidth(1)
            AddFitFunction.SetLineColor(3)
            AddFitFunction.Draw("same")
            hintUpDraw.Draw("c same")
            hintDnDraw.Draw("c same")
            hint.Draw("e3 same")
        elif AddHistogram:
            polfunc.SetLineWidth(1)
            polfunc.Draw("same")
#             hint.Draw("e3 same")
            ptstats.Draw("same")
            hintUpDraw.SetLineColor(4)
            hintDnDraw.SetLineColor(4)
            AddHistogram.SetLineWidth(3)
            AddHistogram.SetLineColor(9)
            AddHistogram.Draw("same")
            hintUpDraw.Draw("c same")
            hintDnDraw.Draw("c same")
            hint.Draw("e3 same")            
        else :
            polfunc.Draw("same")
            hint.Draw("e3 same")
            ptstats.Draw("same")
            hintUpDraw.Draw("c same")
            hintDnDraw.Draw("c same")
        theLeg.Draw()
        
        path = os.path.dirname(PlotName)
        if not os.path.isdir(path):
            os.makedirs(path)
        canvas_controlplot.SaveAs(PlotName+".png")
        canvas_controlplot.SaveAs(PlotName+".pdf")
        
    dict_ = {
        "r":r,
        "hintDnDraw":hintDnDraw,
        "hintUpDraw":hintUpDraw,
    }
    return histogram,hint,polfunc,chi,r.Ndf(),dict_

def F_GetMaxPlusErr(histogram):
    max_ = 0
    for i in range(histogram.GetNbinsX()):
        c = histogram.GetBinContent(i+1)
        e = histogram.GetBinError(i+1)
        if max_ < (c+e):
            max_ = c+e
    return max_

def F_GetHistograms(File,fix = "",Hist=None):
    Histograms = {}
    fn = ROOT.TFile.Open(File)
    for e in fn.GetListOfKeys() : 
        if Hist:
            if e.GetName() not in Hist: 
                continue
        h = fn.Get(e.GetName())
        h_Tmp = h.Clone(e.GetName()+fix)
        Histograms[e.GetName()] = h_Tmp
        h.SetDirectory(0)
        h_Tmp.SetDirectory(0)
    fn.Close()
    return Histograms

def F_GetFunctions(File,Hist=None):
    Functions = {}
    fn = ROOT.TFile.Open(File)
    for e in fn.GetListOfKeys() : 
        if Hist:
            if e.GetName() not in Hist: 
                continue
        h = fn.Get(e.GetName())
        Functions[e.GetName()] = h
    fn.Close()
    return Functions

def F_FValue(rss0,rss1,p1,p2):
    return (rss0-rss1)/p1/(rss1/p2)

def F_PValue(Ftest,p1,p2):
    print Ftest,p1,p2
    PValue = 1.-ROOT.TMath.FDistI(Ftest,p1,p2)
    return PValue

import os
def F_FTest(Histogram,NPol = 11,plotdir="TF/",HistLabel="test"):
    if not os.path.isdir(plotdir) :
        os.makedirs(plotdir)
    Pvalues = {}
    Chis = {} 
    Ndfs = {}
    for i in range(1,NPol):
        Func = "pol"+str(i)
        h_Ratio,_,Fit,chi,Ndf = F_FitHistogram(Histogram,Func,min=0,max=4000,PlotName=plotdir+Func,ColorFunc=2,HistColor=1,HistLabel=HistLabel,NPad=1,Name="Rb",DrawTF=True)
        Chis[i] = chi 
        Ndfs[i] = Ndf
        print i
        if Ndf == 1 :
            break
    for i in range(1,len(Chis)):
        print i
        Chi1 = Chis[i] ; Chi2 = Chis[i+1]
        Ndf1 = Ndfs[i] ; Ndf2 = Ndfs[i+1]
        FValue = F_FValue(Chi1,Chi2,1,Ndf2)
        Pvalues[i] = F_PValue(FValue,1,Ndf2)
    return Pvalues,Chis,Ndfs

def F_FTest2(Histogram,CPol,NPol = 11,plotdir="TF/",HistLabel="test"):
    if not os.path.isdir(plotdir) :
        os.makedirs(plotdir)
    Pvalues = {}
    Chis = {} 
    Ndfs = {}
    for i in range(1,NPol):
        Func = "pol"+str(i)
        h_Ratio,_,Fit,chi,Ndf = F_FitHistogram(Histogram,Func,min=0,max=4000,PlotName=plotdir+Func,ColorFunc=2,HistColor=1,HistLabel=HistLabel,NPad=1,Name="Rb",DrawTF=True)
        Chis[i] = chi 
        Ndfs[i] = Ndf
        print i
        if Ndf == 1 :
            break
    for i in range(1,len(Chis)):
        if i == CPol : 
            continue
        if i < CPol :
            Chi1 = Chis[i] ; Chi2 = Chis[CPol]
            Ndf1 = Ndfs[i] ; Ndf2 = Ndfs[CPol]
            FValue = F_FValue(Chi1,Chi2,CPol-i,Ndf2)
            Pvalues[i] = F_PValue(FValue,CPol-i,Ndf2)
        if i > CPol :
            Chi1 = Chis[CPol] ; Chi2 = Chis[i]
            Ndf1 = Ndfs[CPol] ; Ndf2 = Ndfs[i]
            FValue = F_FValue(Chi1,Chi2,i-CPol,Ndf2)
            Pvalues[i] = F_PValue(FValue,i-CPol,Ndf2)
    return Pvalues,Chis,Ndfs

def F_GetHistograms(File,fix=""):
    Histograms = {}
    fn = ROOT.TFile.Open(File)
    for e in fn.GetListOfKeys(): 
        h = fn.Get(e.GetName())
        h.SetName(e.GetName()+fix)
        Histograms[h.GetName()] = h
        h.SetDirectory(0)
    fn.Close()
    return Histograms

def F_AddLegend(Position,Pairs,TextSize):
    theLeg = TLegend(Position[0], Position[1], Position[2], Position[3],"","NDC")
    theLeg.SetName("theLegend")
    theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42);theLeg.SetTextSize(TextSize)
    theLeg.SetFillColor(0);theLeg.SetBorderSize(0);theLeg.SetLineColor(0);theLeg.SetLineWidth(0);theLeg.SetLineStyle(0);theLeg.SetTextFont(42)
    for pair in Pairs :
        theLeg.AddEntry(pair[0],pair[1],pair[2])
    return theLeg

def F_Store_Histogram(Histograms,OutFile):
    OutFile = os.path.normpath(OutFile)     
    OutputDir = os.path.dirname(OutFile)
    if not os.path.isdir(OutputDir) :
        os.makedirs(OutputDir)     
    outf = ROOT.TFile( OutFile , "recreate")
    for ih in Histograms : 
        Histograms[ih].Write()
    outf.Close()
    print "create %s"%(OutFile)

def F_GetHistograms(File,fix = "",Hist=None):
    Histograms = {}
    fn = ROOT.TFile.Open(File)
    for e in fn.GetListOfKeys() : 
        if Hist:
            if e.GetName() not in Hist: 
                continue
        h = fn.Get(e.GetName())
        h_Tmp = h.Clone(e.GetName()+fix)
        Histograms[e.GetName()] = h_Tmp
        h.SetDirectory(0)
        h_Tmp.SetDirectory(0)
    fn.Close()
    return Histograms

def F_GetFunctions(File,Hist=None):
    Functions = {}
    fn = ROOT.TFile.Open(File)
    for e in fn.GetListOfKeys() : 
        if Hist:
            if e.GetName() not in Hist: 
                continue
        h = fn.Get(e.GetName())
        Functions[e.GetName()] = h
    fn.Close()
    return Functions

def F_SetRangeUser(Histogram,Histograms,min_ = None,max_ = None, UpScale = 1.3, DnScale = 0.7):
    Yaxis = Histogram.GetYaxis()
    if min_ == None :
        min_ = min([Histograms[Name].GetMinimum() for Name in Histograms])*DnScale
    if max_ == None :
        max_ = max([Histograms[Name].GetMaximum() for Name in Histograms])*UpScale
    print min_,max_
    Yaxis.SetRangeUser(min_, max_)
    
def F_AddLegend(Position,Pairs,TextSize):
    theLeg = ROOT.TLegend(Position[0], Position[1], Position[2], Position[3],"","NDC")
    theLeg.SetName("theLegend")
    theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42);theLeg.SetTextSize(TextSize)
    theLeg.SetFillColor(0);theLeg.SetBorderSize(0);theLeg.SetLineColor(0);theLeg.SetLineWidth(0);theLeg.SetLineStyle(0);theLeg.SetTextFont(42)
    for pair in Pairs :
        theLeg.AddEntry(pair[0],pair[1],pair[2])
    return theLeg

def F_GetDataPull(h_data,h_diff):
    h_data.SetBinErrorOption(ROOT.TH1D.kPoisson);
    h_pull = h_data.Clone("pull")
    for i in range(1,h_data.GetNbinsX()+1):
        v = h_diff.GetBinContent(i)/h_data.GetBinError(i)
#         print v
        h_pull.SetBinContent(i,v)
        h_pull.SetBinError(i,1)
    return h_pull

def construct_ratio_histogram(Name,h_Nume,h_Deno,color=1,maxY=100000000):
    h_Ratio=h_Nume.Clone("h_Ratio"+Name)
    h_Ratio.SetLineColor(color)
    h_Ratio.SetLineWidth(3)
    h_Ratio.SetLineStyle(7)
    h_Ratio.SetMarkerColor(color)
    h_Ratio.SetMarkerStyle(8);
    h_Ratio.SetMarkerSize(0);
    h_Ratio.GetYaxis().SetRangeUser( 0 , maxY )
    h_Ratio.GetYaxis().SetNdivisions(504,0)
    h_Ratio.GetYaxis().SetTitle("Ratios   ")
    h_Ratio.GetYaxis().SetTitleOffset(0.35)
    h_Ratio.GetYaxis().SetTitleSize(0.13)
    h_Ratio.GetYaxis().SetTitleSize(0.13)
    h_Ratio.GetYaxis().SetLabelSize(0.11)
    h_Ratio.GetXaxis().SetLabelSize(0.1)
    h_Ratio.GetXaxis().SetTitleOffset(1.0)
    h_Ratio.GetXaxis().SetTitleSize(0.1)
    h_Ratio.SetStats(0)
    for i in range(1,h_Ratio.GetNbinsX()+1,1):
        D  = h_Nume.GetBinContent(i);    eD = h_Nume.GetBinError(i);
        B  = h_Deno.GetBinContent(i); eB = h_Deno.GetBinError(i);
        # print i,D,B
        if D==0: eD=0.92;
        if B<0.1 and eB>=B : eB=0.92; Err= 0.;
        if B!=0.        :Err=TMath.Sqrt( (eD*eD)/(B*B)  +(D*D*eB*eB)/(B*B*B*B)     ); h_Ratio.SetBinContent(i, D/B   );  h_Ratio.SetBinError(i, Err);
        if B==0.        :Err=TMath.Sqrt( (eD*eD)/(eB*eB)+(D*D*eB*eB)/(eB*eB*eB*eB) ); h_Ratio.SetBinContent(i, D/0.92);  h_Ratio.SetBinError(i, Err);
        if D==0 and B==0:                                                             h_Ratio.SetBinContent(i, -1);      h_Ratio.SetBinError(i, 0  );
    return h_Ratio

def F_TF_Mutiply_TH(TF,histogram):
    Name = histogram.GetName()+"_TF1"
    histout = histogram.Clone(Name)
    for i in range(1,histout.GetNbinsX()+1):
        v = histogram.GetBinContent(i)*TF.Eval(histout.GetBinCenter(i))
        histout.SetBinContent(i,v)
    return histout

def F_TFParamUn(histogram,func,nparam,r):
    Name = histogram.GetName()+"_TF%s"%(str(nparam))
    NameUp = Name+"Up"
    NameDown = Name+"Down"
    h_Up = histogram.Clone(NameUp)
    h_Down = histogram.Clone(NameDown)

    funcUp = func.Clone("TFUp")
    funcDown = func.Clone("TFDown")

    value = func.GetParameter(nparam)
    error = r.ParError(nparam)*0.01
    print value,error,error/value

    funcUp.SetParameter(nparam,value+error)
    funcDown.SetParameter(nparam,value-error)

    h_Up = F_TF_Mutiply_TH(funcUp,h_Up)
    h_Down = F_TF_Mutiply_TH(funcDown,h_Down)

    return h_Up,h_Down,funcUp,funcDown

def F_GetParams(r,n,nI=3):
    X = []
    for i in range(n):
        v = r.Parameter(i)
        e = r.ParError(i)
        X.append(
            (v,v-nI*e,v+nI*e)
        )
    return X

def F_RooRealVars(params,VarName="A"):
    Vars = []
    for index,x in enumerate(params):
        Name = VarName+str(index)
        tmp = ROOT.RooRealVar(Name, Name, x[0], x[1], x[2])
        Vars.append(tmp)
    return Vars

def F_RooPolynomial(x,Vars,Name="polyA"):
    poly = ROOT.RooPolynomial(Name, Name, x, RooArgList(*Vars), 0)
    return poly