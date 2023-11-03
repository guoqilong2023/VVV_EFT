import ROOT
import os,glob,re

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

def F_AddLegend(Position,Pairs,TextSize):
    theLeg = ROOT.TLegend(Position[0], Position[1], Position[2], Position[3],"","NDC")
    theLeg.SetName("theLegend")
    theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42);theLeg.SetTextSize(TextSize)
    theLeg.SetFillColor(0);theLeg.SetBorderSize(0);theLeg.SetLineColor(0);theLeg.SetLineWidth(0);theLeg.SetLineStyle(0);theLeg.SetTextFont(42)
    for pair in Pairs :
        theLeg.AddEntry(pair[0],pair[1],pair[2])
    return theLeg

def F_CompareShape(Histograms,Labels,Color1 = [1,2,3,4,6,7],Histograms2=[],Labels2=[],Color2 = [],Yrange = (0,2), plot = "InputShape.png", REGION = None, Ratio_Draw_OP = None):
    
    Hists = []
    for H in Histograms:
        newH = H.Clone(H.GetName()+"_Comptmp")
#         newH.Scale(1/newH.Integral())
        Hists.append(newH)
    
    for index,H in enumerate(Hists):
        H.SetLineWidth(3)
        H.SetLineColor(Color1[index])
        H.SetFillColor(0)
        H.SetMarkerSize(0)
        H.SetStats(0)    
        
    Hists2 = []
    for H in Histograms2:
        newH = H.Clone(H.GetName()+"_Comptmp")
        Hists2.append(newH)
    
    for index,H in enumerate(Hists2):
        H.SetLineWidth(3)
        H.SetLineColor(Color2[index])
        H.SetFillColor(0)
        H.SetMarkerColor(Color2[index])
        H.SetMarkerSize(2)
        if Ratio_Draw_OP == "e0":
            H.SetMarkerStyle(8)
        else:
            H.SetMarkerStyle(1)
        H.SetStats(0)    

    Legends = [(H,Labels[index],"L") for index,H in enumerate(Hists)]
    Legend1 = F_AddLegend((0.55,0.7,0.9,0.96),Legends,0.1)

    max_ = max([H.GetBinContent(H.GetMaximumBin()) for H in Hists])
    canvas = ROOT.TCanvas("canvas","",700,700)
    fPads1 = ROOT.TPad("TFFit_pad1", "", 0.0, 0.6, 1.00, 1.00);
    fPads2 = ROOT.TPad("TFFit_pad2", "", 0.0, 0.00, 1.00, 0.6);
    fPads1.SetTopMargin(0.01);fPads1.SetBottomMargin(0.03);fPads1.SetLeftMargin( 0.10);fPads1.SetRightMargin( 0.03);
    fPads2.SetTopMargin(0.03);fPads2.SetLeftMargin(  0.10 );fPads2.SetRightMargin(0.03);fPads2.SetBottomMargin(0.15);
    fPads1.Draw(); fPads2.Draw(); 
    
    fPads1.cd()
    for index,H in enumerate(Hists):
        if index == 0:
            H.Draw("hist e")
            H.GetYaxis().SetRangeUser( 0 , max_*1.3 )
            H.GetYaxis().SetTitleSize( 0.08 )
            H.GetYaxis().SetTitleOffset(0.6)
            H.GetYaxis().SetLabelSize(0.08)
            H.GetXaxis().SetLabelSize(0.)
        else:
            H.Draw("same hist e")
    Legend1.Draw()
    if REGION:
        RegionTxt       = ROOT.TLatex(0.15,0.88,"%s"%(REGION)            );RegionTxt.SetNDC();RegionTxt.SetTextSize(0.072);  RegionTxt.SetTextFont(42);    RegionTxt.SetLineWidth(2);     
        RegionTxt.Draw("same")
    
    fPads2.cd()
    Legends2 = [(H,Labels2[index],"L") for index,H in enumerate(Hists2)]
    Legend2 = F_AddLegend((0.15,0.75,0.45,0.9),Legends2,0.06)

    for index,H in enumerate(Hists2):
        Draw_op = "hist e"
        if Ratio_Draw_OP:
            Draw_op =   Ratio_Draw_OP
        if index == 0:
            H.Draw(Draw_op)
            H.GetYaxis().SetRangeUser( Yrange[0] , Yrange[1] )
            H.GetYaxis().SetTitle("Ratio")
            H.GetYaxis().SetTitleSize( 0.05 )
            H.GetYaxis().SetTitleOffset(0.9)
            H.GetYaxis().SetLabelSize(0.06)
            H.GetXaxis().SetTitle("M_{JJJ}*")
            H.GetXaxis().SetLabelSize(0.06)
            H.GetXaxis().SetTitleOffset(0.95)
        else:
            H.Draw("same %s"%(Draw_op))
    axis1=ROOT.TLine( H.GetBinLowEdge(1),1,H.GetBinLowEdge(H.GetNbinsX()+1),1); axis1.SetLineColor(1); axis1.SetLineWidth(1); axis1.Draw("same");
    Legend2.Draw()
    
    canvas.SaveAs(plot)
    canvas.SaveAs(plot.replace(".png",".pdf"))
    return None


Version = "V4"
Path = "/eos/user/q/qiguo/www/VVV/Limit/%s/"%(Version)
if not os.path.exists(Path):
    os.makedirs(Path)
SystFile = "/eos/user/q/qiguo/www/VVV/Limit/%s/Syst.p"%(Version)
InputRFile = "/eos/user/q/qiguo/www/VVV/Limit/%s/LimitInput.root"%(Version)
BKGInputRFile = "/eos/user/q/qiguo/www/VVV/Limit/%s/BKGLimitInput.root"%(Version)
OuputDir = "/eos/user/q/qiguo/www/VVV/Limit/%s/"%(Version)
plotPath = "/eos/user/q/qiguo/www/VVV/Limit/%s/SignalUnPlots/"%(Version)

Hists = {}
BKGRF = ROOT.TFile(BKGInputRFile)
for e in BKGRF.GetListOfKeys() : 
    h = BKGRF.Get(e.GetName())
    Hists[e.GetName()] = h
    h.SetDirectory(0)
BKGRF.Close()

# Signal_Systs_Path = "/eos/user/q/qiguo/WVV_EFT/plot/basic/V2/Signal_Systs/ROOTFiles/*.root"
Signal_Systs_Path = "/eos/user/q/qiguo/WVV_EFT/plot/basic/V2/Signal_Systs/ROOTFiles/*_*.root"
Files = glob.glob(Signal_Systs_Path)
for index,File in enumerate(Files):
    TmpRF = ROOT.TFile(File)
    for e in TmpRF.GetListOfKeys() : 
        h = TmpRF.Get(e.GetName())
        if e.GetName() not in Hists:
            Hists[e.GetName()] = h
            h.SetDirectory(0)
    TmpRF.Close()
    if index%10 == 0:
        print "done",index

# F_Store_Histogram(Hists,"/eos/user/q/qiguo/www/VVV/Limit/V4/LimitInput_New3.root")
F_Store_Histogram(Hists,"/eos/user/q/qiguo/www/VVV/Limit/V4/LimitInput_test.root")

UNs = []
strips = ["Up", "Down", "Dn", "up", "dn"]
searchUN = re.compile(r"h_WWW_Dim6_cW_0p3_MJJlv_(.*)")
for h in Hists:
    if searchUN.search(h):
        UN = searchUN.search(h).groups()[0]
        for s in strips:
            UN = UN.strip(s)
        UNs.append(UN)
UNs = list(set(UNs))
    
processes = [
    "WWW_Dim6",
    "WWZ_Dim6",
    "WZZ_Dim6",
    "ZZZ_Dim6",
]
All_OPs = {}
OPs = [
    "cW",
    "cHbox","cHDD","cHW","cHB","cHWB","cHl3","cHq1","cHq3","cHu","cHd","cll1"
] 
termsin = ["sm","quad","sm_lin_quad"]
terms = ["sm","quad","lin"]

for OP in OPs:
    for UN in UNs:
        for process in processes:
            # get names
            HistNames = list(Hists.keys())
            HistName = "h_sm_{process}_{OP}_{UN}".format(
                process = process,
                OP = OP,
                UN = UN,
            )
            centername = "h_sm_{process}_{OP}".format(
                process = process,
                OP = OP,
            )
            for name in HistNames:
                if HistName in name:
                    if "up" in name.lower():
                        upname = name
                    if ("dn" in name.lower()) or ("down" in name.lower()):
                        dnname = name

            # get linear term
            h_c_q = Hists.get(centername.replace("sm","quad"))
            h_u_q = Hists.get(upname.replace("sm","quad"))
            h_d_q = Hists.get(dnname.replace("sm","quad"))

            h_c_s = Hists.get(centername.replace("sm","sm"))
            h_u_s = Hists.get(upname.replace("sm","sm"))
            h_d_s = Hists.get(dnname.replace("sm","sm"))

            h_c_l = Hists.get(centername.replace("sm","sm_lin_quad")).Clone(centername.replace("sm","lin"))
            h_u_l = Hists.get(upname.replace("sm","sm_lin_quad")).Clone(upname.replace("sm","lin"))
            h_d_l = Hists.get(dnname.replace("sm","sm_lin_quad")).Clone(dnname.replace("sm","lin"))

            h_c_l.Add(h_c_q,-1);h_c_l.Add(h_c_s,-1)
            h_u_l.Add(h_u_q,-1);h_u_l.Add(h_u_s,-1)
            h_d_l.Add(h_d_q,-1);h_d_l.Add(h_d_s,-1)

            # make plot
            for term in ["s","l","q"]:
                h_S_Center = eval("h_c_%s"%(term)).Clone("h_S_Center")
                h_Ratio_S_Up = eval("h_u_%s"%(term)).Clone("h_Ratio_S_Up")
                h_Ratio_S_Dn = eval("h_d_%s"%(term)).Clone("h_Ratio_S_Dn")
                h_Ratio_S_Up.Divide(h_S_Center)
                h_Ratio_S_Dn.Divide(h_S_Center)

                PlotName = "%s_%s_%s_%s.png"%(process, OP, term, UN)
                plotPathsub = os.path.join(plotPath,OP,UN) 
                if not os.path.isdir(plotPathsub):
                    os.makedirs(plotPathsub)
                os.system("cp /afs/cern.ch/work/q/qiguo/public/gKK/plot/v5/index.php %s"%(plotPathsub))
                Plot = os.path.join(plotPathsub,PlotName)
                print Plot

                if eval("h_c_%s"%(term)).Integral() != 0:
                    RatioUp = str(round(eval("h_u_%s"%(term)).Integral()/eval("h_c_%s"%(term)).Integral(),3))
                    RatioDn = str(round(eval("h_d_%s"%(term)).Integral()/eval("h_c_%s"%(term)).Integral(),3))
                else:
                    RatioUp = ""
                    RatioDn = ""

                F_CompareShape(
                [eval("h_c_%s"%(term)),eval("h_u_%s"%(term)),eval("h_d_%s"%(term))],["Center","%s (up) "%(UN)+RatioUp,"%s (down) "%(UN)+RatioDn],[1,2,4,],
                [h_Ratio_S_Up,h_Ratio_S_Dn],["%s (up)"%(UN),"%s (down)"%(UN),],[2,4],
                Yrange = (0,2.),
                plot = Plot,
                REGION = "%s, %s, %s"%(process.replace("_"," "), OP, term),
                )



# import re
# def F_ExtractHistogram(Hists,TermName):
#     # get un name
#     # get center,up,dn histogram
#     # return information 
#     return h_c,h_up,h_dn,unname,process

# def F_PlotUncertainty(h_c,h_up,h_dn,unname,termname):

    
# def F_Uncertainty_Table(h_c,h_up,h_dn,unname,termname):


# {
#     "WWW_Dim6":{
#         "cW":{
#             ("jer","(sm,bin1)"): 0.01,
#         }
#     }
# }
