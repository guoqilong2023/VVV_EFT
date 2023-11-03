import re
processes = ["WWW_Dim6","WWZ_Dim6","WZZ_Dim6","ZZZ_Dim6",]
Hists = {}
searchHists = re.compile(r"h_(.*)_Dim6_cW_0p3_MJJlv(.*)")
HistName = "h_%s_cW_0p3_MJJlv"

Files = glob.glob("/eos/user/q/qiguo/WVV_EFT/plot/basic/V2/Signal_Systs/ROOTFiles/*_cW_*.root")
for File in Files:
    TmpRF = ROOT.TFile(File)
    for e in TmpRF.GetListOfKeys() : 
        h = TmpRF.Get(e.GetName())
        if e.GetName() not in Hists:
            Hists[e.GetName()] = h
            h.SetDirectory(0)
    TmpRF.Close()

SelectedHists = {}
for h in Hists:
    if searchHists.search(h):
        SelectedHists[h] = Hists[h]

UNs = []
searchUN = re.compile(r"h_WWW_Dim6_cW_0p3_MJJlv_(.*)")
for h in SelectedHists:
    if searchUN.search(h):
        UNs.append(searchUN.search(h).groups()[0])

VVVCenter = SelectedHists[HistName%(processes[0])].Clone(HistName%("VVV"))
for p in processes[1:]:
    VVVCenter.Add(SelectedHists[HistName%(p)])

VVVSysts = {}
for UN in UNs:
    Name = "VVV"+UN+"Up"
    VVV = SelectedHists[HistName%(processes[0])+"_"+UN+"Up"].Clone(HistName%("VVV")+UN+"Up")
    for p in processes[1:]:
        VVV.Add(SelectedHists[HistName%(p)+"_"+UN+"Up"])
    VVVSysts[Name] = VVV

    Name = "VVV"+UN+"Dn"
    VVV = SelectedHists[HistName%(processes[0])+"_"+UN+"Dn"].Clone(HistName%("VVV")+UN+"Dn")
    for p in processes[1:]:
        VVV.Add(SelectedHists[HistName%(p)+"_"+UN+"Dn"])
    VVVSysts[Name] = VVV
    
nbin = VVVCenter.GetNbinsX()

for i in range(1,nbin+1):
    min_ = histogram.GetBinLowEdge(i)
    max_ = histogram.GetBinLowEdge(i+1)
    value = VVVCenter.GetBinContent(i)
    error = VVVCenter.GetBinError(i)
    systs_square = 0
    for UN in UNs:
        Name = "VVV"+UN+"Up"
        up = VVVSysts[Name].GetBinContent(i)
        Name = "VVV"+UN+"Dn"
        dn = VVVSysts[Name].GetBinContent(i)
        systs_square += (up-dn)**2
    systs = systs_square**0.5
    print "%s-%s"%(min_,max_),value,error,systs