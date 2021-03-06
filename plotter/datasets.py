import ROOT
import copy
import auxiliary as aux
import os.path
import style
from rwthColors import rwth


path = "../histogramProducer/"

class Dataset:
    names = []
    files = []
    xsecs = []
    lname = []

    color = None
    label = ""

    def __radd__(self, dset):
        if not dset: return self

        out = copy.deepcopy(self)

        out.names.extend(dset.names)
        out.files.extend(dset.files)
        out.xsecs.extend(dset.xsecs)
        out.label += " + "+dset.label
        out.lname.extend(dset.lname)

        if not self.color:
            out.color = dset.color
        return out

    def __init__(self, n, xsec=-1, col=ROOT.kBlack, fullname="",):
        fname = path + n + "_combiHists.root"
        if xsec == -1: xsec = aux.getXsecFromName(n)
        self.names = [ n ]
        self.files = [ fname ]
        self.xsecs = [ xsec ]
        self.color = col
        self.label = n
        self.lname = [ fullname ]

    def mcm(self):
        for fullname in self.lname:
            print "https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name="+fullname


    def __str__(self):
        return "Dataset: " + self.label + "\ncolor: " + str(self.color) + \
            "\nnames: "+", ".join(self.names) + \
            "\nfiles: "+", ".join(self.files) + \
            "\nxsecs: "+", ".join(str(i) for i in self.xsecs)

    def getHist(self, name):
        h0 = None
        for i in range(len(self.files)):
            h = aux.getFromFile(self.files[i], name)
            if isinstance(h, ROOT.TH1):
                if self.xsecs[i]:
                    if style.additionalPoissonUncertainty:
                        aux.addPoissonUncertainty(h)
                    if type(self.xsecs[i]) == str:
                        mass =int(name.split("/")[0].split("_")[0])
                        xsec = aux.getXsecInfoSMS(mass, self.xsecs[i])[0]
                    else: xsec = self.xsecs[i]
                    h.Scale(aux.intLumi * xsec)
                h.SetLineColor(self.color)
                h.SetMarkerColor(self.color)
            if h0: h0.Add(h)
            else: h0 = h
        return h0

    def getHistFromTree(self, variable, weight, nBins, tname="simpleTree"):
        h0 = None
        for i in range(len(self.files)):
            tree = ROOT.TChain(tname)
            tree.AddFile(self.files[i])
            h = aux.createHistoFromTree(tree, variable, weight, nBins)
            if isinstance(h, ROOT.TH1):
                if self.xsecs[i]:
                    if style.additionalPoissonUncertainty:
                        aux.addPoissonUncertainty(h)
                    if type(self.xsecs[i]) == str:
                        mass = int(name.split("/")[0].split("_")[0])
                        self.xsecs[i] = aux.getXsecInfoSMS(mass, self.xsecs[i])[0]
                    h.Scale(aux.intLumi * self.xsecs[i])
                h.SetLineColor(self.color)
                h.SetMarkerColor(self.color)
            if h0: h0.Add(h)
            else: h0 = h
        #if style.divideByBinWidth: h0.Scale(1., "width")
        return h0

data = Dataset("SinglePhoton_Run2016B-03Feb2017_ver2-v2", 0, ROOT.kBlack) \
    + Dataset("SinglePhoton_Run2016C-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("SinglePhoton_Run2016D-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("SinglePhoton_Run2016E-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("SinglePhoton_Run2016F-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("SinglePhoton_Run2016G-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("SinglePhoton_Run2016H-03Feb2017_ver2-v1", 0, ROOT.kBlack) \
    + Dataset("SinglePhoton_Run2016H-03Feb2017_ver3-v1", 0, ROOT.kBlack)
data.label = "Data"

dataHt = Dataset("JetHT_Run2016B-03Feb2017_ver2-v2", 0, ROOT.kBlack) \
    + Dataset("JetHT_Run2016C-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("JetHT_Run2016D-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("JetHT_Run2016E-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("JetHT_Run2016F-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("JetHT_Run2016G-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("JetHT_Run2016H-03Feb2017_ver2-v1", 0, ROOT.kBlack) \
    + Dataset("JetHT_Run2016H-03Feb2017_ver3-v1", 0, ROOT.kBlack)
dataHt.label = "Data"


###############################################################################
# Simulation
###############################################################################

###############################################################################
# GJet
gjets40dr = Dataset("GJets_DR-0p4_HT-40To100", 17420.0, ROOT.kCyan+4, "GJets_DR-0p4_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
gjets100dr = Dataset("GJets_DR-0p4_HT-100To200", 5383.0, ROOT.kCyan-1, "GJets_DR-0p4_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
gjets200dr = Dataset("GJets_DR-0p4_HT-200To400", 1176.0, ROOT.kCyan+2, "GJets_DR-0p4_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
gjets400dr = Dataset("GJets_DR-0p4_HT-400To600", 132.1, ROOT.kCyan+3, "GJets_DR-0p4_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
gjets600dr = Dataset("GJets_DR-0p4_HT-600ToInf", 44.32, ROOT.kCyan+2, "GJets_DR-0p4_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
gjets_dr = gjets40dr + gjets100dr + gjets200dr + gjets400dr + gjets600dr
gjets_dr.label = "#gamma+jet"

###############################################################################
# Multijet
qcd50 = Dataset("QCD_HT50to100", 246400000, ROOT.kBlue+1, "QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
qcd100 = Dataset("QCD_HT100to200", 27990000, ROOT.kBlue+1, "QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
qcd200 = Dataset("QCD_HT200to300", 1712000, ROOT.kBlue+2, "QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
qcd300 = Dataset("QCD_HT300to500", 347700, ROOT.kBlue+3, "QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
qcd500 = Dataset("QCD_HT500to700", 32100, ROOT.kBlue+4, "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
qcd700 = Dataset("QCD_HT700to1000", 6831, ROOT.kBlue+3, "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
qcd1000 = Dataset("QCD_HT1000to1500", 1207, ROOT.kBlue+2, "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
qcd1500 = Dataset("QCD_HT1500to2000", 119.9, ROOT.kBlue+1, "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
qcd2000 = Dataset("QCD_HT2000toInf", 25.24,  ROOT.kBlue, "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
qcd = qcd500 + qcd700 + qcd1000 + qcd1500 + qcd2000
qcd.label = "Multijet"

gqcd = gjets_dr + qcd
gqcd.label = "(#gamma)+jet"

###############################################################################
# WJet
wjets100 = Dataset("WJetsToLNu_HT-100To200", 1345., ROOT.kRed-6, "WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
wjets200 = Dataset("WJetsToLNu_HT-200To400", 359.7, ROOT.kRed-5, "WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
wjets400 = Dataset("WJetsToLNu_HT-400To600", 48.91, ROOT.kRed-4, "WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
wjets600 = Dataset("WJetsToLNu_HT-600To800", 12.05, ROOT.kRed-3, "WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
wjets800 = Dataset("WJetsToLNu_HT-800To1200", 5.501, ROOT.kRed-2, "WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
wjets1200 = Dataset("WJetsToLNu_HT-1200To2500", 1.329, ROOT.kRed-1, "WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
wjets2500 = Dataset("WJetsToLNu_HT-2500ToInf", 0.03216, ROOT.kRed, "WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
wjetsSamples = [ wjets200, wjets400, wjets600, wjets800, wjets1200, wjets2500 ]
for ds in wjetsSamples: ds.xsecs = [ ds.xsecs[0] * 1.21 ] # k-factor
wjets = sum(wjetsSamples)
wjets.label = "W#rightarrowl#nu"

###############################################################################
# ZJet
znunu100 = Dataset("ZJetsToNuNu_HT-100To200", 280.47, ROOT.kMagenta-2 , "ZJetsToNuNu_HT-100To200_13TeV-madgraph")
znunu200 = Dataset("ZJetsToNuNu_HT-200To400", 77.67, ROOT.kMagenta-1 , "ZJetsToNuNu_HT-200To400_13TeV-madgraph")
znunu400 = Dataset("ZJetsToNuNu_HT-400To600", 10.73, ROOT.kMagenta+4 , "ZJetsToNuNu_HT-400To600_13TeV-madgraph")
znunu600 = Dataset("ZJetsToNuNu_HT-600To800", 2.559, ROOT.kMagenta+3 , "ZJetsToNuNu_HT-600To800_13TeV-madgraph")
znunu800 = Dataset("ZJetsToNuNu_HT-800To1200", 1.1796, ROOT.kMagenta+2 , "ZJetsToNuNu_HT-800To1200_13TeV-madgraph")
znunu1200 = Dataset("ZJetsToNuNu_HT-1200To2500", 0.28833, ROOT.kMagenta+1 , "ZJetsToNuNu_HT-1200To2500_13TeV-madgraph")
znunu2500 = Dataset("ZJetsToNuNu_HT-2500ToInf", 0.006945, ROOT.kMagenta+0 , "ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph")
znunu600Inf = Dataset("ZJetsToNuNu_HT-600ToInf", 4.116, ROOT.kMagenta+0 , "ZJetsToNuNu_HT-600ToInf_13TeV-madgraph")
znunuSamples = znunu100, znunu200, znunu400, znunu600, znunu800, znunu1200, znunu2500
znunuSamples = znunu400, znunu600, znunu800, znunu1200, znunu2500
for ds in znunuSamples: ds.xsecs = [ ds.xsecs[0] * 1.23 ]
znunu = sum(znunuSamples)
znunu.label = "Z#rightarrow#nu#nu"

###############################################################################
# TTbar
# https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns#TTbar
ttjets_nlo = Dataset("TTJets-amcatnloFXFX", 6.675e+02,  ROOT.kRed+2, "TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8")
ttjets_nlo.xsecs[0] = 831.76 # NLO; Uncertaniy: +19.77 -29.20 +35.06 -35.06
ttjets = Dataset("TTJets-madgraphMLM", 5.098e+02,  ROOT.kRed+2, "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
ttjets.xsecs[0] = 831.76 #+19.77 -29.20 +35.06 -35.06
ttjets.label = "t#bar{t}"

ttjets0 = Dataset("TTJets_HT-0to600", 5.098e+02, ROOT.kRed+2, "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
ttjets600 = Dataset("TTJets_HT-600to800_ext", 1.61, ROOT.kRed+2, "TTJets_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
ttjets800 = Dataset("TTJets_HT-800to1200_ext", 0.663, ROOT.kRed+2, "TTJets_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
ttjets1200 = Dataset("TTJets_HT-1200to2500_ext", 0.12, ROOT.kRed+2, "TTJets_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
ttjets2500 = Dataset("TTJets_HT-2500toInf_ext", 0.00143, ROOT.kRed+2, "TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
ttjets_ht_samples = ttjets0, ttjets600, ttjets800, ttjets1200, ttjets2500
for ds in ttjets_ht_samples: ds.xsecs = [ ds.xsecs[0] * 831.76/5.098e+02 ] # k-factor
ttjets_ht = sum(ttjets_ht_samples)
ttjets_ht.label = "t#bar{t}"

###############################################################################
# WG
# k-factor from http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2016/078 EXO-16-014
wg40 = Dataset("WGJets_MonoPhoton_PtG-40to130", 12.7, ROOT.kRed, "WGJets_MonoPhoton_PtG-40to130_TuneCUETP8M1_13TeV-madgraph")
#wg40.xsecs[0] *= 1.125/0.6565 # NLO, estimated for pt>130
wg130 = Dataset("WGJets_MonoPhoton_PtG-130", 0.6565, ROOT.kRed-1, "WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph")
#wg130.xsecs[0] = 1.125 # NLO
wg40.xsecs[0] *= 1.34 # from monophotons
wg130.xsecs[0] *= 1.34

wg = wg40 + wg130
wg.label = "#gammaW"
wg_nlo = Dataset("WGToLNuG-amcatnloFXFX_ext", 512.1, ROOT.kRed-2, "WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8")
wg_nlo.label = "#gammaW (NLO)"
wg130_nlo = Dataset("WGToLNuG_PtG-130-amcatnloFXFX", 1.125, ROOT.kRed-3, "WGToLNuG_PtG-130_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8")
wg130_nlo.label = "#gammaW (NLO), p^{#gamma}_{T}>130"

###############################################################################
# ZG
# k-factor from http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2016/078 EXO-16-014
zg40 = Dataset("ZNuNuGJets_MonoPhoton_PtG-40to130", 2.789*1.39, ROOT.kMagenta+1, "ZNuNuGJets_MonoPhoton_PtG-40to130_TuneCUETP8M1_13TeV-madgraph")
zg130 = Dataset("ZNuNuGJets_MonoPhoton_PtG-130", 0.1832*1.39, ROOT.kMagenta+2, "ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph")
zg = zg40 + zg130
zg.label = "#gammaZ(#nu#nu)"

zg_nlo = Dataset("ZGTo2NuG", 27.99, ROOT.kMagenta-3, "ZGTo2NuG_PtG-130_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8")
zg_nlo.label = "#gammaZ (NLO)"
zg130_nlo = Dataset("ZGTo2NuG_PtG-130", 0.2762, ROOT.kMagenta-1, "ZGTo2NuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8")
zg130_nlo.label = "#gammaZ (NLO), p^{#gamma}_{T}>130"

###############################################################################
# TTG
ttg = Dataset("TTGJets", 3.697, ROOT.kOrange, "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8")
ttg.label = "#gammat#bar{t}"

tg = Dataset("TGJets_amcatnlo_madspin", 2.967, ROOT.kOrange+2, "TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8")
tg.label = "#gammat"

###############################################################################
# Signal samples
t5wg_1600_100 = Dataset("SMS-T5Wg_1600_100", 0.00810078, ROOT.kRed, "")
t5wg_1600_100.label = "T5Wg 1600 100"
t5wg_1600_1500 = Dataset("SMS-T5Wg_1600_1500", 0.00810078, ROOT.kRed+4, "")
t5wg_1600_1500.label = "T5Wg 1600 1500"
t5wg_2000_100 = Dataset("SMS-T5Wg_2000_100", 0.000981077, ROOT.kRed+4, "")
t5wg_2000_100.label = "T5Wg 2000 100"
t5wg_1750_1700 = Dataset("SMS-T5Wg_1750_1700", 0.00359842, ROOT.kRed+4, "")
t5wg_1750_1700 = Dataset("SMS-T5Wg_1750_1700", 0.00359842, ROOT.kRed+4, "")
t5wg_1600_800 = Dataset("SMS-T5Wg_1600_800", 0.00810078, ROOT.kRed+4, "")
t5wg_1000_100 = Dataset("SMS-T5Wg_1000_100", 0.325388, ROOT.kRed+4, "")
t5wg_1000_100.label = "T5Wg 1000 100"

t6gg_1750_1650 = Dataset("SMS-T6gg_1750_1650", 0.000646271, ROOT.kRed+4, "")
t6gg_1750_1650.label = "T6gg 1750 1650"
t6gg_1300_600 = Dataset("SMS-T6gg_1300_600", 0.0086557, ROOT.kRed+4, "")
t6gg_1100_600 = Dataset("SMS-T6gg_1100_600", 0.0313372, ROOT.kRed+4, "")

t5wg = Dataset("SMS-T5Wg", "data/xSec_SMS_Gluino_13TeV.pkl")
t6wg = Dataset("SMS-T6Wg", "data/xSec_SMS_Squark_13TeV.pkl")
t5wg_ext = Dataset("SMS-T5Wg_mGo2150To2500", "data/xSec_SMS_Gluino_13TeV.pkl")
t6wg_ext = Dataset("SMS-T6Wg_mSq1850To2150", "data/xSec_SMS_Squark_13TeV.pkl")
tching = Dataset("SMS-TChiNG_BF50N50G", "data/xSec_SMS_TChiNG_13TeV_interpolated_5.pkl")
#ggm1 = Dataset("GGM_GravitinoLSP_M1-200to1500_M2-200to1500")

tchiwg_700 = Dataset("SMS-TChiWG_700", 9.51032/1000, ROOT.kRed+4, "")

final_zg = zg + znunu
final_zg.color = rwth.myRed
final_zg.label = "#gammaZ"
final_wg = wg + wjets
final_wg.color = rwth.myOrange
final_wg.label = "#gammaW"
final_tg = ttjets_ht+ttg
final_tg.color = rwth.myBlue
final_tg.label = "#gammat#bar{t}"

final_z = znunu
final_z.color = rwth.myRed
final_z.label = "#it{Z}(#nu#nu)"
final_w = wjets
final_w.color = rwth.myOrange
final_w.label = "#it{W}(l#nu)"
final_t = ttjets_ht
final_t.color = rwth.myBlue
final_t.label = "#it{t#bar{t}}"

