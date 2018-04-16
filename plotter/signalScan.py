#!/usr/bin/env python2

from include import *
import DatacardParser
import multiprocessing
import glob

from finalPredictions import metHist

def proceedWithWeakScan(outputDir, dset, selection, saveName, combi, onlyHigh):
    scanName = "{}_{}_{}".format(dset.label, selection, saveName)
    if combi: scanName += "_"+combi
    if onlyHigh: scanName += "_highHTG"
    xsecFile = dset.xsecs[0]
    scanRes = {}
    for fname in glob.glob("{}/*.txt.limit".format(outputDir)):
        m = re.match(".*[^\d](\d+)_(\d+).txt.limit", fname)
        with open(fname) as f: scanRes[int(m.group(1))] = limitTools.infoFromOut(f.read())

    defaultGr = ROOT.TGraph(len(scanRes))
    graphs = dict( (x,defaultGr.Clone(x)) for x in ["obs","exp","exp1up","exp1dn","exp2up","exp2dn"] )
    obsGr = ROOT.TGraph()
    xsecGr = ROOT.TGraph()
    xsecGrUp = ROOT.TGraph()
    xsecGrDn = ROOT.TGraph()
    exp1sigma = ROOT.TGraphAsymmErrors()
    exp2sigma = ROOT.TGraphAsymmErrors()
    for i, m in enumerate(sorted(scanRes)):
        xsec, xsec_unc = aux.getXsecInfoSMS(m, xsecFile)
        xsecGr.SetPoint(i, m, xsec)
        xsecGrUp.SetPoint(i, m, xsec*(1+xsec_unc))
        xsecGrDn.SetPoint(i, m, xsec*(1-xsec_unc))
        obsGr.SetPoint(i, m, xsec*scanRes[m]["obs"])
        expR = scanRes[m]["exp"]
        exp1sigma.SetPoint(i, m, xsec*expR)
        exp2sigma.SetPoint(i, m, xsec*expR)
        exp1sigma.SetPointEYhigh(i, xsec*(scanRes[m]["exp1up"]-expR))
        exp2sigma.SetPointEYhigh(i, xsec*(scanRes[m]["exp2up"]-expR))
        exp1sigma.SetPointEYlow(i, xsec*(expR-scanRes[m]["exp1dn"]))
        exp2sigma.SetPointEYlow(i, xsec*(expR-scanRes[m]["exp2dn"]))
        for name in graphs:
            graphs[name].SetPoint(i, m, xsec*scanRes[m][name] )
    writeDict(graphs, outputDir+"/Graphs2d.root")

    # beautify
    obsGr.SetLineWidth(2)

    for g in xsecGr, xsecGrUp, xsecGrDn:
        g.SetLineColor(ROOT.kBlue)
    xsecGrUp.SetLineStyle(2)
    xsecGrDn.SetLineStyle(2)

    exp2sigma.SetFillColor(ROOT.kOrange)
    exp2sigma.SetLineColor(exp2sigma.GetFillColor())
    exp1sigma.SetFillColor(ROOT.kGreen+1)
    exp1sigma.SetLineColor(2)
    exp1sigma.SetLineStyle(2)
    exp1sigma.SetLineWidth(2)

    lsp_ = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{}}#scale[0.85]{_{1}}"
    lsp_0 = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
    lsp_pm = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"

    scanParticle = lsp_0 if "NG" in scanName else lsp_
    exp2sigma.SetTitle(";m_{{{}}} (GeV);95% CL cross section upper limit (pb)".format(scanParticle))
    h = exp2sigma.GetHistogram()
    h.GetXaxis().SetRangeUser(300,1300)
    h.SetMaximum(2)
    h.SetMinimum(10e-4)
    h.Draw()

    # draw
    exp2sigma.Draw("3")
    exp1sigma.Draw("3")
    exp1sigma.Draw("xc")

    xsecGr.Draw("xc")
    xsecGrUp.Draw("xc")
    xsecGrDn.Draw("xc")

    obsGr.Draw("xcp")
    h.Draw("same")

    # legend
    exp1sigmaClone = exp1sigma.Clone()
    exp1sigmaClone.SetLineColor(exp1sigma.GetFillColor())
    exp1sigmaClone.SetLineStyle(1)
    leg = ROOT.TLegend(.45,.59,.94,.92)
    leg.SetFillStyle(0)
    leg.AddEntry(obsGr, "Observed limit", "l")
    leg.AddEntry(exp2sigma, "Expected limit #pm 1(2) #sigma_{experiment}", "f")
    leg.AddEntry(xsecGr, "Signal cross section #pm #sigma_{theory}", "l")
    leg.Draw()

    leg2 = ROOT.TLegend(.45,.66,.94,.85)
    leg2.SetFillStyle(0)
    leg2.AddEntry(None, "", "")
    leg2.AddEntry(exp1sigmaClone, "", "f")
    leg2.AddEntry(None, "", "")
    leg2.Draw()

    leg3 = leg.Clone()
    leg3.Clear()
    leg3.AddEntry(None, "", "")
    leg3.AddEntry(exp1sigma, "", "l")
    leg3.AddEntry(None, "", "")
    leg3.Draw()

    leg4 = ROOT.TLegend(.45,.63,.94,.662)
    leg4.SetFillStyle(0)
    leg4.AddEntry(xsecGrUp, "", "l")
    leg4.AddEntry(xsecGrUp, "", "l")
    leg4.Draw()

    t = ROOT.TLatex()
    if "WG" in scanName:
        info = "pp#rightarrow%s(#tilde{G}#gamma) %s(#tilde{G}W^{#pm})"%(lsp_0,lsp_pm)
        t.DrawLatexNDC(.18,.18, info)
    if "NG" in scanName:
        t.DrawLatexNDC(.18,.3, "pp#rightarrow{1}{1}/{1}{0}".format(lsp_0,lsp_pm))
        t.DrawLatexNDC(.18,.25, "{}#rightarrow{}+soft".format(lsp_pm,lsp_0))
        t.DrawLatexNDC(.18,.18, "{0}(#tilde{{G}}#gamma) {0}(#tilde{{G}}H/Z)".format(lsp_0))
    aux.Label()
    ROOT.gPad.SetLogy()
    aux.save("{}_limit".format(scanName), log=False, endings=[".pdf", ".root"])

def getPointFromDir(name):
    m = re.match("(\d*)_(\d*)", name)
    return int(m.group(1)), int(m.group(2))

def writeDict(d, filename):
    f = ROOT.TFile(filename, "recreate")
    for name, ob in d.iteritems():
        if ob: ob.Write(name)
    f.Close()

def readDict( filename ):
    f = ROOT.TFile( filename )
    tmpDir = f.GetDirectory( path )
    d = {}
    for element in tmpDir.GetListOfKeys():
        obj = element.ReadObj()
        obj = ROOT.gROOT.CloneObject( obj )
        d[element.GetName()] = obj
    return d

def writeSMSLimitConfig(infile, configName):
    pubStr = "Private Work"
    text = """
HISTOGRAM {0} obs_hist
EXPECTED {0} exp exp1up exp1dn kRed kOrange
OBSERVED {0} obs obs1up obs1dn kBlack kGray
PRELIMINARY {2}
LUMI {1:.1f}
ENERGY 13
""".format(infile,aux.intLumi/1e3,pubStr)
    with open(configName, "w+") as f:
        f.write(text)

def getXsecFile(name):
    xsecFile = ""
    if "T5" in name: xsecFile = "data/xSec_SMS_Gluino_13TeV.pkl"
    elif "T6" in name: xsecFile = "data/xSec_SMS_Squark_13TeV.pkl"
    elif "TChiWG" in name: xsecFile = "data/xSec_SMS_N2C1_13TeV.pkl"
    elif "TChiNG" in name: xsecFile = "data/xSec_SMS_TChiNG_13TeV.pkl"
    elif "GMSB" in name: xsecFile = "data/xSec_SMS_TChiNG_13TeV_interpolated_5.pkl"
    else: print "Do not know which cross section belongs to", name
    return xsecFile

def getMultiScanName(inputSignal):
    return os.path.basename(inputSignal).split("_")[0][4:]

def getScanName(inputSignal, combi):
    return getMultiScanName(inputSignal).replace("Wg", combi)


def getSignalRegionHist(dset, name):
    binning = [350, 450, 600, 900]
    hlow = metHist(dset, name, binning, True)
    high = metHist(dset, name, binning, False)
    hOut = ROOT.TH1F("", "", 6, 0, 6)
    for i in range(1,4):
        hOut.SetBinContent(i, hlow.GetBinContent(i))
        hOut.SetBinError(i, hlow.GetBinError(i))
    for i in range(1,4):
        hOut.SetBinContent(i+3, high.GetBinContent(i))
        hOut.SetBinError(i+3, high.GetBinError(i))
    # set number of unselected events
    hOut.SetBinContent(0, dset.getHist(name).Integral(0,-1, 0, -1) - hOut.Integral(0,-1))

    return hOut


def getSignalUncertainties(dset, selection, point, combi):
    out = {}

    hNominal = getSignalRegionHist(dset, "{}/{}/nominal".format(point, selection))

    hGenMet = getSignalRegionHist(dset, "{}/{}/systematics/gen".format(point, selection))
    out["genMet"] = [1.+abs(hNominal.GetBinContent(b)-hGenMet.GetBinContent(b))/2/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,7)]

    hPuUp = getSignalRegionHist(dset, "{}/{}/systematics/nopuUp".format(point, selection))
    hPuDn = getSignalRegionHist(dset, "{}/{}/systematics/nopuDn".format(point, selection))

    hPuUp.Scale(hNominal.Integral(0,-1)/hPuUp.Integral(0,-1))
    hPuDn.Scale(hNominal.Integral(0,-1)/hPuDn.Integral(0,-1))
    out["pu"] = [1.+abs(hPuUp.GetBinContent(b)-hPuDn.GetBinContent(b))/2/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,7)]

    hJesUp = getSignalRegionHist(dset, "{}/{}/systematics/jesUp".format(point, selection))
    hJesDn = getSignalRegionHist(dset, "{}/{}/systematics/jesDn".format(point, selection))
    hJerUp = getSignalRegionHist(dset, "{}/{}/systematics/jerUp".format(point, selection))
    hJerDn = getSignalRegionHist(dset, "{}/{}/systematics/jerDn".format(point, selection))
    out["jes"] = [1.+math.sqrt((abs(hJesUp.GetBinContent(b)-hJesDn.GetBinContent(b))/2/hNominal.GetBinContent(b))**2 \
        + (abs(hJerUp.GetBinContent(b)-hJerDn.GetBinContent(b))/2/hNominal.GetBinContent(b))**2) if hNominal.GetBinContent(b) else 1 for b in range(1,7)]

    hIsr = getSignalRegionHist(dset, "{}/{}/systematics/isr".format(point, selection))
    hIsrUp = getSignalRegionHist(dset, "{}/{}/systematics/isrUp".format(point, selection))
    hIsrDn = getSignalRegionHist(dset, "{}/{}/systematics/isrDn".format(point, selection))
    hIsr.Scale(hNominal.Integral(0,-1)/hIsr.Integral(0,-1))
    hIsrUp.Scale(hNominal.Integral(0,-1)/hIsrUp.Integral(0,-1))
    hIsrDn.Scale(hNominal.Integral(0,-1)/hIsrDn.Integral(0,-1))
    out["isr"] = [1.+abs(hIsrUp.GetBinContent(b)-hIsrDn.GetBinContent(b))/2/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,7)]
    if "TChi" in dset.label or "GMSB" in dset.label: # weak scan
        out["isr"] = [1.+abs(hIsr.GetBinContent(b)-hNominal.GetBinContent(b))/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,7)]

    scaleHists = dict([(i,getSignalRegionHist(dset, "{}/{}/weights/{}".format(point, selection, i))) for i in range(1,9)])
    for hname, h in scaleHists.iteritems():
        h.Scale(hNominal.Integral(0,-1)/h.Integral(0,-1))
    out["scale"] = [1.+(max([h.GetBinContent(b) for h in scaleHists.values()])-min([h.GetBinContent(b) for h in scaleHists.values()]))/2./hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,7)]


    jCRScale = 0.88 # todo: update
    jCRNorm = 1e-1 # todo: update, and also different for low and high htg
    hCont = getSignalRegionHist(dset, "{}/{}/jCR/{:.6f}".format(point, selection, jCRScale))
    hCont.Scale(jCRNorm)
    cont = [hCont.GetBinContent(b) for b in range(1,7)]

    if combi == "gg":
        hNominal = getSignalRegionHist(dset, "{}/{}/nominalGG".format(point, selection))
    elif combi == "wg":
        hNominal = getSignalRegionHist(dset, "{}/{}/nominalWG".format(point, selection))
    elif combi:
        print "ERROR: Combination '{}' is not implemented".format(combi)

    acc = [hNominal.GetBinContent(b) for b in range(1,7)]
    out["stat"] = [1.+hNominal.GetBinError(b)/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,7)]

    return acc, cont, out

def writeDataCards(outputDir, dset, selection, dataDatacard, combi, onlyHigh):
    dirs = aux.getDirNames(dset.files[0])
    #dirs = ["1600_100"] # cross check
    dc = limitTools.MyDatacard(dataDatacard)
    binNames = ["binlowEMHT_24", "binlowEMHT_25", "binlowEMHT_26", "binhighEMHT_24", "binhighEMHT_25", "binhighEMHT_26"]
    if onlyHigh:
        binNames = ["binhighEMHT_24", "binhighEMHT_25", "binhighEMHT_26"]
        dc.removeBin("binlowEMHT_24")
        dc.removeBin("binlowEMHT_25")
        dc.removeBin("binlowEMHT_26")

    for d in dirs:
        m1, m2 = getPointFromDir(d)
        if dset == t5wg and m1 < 1400: continue
        if dset == t6wg and m1 < 1100: continue
        acc, cont, syst = getSignalUncertainties(dset, selection, d, combi)
        # subtract signal contamination
        acc = [a-b for a,b in zip(acc,cont)]

        if onlyHigh:
            acc = acc[3:]
            cont = cont[3:]
            for a,b in syst.iteritems():
                syst[a] = b[3:]
        # Assume 1/4 gg, 1/2 wg and 1/4 ww in the scan, correct for those fractions
        if combi == "gg": acc = [4.*a for a in acc]
        elif combi == "wg": acc = [2.*a for a in acc]
        elif combi: print "ERROR: Do not know what to do with combination '{}'".format(combi)

        systs = {}
        for sName, valueList in syst.iteritems():
            if sName == "stat":
                for binName, unc in zip(binNames, valueList):
                    systs["signalStat_{}".format(binName)] = {binName:unc}
            else:
                systs[sName] = dict(zip(binNames, valueList))

        dc.newSignal( dict(zip(binNames, acc)), systs)
        dcName = "{}/{}.txt".format(outputDir, d)
        dc.write(dcName)

def callMultiCombine(outputDir):
    files = glob.glob("{}/*.txt".format(outputDir))
    p = multiprocessing.Pool()
    p.map(limitTools.callCombine, files)

def callMultiSignificance(outputDir):
    files = glob.glob("{}/*.txt".format(outputDir))
    p = multiprocessing.Pool()
    p.map(limitTools.callCombineSignificance, files)

def clearWrongCombineOutputs(outputDir):
    files = glob.glob("{}/*.txt.limit".format(outputDir))
    for fname in files:
        with open(fname) as f:
            rInfo = limitTools.infoFromOut(f.read())
        if rInfo["obs"] == 0:
            os.remove(fname)

def build2dGraphsLimit(outputDir):
    files = glob.glob("{}/*.txt.limit".format(outputDir))
    defaultGr = ROOT.TGraph2D(len(files))
    graphs = dict( (x,defaultGr.Clone(x)) for x in ["obs","exp","exp1up","exp1dn","exp2up","exp2dn"] )
    for g in graphs.values(): g.SetDirectory(0)
    for ifile, _file in enumerate(files):
        m = re.match(".*[^\d](\d+)_(\d+).txt.limit", _file)
        m1 = int(m.group(1))
        m2 = int(m.group(2))
        with open(_file) as f:
            rInfo = limitTools.infoFromOut(f.read())
        for name, gr in graphs.iteritems():
            graphs[name].SetPoint(ifile, m1, m2, rInfo[name])
    writeDict(graphs, outputDir+"/saved_graphs2d_limit.root")
    return graphs

def latexScanName(scanName):
    if scanName == "T6gg":
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        return "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma#tilde{G}"%(lsp_s,lsp_s)
    elif scanName == "T6Wg":
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-3.8]{#scale[0.85]{_{1}}}  "
        return "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}"%(lsp_s,lsp_s)
    elif scanName == "T5gg":
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        return "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma#tilde{G}"%(lsp_s,lsp_s)
    elif scanName == "T5Wg":
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-3.8]{#scale[0.85]{_{1}}}  "
        return "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}"%(lsp_s,lsp_s)
    return scanName

def drawSignalScanHist(h, scanName, saveName):
    smsScan = sms.sms(scanName)
    style.style2d()
    h.SetTitle("")
    if "T5" in scanName: h.GetXaxis().SetTitle("m_{#tilde{g}} (GeV)")
    if "T6" in scanName: h.GetXaxis().SetTitle("m_{#tilde{q}} (GeV)")
    if "gg" in scanName: h.GetYaxis().SetTitle("m_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)")
    if "Wg" in scanName: h.GetYaxis().SetTitle("m_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-3.8]{#scale[0.85]{_{1}}}  } (GeV)")
    hname = h.GetName()
    if hname == "significance":
        h.SetMinimum(-3)
        h.SetMaximum(3)
        style.setPaletteBWR()
        #ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
        h.GetXaxis().SetRangeUser(smsScan.Xmin,smsScan.Xmax)
        h.GetZaxis().SetNdivisions(6,0,0)
        h.GetZaxis().SetTitle("Significance (s.d.)")
    elif hname == "xsec":
        h.Scale(1000)
        h.GetZaxis().SetTitle("signal cross section (fb)")
    elif "sigComparisonAllPhi" in hname:
        var, emhtSel, metSel = h.GetName().split("_")
        emhtSel = emhtSel.replace("binlowEMHT", "#it{H}_{T}^{#gamma}<2TeV").replace("binhighEMHT", "#it{H}_{T}^{#gamma}>2TeV")
        metSel = metSel.replace("24", "350<#it{p}_{#lower[-.2]{T}}^{miss}<450GeV").replace("25", "450<#it{p}_{#lower[-.2]{T}}^{miss}<600GeV").replace("26", "#it{p}_{#lower[-.2]{T}}^{miss}>600GeV")
        h.GetZaxis().SetTitle("#frac{{Signif. |#Delta#Phi|>0.3}}{{Significance}} ({} {})".format(emhtSel, metSel))
        h.SetTitleSize(h.GetZaxis().GetTitleSize()*.9, "z")
        h.SetMinimum(.5)
        h.SetMaximum(1.5)
        style.setPaletteBWR()
    elif "signalContamination" in hname:
        var, emhtSel, metSel = h.GetName().split("_")
        emhtSel = emhtSel.replace("binlowEMHT", "#it{H}_{T}^{#gamma}<2TeV").replace("binhighEMHT", "#it{H}_{T}^{#gamma}>2TeV")
        metSel = metSel.replace("24", "350<#it{p}_{#lower[-.2]{T}}^{miss}<450GeV").replace("25", "450<#it{p}_{#lower[-.2]{T}}^{miss}<600GeV").replace("26", "#it{p}_{#lower[-.2]{T}}^{miss}>600GeV")
        h.GetZaxis().SetTitle("Sig. cont. ({} {}) (%)".format(emhtSel, metSel))
        h.Scale(100)
    else:
        uncerts = {"isr": "ISR unc.", "jes": "Jet energy unc.", "pu": "Pile-up unc.", "scale": "Renorm.+Factor. scale unc.", "genMet": "FastSim #it{p}_{T}^{miss} unc.", "totalExp": "exp. unc.", "signalStat": "stat. unc."}
        var, emhtSel, metSel = h.GetName().split("_")
        emhtSel = emhtSel.replace("binlowEMHT", "#it{H}_{T}^{#gamma}<2TeV").replace("binhighEMHT", "#it{H}_{T}^{#gamma}>2TeV")
        metSel = metSel.replace("24", "350<#it{p}_{#lower[-.2]{T}}^{miss}<450GeV").replace("25", "450<#it{p}_{#lower[-.2]{T}}^{miss}<600GeV").replace("26", "#it{p}_{#lower[-.2]{T}}^{miss}>600GeV")
        var = uncerts[var] if var in uncerts else var
        h.Scale(100)
        h.GetZaxis().SetTitle("{} ({} {}) (%)".format(var, emhtSel, metSel))
        if var in uncerts or var in uncerts.values():
            h.SetMinimum(0)
            #h.SetMaximum(15)
    h.GetZaxis().SetTitleOffset(1.42)
    h.GetYaxis().SetTitleOffset(h.GetYaxis().GetTitleOffset()*.95)
    h.SetLabelSize(h.GetLabelSize()*0.85, "xyz")

    c = ROOT.TCanvas()
    h.Draw("colz")
    aux.drawDiagonal(h)
    if hname == "significance":
        arxivText = "arXiv:1707.06193"
        text = ROOT.TLatex()
        text.SetTextFont(82)
        text.SetTextSize(0.85*text.GetTextSize())
        text.DrawTextNDC(.19, .84, arxivText)
        l = aux.Label2D(info="#scale[.76]{{{}}}".format(latexScanName(scanName)), status="Supplementary", drawAll=False)
        l.lum = ROOT.TLatex( .51, .95, "35.9 fb^{-1} (13 TeV)" )
        l.draw()
    else:
        sn1, sn2, sn3 = latexScanName(scanName).split(",")
        label = aux.Label2D()
        label.info.DrawLatexNDC(label.info.GetX(), label.info.GetY()-.02, ",".join([sn1,sn2]))
        label.info.DrawLatexNDC(label.info.GetX(), label.info.GetY()-.08, " "+sn3.strip())

    aux.save("{}_{}".format(scanName,saveName))
    style.defaultStyle()

def drawLimitInput1d(outputDir, scanName, xsecFile):
    files = glob.glob("{}/*.txt".format(outputDir))
    files = sorted(files, key=lambda f: int(f.split("_")[-2]))
    dc = limitTools.MyDatacard(files[0])
    defaultGr = ROOT.TGraph(len(files))
    graphs = dict((x,defaultGr.Clone(x)) for x in ["Acceptance_"+b for b in dc.bins]+["xsec"])
    uncerts = "isr", "jes", "pu", "scale", "genMet", "signalStat", "totalExp"
    graphs.update(dict([(x+"_"+y,defaultGr.Clone(x+"_"+y)) for x in uncerts for y in dc.bins]))
    for ifile, f in enumerate(files):
        m = re.match(".*_(\d+)_(\d+).txt", f)
        m1 = int(m.group(1))
        xsec = aux.getXsecInfoSMS(m1, xsecFile)[0]
        graphs["xsec"].SetPoint(ifile, m1, xsec)
        dc = limitTools.MyDatacard(f)
        systDict = dict([(l[0],l) for l in dc.systs])
        for b in dc.bins:
            graphs["Acceptance_"+b].SetPoint(ifile, m1, dc.exp[b]["signal"]/(xsec*aux.intLumi))
            for unc in uncerts:
                if unc == "signalStat":
                    value = systDict[unc+"_"+b][4][b]["signal"] if unc+"_"+b in systDict else 1
                elif unc == "totalExp":
                    totalUncertsSources = "dataMC", "genMet", "isr", "pu", "jes", "scale", "signalStat_"+b, "trigger", "lumi"
                    value = math.sqrt(sum([(systDict[s][4][b]["signal"]-1)**2 for s in totalUncertsSources]))+1 if "signalStat_"+b in systDict else 1
                else:
                    value = systDict[unc][4][b]["signal"]
                graphs[unc+"_"+b].SetPoint(ifile, m1, value-1)
    for name, gr in graphs.iteritems():
        c = ROOT.TCanvas()
        var, emhtSel, metSel = name.split("_")
        emhtSel = emhtSel.replace("binlowEMHT", "#it{H}_{T}^{#gamma}<2TeV").replace("binhighEMHT", "#it{H}_{T}^{#gamma}>2TeV")
        metSel = metSel.replace("24", "350<#it{p}_{#lower[-.2]{T}}^{miss}<450GeV").replace("25", "450<#it{p}_{#lower[-.2]{T}}^{miss}<600GeV").replace("26", "#it{p}_{#lower[-.2]{T}}^{miss}>600GeV")
        gr.SetTitle(";#it{m}_{#chi} (GeV);uncertainty (%)"+var)
        gr.Draw()
        l = aux.Label(info=scanName+" "+metSel+emhtSel)
        aux.save("uncertainty_{}_{}".format(scanName,name), log=False)

def drawLimitInput(outputDir, scanName):
    if "TChi" in scanName: return drawLimitInput1d(outputDir, scanName, xsecFile)
    files = glob.glob("{}/*.txt".format(outputDir))
    dc = limitTools.MyDatacard(files[0])
    defaultGr = ROOT.TGraph2D(len(files))
    graphs = dict((x,defaultGr.Clone(x)) for x in ["Acceptance_"+b for b in dc.bins]+["xsec"])
    uncerts = "isr", "jes", "pu", "scale", "genMet", "signalStat", "totalExp"
    graphs.update(dict([(x+"_"+y,defaultGr.Clone(x+"_"+y)) for x in uncerts for y in dc.bins]))
    for ifile, f in enumerate(files):
        m = re.match(".*_(\d+)_(\d+).txt", f)
        m1 = int(m.group(1))
        m2 = int(m.group(2))
        xsec = aux.getXsecInfoSMS(m1, xsecFile)[0]
        graphs["xsec"].SetPoint(ifile, m1, m2, xsec)
        dc = limitTools.MyDatacard(f)
        systDict = dict([(l[0],l) for l in dc.systs])
        for b in dc.bins:
            graphs["Acceptance_"+b].SetPoint(ifile, m1, m2, dc.exp[b]["signal"]/(xsec*aux.intLumi))
            for unc in uncerts:
                if unc == "signalStat":
                    value = systDict[unc+"_"+b][4][b]["signal"] if unc+"_"+b in systDict else 1
                elif unc == "totalExp":
                    totalUncertsSources = "dataMC", "genMet", "isr", "pu", "jes", "scale", "signalStat_"+b, "trigger", "lumi"
                    value = math.sqrt(sum([(systDict[s][4][b]["signal"]-1)**2 for s in totalUncertsSources]))+1 if "signalStat_"+b in systDict else 1
                else:
                    value = systDict[unc][4][b]["signal"]
                graphs[unc+"_"+b].SetPoint(ifile, m1, m2, value-1)
    for name, gr in graphs.iteritems():
        h = gr.GetHistogram()
        drawSignalScanHist(h, scanName, name)

def getXsecLimitHistDelaunay(gr):
    grScaled = scaleObsWithXsec(gr)
    grScaled.SetNpx(500)
    grScaled.SetNpy(500)
    h = grScaled.GetHistogram()
    h = h.Clone(aux.randomName())
    return h

def build2dGraphs(outputDir):
    build2dGraphsLimit(outputDir)

def getObsUncertainty(gr2d, xsecFile):
    gr2dup = gr2d.Clone(aux.randomName())
    gr2ddn = gr2d.Clone(aux.randomName())
    gr2dup.SetDirectory(0)
    gr2ddn.SetDirectory(0)
    points = [ (gr2d.GetX()[i], gr2d.GetY()[i], gr2d.GetZ()[i]) for i in range(gr2d.GetN()) ]
    for ip, (x,y,z) in enumerate(points):
        xsec, unc = aux.getXsecInfoSMS(x, xsecFile)
        gr2dup.SetPoint(ip, x, y, z*(1-unc/100))
        gr2ddn.SetPoint(ip, x, y, z*(1+unc/100))
    obsUp = limitTools.getContour(gr2dup)
    obsUp.SetName("obs1up")
    obsDn = limitTools.getContour(gr2ddn)
    obsDn.SetName("obs1dn")
    return {"obs1up": obsUp, "obs1dn": obsDn}

def interpolateAlongY(h2):
    h2Out = h2.Clone(aux.randomName())
    for xbin, ybin in aux.loopH(h2):
        if not h2.GetBinContent(xbin,ybin):
            ybinUp, cUp = 0, 0
            for ybinUp in range(ybin,h2.GetNbinsY()):
                cUp = h2.GetBinContent(xbin,ybinUp)
                if cUp: break
            ybinDn, cDn = 0, 0
            for ybinDn in range(ybin,0, -1):
                cDn = h2.GetBinContent(xbin,ybinDn)
                if cDn: break
            if not cUp or not cDn: continue
            h2Out.SetBinContent(xbin,ybin, ((cUp-cDn)*ybin+(cDn*ybinUp-cUp*ybinDn))/(ybinUp-ybinDn)) # linear interpolation
    return h2Out

def interpolateHoles(h2):
    for xbin, ybin in aux.loopH(h2):
        c = h2.GetBinContent(xbin, ybin)
        if abs(c)>1e-5: continue
        cNord = h2.GetBinContent(xbin, ybin+1)
        cSout = h2.GetBinContent(xbin, ybin-1)
        cWest = h2.GetBinContent(xbin+1, ybin)
        cEast = h2.GetBinContent(xbin-1, ybin)
        if cNord and cSout and cWest and cEast:
            h2.SetBinContent(xbin, ybin, sum([cNord,cSout,cWest,cEast])/4.)
    return h2

def getXsecLimitHist( gr2d, h, xsecFile ):
    h.SetDirectory(0) # or the next line will overwrite the hist?
    points = [ (gr2d.GetX()[i], gr2d.GetY()[i], gr2d.GetZ()[i]) for i in range(gr2d.GetN()) ]
    for x,y,z in points:
        xsec = aux.getXsecInfoSMS(x, xsecFile)[0]
        h.SetBinContent(h.FindBin(x,y), z*xsec )
    return h

def scaleGraph2dByXsec( gr2d, xsecFile ):
    out = gr2d.Clone(aux.randomName())
    for i in range(gr2d.GetN()):
        x = gr2d.GetX()[i]
        z = gr2d.GetZ()[i]
        xsec = aux.getXsecInfoSMS(x, xsecFile)[0]
        out.SetPoint(i, x, gr2d.GetY()[i], xsec*z)
    return out

def getHistForModel( model ):
    h = ROOT.TH2F()
    if "T5" in model: h = ROOT.TH2F("", "", 35, 775, 2525, 500, 0, 2500)
    elif "T6" in model: h = ROOT.TH2F("", "", 24, 975, 2175, 220, 0, 2200)
    else: print "Not specified model", model
    smsScan = sms.sms(model)
    h.SetTitle("{};{};{};".format(model, smsScan.sParticle, smsScan.LSP))
    h.SetMinimum(0)
    return h

def smoothContour(gr, neighbors=5, sigma=.5):
    fgaus = ROOT.TF1("fgaus", "gaus", -10, 10)
    fgaus.SetParameters(1,0,1)
    weights = [fgaus.Eval(i*sigma) for i in range(neighbors)]
    out = gr.Clone(aux.randomName())
    out.Set(0)
    n = gr.GetN()
    Xs = [gr.GetX()[i] for i in range(n)]
    Ys = [gr.GetY()[i] for i in range(n)]
    n = Ys.index(max(Ys))+1
    Xs = Xs[0:n]
    Ys = Ys[0:n]
    for i, (x, y) in enumerate(zip(Xs,Ys)):
        pNeigh = min(neighbors, i+1, n-i)
        newX, ws, newY = 0, 0, 0
        for j in range(pNeigh):
            if j:
                newX += (Xs[i-j]+Xs[i+j])*weights[j]
                newY += (Ys[i-j]+Ys[i+j])*weights[j]
                ws += 2*weights[j]
            else:
                newX += x*weights[0]
                newY += y*weights[0]
                ws += weights[0]
        out.SetPoint(i, newX/ws, newY/ws)
    return out

def build1dGraphs(outputDir, dset):
    graphs = readDict(outputDir+"/saved_graphs2d_limit.root")
    toDraw = dict( [(name,limitTools.getContour(gr)) for name,gr in graphs.iteritems() ] )
    toDraw.update(getObsUncertainty(graphs["obs"], dset.xsecs[0]))
    toDraw = dict([(n,smoothContour(gr)) for n,gr in toDraw.iteritems()])
    #toDraw["obs_hist"] = getXsecLimitHistDelaunay(graphs["obs"])
    toDraw["obs_gr_scaled"] = scaleGraph2dByXsec(graphs["obs"], dset.xsecs[0])
    toDraw["obs_hist_gaps"] = getXsecLimitHist( graphs["obs"], getHistForModel(dset.label), dset.xsecs[0] )
    toDraw["obs_hist"] = interpolateAlongY(toDraw["obs_hist_gaps"])
    writeDict(toDraw, outputDir+"/saved_graphs1d_limit.root")

def signalScan(dset, selection, dataDatacard, saveName="", combi="", onlyHigh=False):

    outputDir = "limitCalculations/{}_{}".format(dset.label, selection)
    if saveName: outputDir += "_" + saveName
    if combi: outputDir += "_" + combi
    if onlyHigh: outputDir += "_highHtg"
    if not os.path.isdir(outputDir): os.mkdir(outputDir)
    if True and False:
        writeDataCards(outputDir, dset, selection, dataDatacard, combi, onlyHigh)
        if dset == t5wg: writeDataCards(outputDir, t5wg_ext, selection, dataDatacard, combi, onlyHigh)
        if dset == t6wg: writeDataCards(outputDir, t6wg_ext, selection, dataDatacard, combi, onlyHigh)
    #drawLimitInput(outputDir, dset.label)
    #callMultiCombine(outputDir)
    ##clearWrongCombineOutputs(outputDir)
    if "TChi" in dset.label: return proceedWithWeakScan(outputDir, dset, selection, saveName, combi, onlyHigh)
    build2dGraphs(outputDir)
    build1dGraphs(outputDir, dset)
    writeSMSLimitConfig(outputDir+"/saved_graphs1d_limit.root", outputDir+"/smsPlotterer.cfg")
    subprocess.call(["python2", "smsPlotter/python/makeSMSplots.py", outputDir+"/smsPlotterer.cfg", "plots/limits_{}_".format(outputDir.replace("limitCalculations/", ""))])

if __name__ == "__main__":
    #signalScan(t5wg, "original", "dataCards/final_original.txt", "test")
    #signalScan(t5wg, "original", "dataCards/final_original.txt", "test", "gg")
    #signalScan(t6wg, "original", "dataCards/final_original.txt", "test")
    #signalScan(t6wg, "original", "dataCards/final_original.txt", "test", "gg")

    #signalScan(t5wg, "di_cleaned", "dataCards/final_di_cleaned.txt", "test")
    #signalScan(t5wg, "lep_cleaned", "dataCards/final_lep_cleaned.txt", "test")
    #signalScan(t5wg, "dilep_cleaned", "dataCards/final_dilep_cleaned.txt", "test")
    #signalScan(t5wg, "st_cleaned", "dataCards/final_st_cleaned.txt", "test")
    #signalScan(t5wg, "all_cleaned", "dataCards/final_all_cleaned.txt", "test")
    signalScan(t5wg, "dilep_cleaned", "dataCards/final_dilep_cleaned.txt", "test", onlyHigh=True)

    # repeat with tching
    #signalScan(tching, "original", "dataCards/final_original.txt", "test")
    #signalScan(tching, "dilep_cleaned", "dataCards/final_dilep_cleaned.txt", "test")
    #signalScan(tching, "all_cleaned", "dataCards/final_all_cleaned.txt", "test")
    #signalScan(tching, "dilep_cleaned", "dataCards/final_dilep_cleaned.txt", "test", onlyHigh=True)
