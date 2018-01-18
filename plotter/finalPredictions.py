#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

def metHist(dataset, name, binning=None, lowEmht=True):
    h = dataset.getHist(name)
    if not h: return
    if isinstance(h, ROOT.TH2):
        if lowEmht: h = h.ProjectionX(aux.randomName(), 1, 1)
        else:       h = h.ProjectionX(aux.randomName(), 2, 2)
    if binning: h = aux.rebin(h, binning)
    aux.appendFlowBin(h, False)
    h.SetYTitle(aux.getYAxisTitle(h))
    return h

def pdfUncertainty(dataset, lowEmht, dirName, nBins, saveName=""):
    hNominal = metHist(dataset, "0_0/{}/nominal".format(dirName), nBins, lowEmht)
    hList = []
    for iWeight in range(9,110):
        h = metHist(dataset, "0_0/{}/weights/{}".format(dirName, iWeight), nBins, lowEmht)
        hList.append(h.Clone(aux.randomName()))
    hUp, hDn = aux.getEnvelopeHists(hList)
    hSys = aux.getSystFromEnvelopes(hNominal, hUp, hDn)
    hSys = aux.getSystFromVariance(hNominal, hList)

    if saveName:
        c = ROOT.TCanvas()
        m = multiplot.Multiplot()
        hNominal.SetLineColor(ROOT.kBlack)
        hNominal.SetLineWidth(2)
        m.add(hNominal, "Nominal")
        for h in hList:
            h.SetLineColor(ROOT.kRed)
            h.darwOption_ = "hist"
            m.add(h, "100 PDF replicas" if h == hList[0] else "")
        m.Draw()
        denominator = hNominal.Clone(aux.randomName())
        for b in aux.loopH(denominator): denominator.SetBinError(b,0)
        r = ratio.Ratio("PDF uncert.", hNominal, denominator, hSys)
        r.draw(.5,1.5)
        l = aux.Label(sim=True, info=dataset.label)
        aux.save("pdfUncertainty_"+saveName, normal=False)
    return hSys

def puUncertainty(dataset, lowEmht, dirName, nBins, saveName=""):
    hNominal = metHist(dataset, "0_0/{}/nominal".format(dirName), nBins, lowEmht)
    hUp = metHist(dataset, "0_0/{}/systematics/puUp".format(dirName), nBins, lowEmht)
    hDn = metHist(dataset, "0_0/{}/systematics/puDn".format(dirName), nBins, lowEmht)
    hSys = aux.getSystFromEnvelopes(hNominal, hUp, hDn)

    if saveName:
        c = ROOT.TCanvas()
        hSys.SetFillColor(ROOT.kRed)
        hNominal.SetLineColor(ROOT.kBlack)
        hNominal.SetLineWidth(2)
        hSys.Draw("e2")
        hNominal.Draw("e hist same")
        leg = ROOT.TLegend(.50,.82,.94,.915)
        leg.SetFillStyle(0)
        leg.AddEntry(hSys, "Nominal #pm pileup uncert.", "f")
        leg.Draw()
        leg2 = leg.Clone()
        leg2.Clear()
        leg2.AddEntry(hNominal, "", "l")
        leg2.Draw()
        denominator = hNominal.Clone(aux.randomName())
        for b in aux.loopH(denominator): denominator.SetBinError(b,0)
        r = ratio.Ratio("pileup uncert.", hNominal, denominator, hSys)
        r.draw(.8,1.2)
        l = aux.Label(sim=True, info=dataset.label)
        aux.save("puUncertainty_"+saveName, normal=False)
    return hSys

def scaleUncertainty(dataset, lowEmht, dirName, nBins, saveName=""):
    hNominal = metHist(dataset, "0_0/{}/nominal".format(dirName), nBins, lowEmht)
    hList = []
    for iWeight in range(9):
        h = metHist(dataset, "0_0/{}/weights/{}".format(dirName, iWeight), nBins, lowEmht)
        hList.append(h.Clone(aux.randomName()))
    hUp, hDn = aux.getEnvelopeHists(hList)
    hSys = aux.getSystFromEnvelopes(hNominal, hUp, hDn)

    if saveName:
        c = ROOT.TCanvas()
        m = multiplot.Multiplot()
        hNominal.SetLineColor(ROOT.kBlack)
        m.add(hNominal)
        scales = [(1,1), (1,2), (1,.5), (2,1), (2,2), (2,.5), (.5,1), (.5,2), (.5,.5)]
        for iWeight, (mur, muf) in enumerate(scales):
            h = hList[iWeight]
            h.drawOption_ = "hist"
            if iWeight/3 == 0: h.SetLineColor(ROOT.kBlack)
            if iWeight/3 == 1: h.SetLineColor(ROOT.kRed)
            if iWeight/3 == 2: h.SetLineColor(ROOT.kBlue)
            if iWeight%3 == 1: h.SetLineStyle(2)
            if iWeight%3 == 2: h.SetLineStyle(3)
            m.add(h, "#mu_{{r}}={} #mu_{{f}}={}".format(mur,muf))
        m.Draw()
        denominator = hNominal.Clone(aux.randomName())
        for b in aux.loopH(denominator): denominator.SetBinError(b,0)
        r = ratio.Ratio("scale uncert.", hNominal, denominator, hSys)
        r.draw(.5,1.5)
        l = aux.Label(sim=True, info=dataset.label)
        aux.save("scaleUncertainty_"+saveName, normal=False)
    return hSys


def jecUncertainty(dataset, lowEmht, dirName, nBins, saveName=""):
    hNominal = metHist(dataset, "0_0/{}/nominal".format(dirName), nBins, lowEmht)
    hJesUp = metHist(dataset, "0_0/{}/systematics/jesUp".format(dirName), nBins, lowEmht)
    hJesDn = metHist(dataset, "0_0/{}/systematics/jesDn".format(dirName), nBins, lowEmht)
    hJerUp = metHist(dataset, "0_0/{}/systematics/jerUp".format(dirName), nBins, lowEmht)
    hJerDn = metHist(dataset, "0_0/{}/systematics/jerDn".format(dirName), nBins, lowEmht)

    hJesSys = aux.getSystFromEnvelopes(hNominal, hJesUp, hJesDn)
    hJerSys = aux.getSystFromEnvelopes(hNominal, hJerUp, hJerDn)
    hSys = aux.addUncertaintiesQuadratic([hJesSys,hJerSys])

    if saveName:
        c = ROOT.TCanvas()
        m = multiplot.Multiplot()
        hNominal.SetLineColor(ROOT.kBlack)
        m.add(hNominal)
        for h in hJesUp, hJesDn:
            h.SetLineColor(ROOT.kRed)
            h.drawOption_ = "hist"
        for h in hJerUp, hJerDn:
            h.SetLineColor(ROOT.kBlue)
            h.drawOption_ = "hist"
        m.add(hJesUp, "Scale +")
        m.add(hJesDn, "Scale -")
        m.add(hJerUp, "Resol. +")
        m.add(hJerDn, "Resol. -")
        m.Draw()
        denominator = hNominal.Clone(aux.randomName())
        for b in aux.loopH(denominator): denominator.SetBinError(b,0)
        r = ratio.Ratio("jet uncert.", hNominal, denominator, hSys)
        r.draw(.8,1.2)
        l = aux.Label(sim=True, info=dataset.label)
        aux.save("jecUncertainty_"+saveName, normal=False)
    return hSys


def plotOverlayedPredictions(filename):
    c = ROOT.TCanvas()
    f = ROOT.TFile(filename)
    hdir = f.Get("dir")
    hdir.GetXaxis().SetRangeUser(0,100)
    aux.drawOpt(hdir, "data")
    m = multiplot.Multiplot()
    isData = "data" in filename or "_final_" in filename
    if isData: m.add(hdir, "Data")
    else: m.add(hdir, "Direct simulation")
    histsForRatio = []
    colors = [rwth.magenta, rwth.lila, rwth.mayGreen, rwth.orange, rwth.red]
    selScales = ["1.0", "0.89", "0.9", "0.91", "0.92"]
    if "highEMHT" in filename:
        selScales = ["1.0", "0.8", "0.83", "0.86", "0.9"]
        m.maximum = 50
    elif "lowEMHT" in filename:
        m.maximum = 8000
    for iScale, scale in enumerate(selScales):
        h = f.Get(scale)
        if not isinstance(h, ROOT.TH1): continue
        h.GetXaxis().SetRangeUser(0,100)
        h.Scale(hdir.Integral()/h.Integral())
        h.SetLineColor(colors[iScale])
        h.SetMarkerColor(colors[iScale])
        if scale == "1.0":
            h.SetLineColor(ROOT.kBlue)
            m.add(h, "Jet CR")
        else:
            m.add(h, "Jet CR, scaled by "+scale)
        histsForRatio.append(h.Clone(aux.randomName()))
    m.Draw()
    r = ratio.Ratio("#frac{Scaled CR}{Data}", histsForRatio[0], hdir)
    for hist in [ r.ratio, r.ratioSys, r.ratioStat, r.totalUncert ]:
        hist.SetTitleOffset(1.4, "Y")

    if "highEMHT" in filename:
        r.draw(0.8, 1.2)
    else:
        r.draw()
    for h in histsForRatio[1:]:
        h.Divide(hdir)
        h.Draw("same")
    aux.Label(sim=not isData)
    l = ROOT.TLatex()
    if "lowEMHT" in filename: l.DrawLatexNDC(.3, .4, "#scale[.7]{#it{H}_{T}^{#gamma}<2TeV}")
    if "highEMHT" in filename: l.DrawLatexNDC(.2, .75, "#scale[.7]{#it{H}_{T}^{#gamma}>2TeV}")
    aux.save(filename.replace(".root","_overlay").replace("/","_"))

def getDirNames(filename, path):
    f = ROOT.TFile(filename)
    p = f.GetDirectory(path)
    return [key.GetName() for key in p.GetListOfKeys()]

def gjetPredictionHist(dirHist, preSet, subSet, nBins, lowEmht=True, saveName="", preDirJet="original"):

    # modify dir hist
    dirHist = dirHist.Clone("dir")
    maxBin = dirHist.FindFixBin(100)-1 # do not take met=100 bin
    dirInt, dirIntErr = aux.integralAndError(dirHist, 1, maxBin)
    for bin in [0]+range(maxBin+1, dirHist.GetNbinsX()+2):
        dirHist.SetBinContent(bin, 0)
        dirHist.SetBinError(bin, 0)


    scales = getDirNames(preSet.files[0], "0_0/{}/jCR".format(preDirJet))
    gr = ROOT.TGraph()
    for iscale, scale in enumerate(scales):
        preHist = metHist(preSet, "0_0/{}/jCR/{}".format(preDirJet, scale), nBins, lowEmht)
        if subSet:
            mcPreHist = metHist(subSet, "0_0/{}/jCR/{}".format(preDirJet, scale), nBins, lowEmht)
            preHist.Add(mcPreHist, -1)
            preHist.SetName(scale)
        for bin in [0]+range(maxBin+1, dirHist.GetNbinsX()+2):
            preHist.SetBinContent(bin, 0)
            preHist.SetBinError(bin, 0)
        preHist.Scale(dirInt/preHist.Integral(1,maxBin))
        chi2 = dirHist.Chi2Test(preHist, "OF UF CHI2 WW")
        c = ROOT.TCanvas()
        dirHist.GetXaxis().SetRangeUser(0,100)
        preHist.GetXaxis().SetRangeUser(0,100)
        dirHist.Draw()
        preHist.SetLineColor(2)
        preHist.Draw("same")
        r = ratio.Ratio("#gamma/Jet", dirHist,preHist)
        r.draw(.95,1.05)
        l = aux.Label(info="Scale={}  #chi^{{2}}={}".format(scale, chi2))
        #aux.save(saveName+"_scale{}percent".format(int(float(scale)*100)), "savedFitPredictions/",log=False)
        gr.SetPoint(iscale, float(scale), chi2)
    c = ROOT.TCanvas()
    points = [(i, gr.GetX()[i],gr.GetY()[i]) for i in range(gr.GetN())]
    ys = [gr.GetY()[i] for i in range(gr.GetN())]
    minIndex = ys.index(min(ys))
    sidePoints = 3
    if "final_highEMHT" in saveName: sidePoints = 8
    gr2 = ROOT.TGraph()
    for iNew,iOld in enumerate(range(max(0,minIndex-sidePoints),min(gr.GetN(),minIndex+sidePoints+1))):
        gr2.SetPoint(iNew, points[iOld][1], points[iOld][2])
    if not "ee_highEMHT" in saveName and not "chi2plot" in saveName:
        gr = gr2

    gr.SetTitle(";Scale;#chi^{2}")
    gr.SetMarkerStyle(20)
    maximum = max([gr.GetY()[i] for i in range(gr.GetN())])
    gr.SetMaximum(1.7*maximum)
    gr.Draw("ap")
    gr.Fit("pol2", "Q")
    fitFunc = gr.GetFunction("pol2")
    fitFunc.SetLineColor(ROOT.kRed)
    parameters = fitFunc.GetParameters()
    a0, a1, a2 = parameters[0], parameters[1], parameters[2]
    chi2AtMin = a0 - a1**2/(4*a2)
    fitScale = -a1/(2*a2)
    deltaChi2 = 1 # change chi2 by this value for the uncertainty
    fitErr = aux.sqrt(deltaChi2/a2) if a2>0 else 0
    errFunc = fitFunc.Clone()
    errFunc.SetRange(fitScale-fitErr, fitScale+fitErr)
    errFunc.SetFillColorAlpha(errFunc.GetLineColor(), .5)
    errFunc.SetFillStyle(1001)

    errFunc.Draw("FC same")
    text = ROOT.TLatex()
    text.SetTextSize(text.GetTextSize()*0.8)
    text.DrawLatexNDC(.2, .75, "Best scale = {:.3f} #pm {:.3f} (stat.) #pm {:.3f} (syst.)".format(fitScale,fitErr,abs(1-fitScale)))
    text.DrawLatexNDC(.4, .67, "#chi^{{2}}/NDF = {:.1f} / {}".format(chi2AtMin, maxBin-1))
    leg = ROOT.TLegend(.54,.82,.94,.92)
    leg.SetFillStyle(0)
    leg.AddEntry(errFunc, "Statistical uncertainty", "f")
    leg.Draw()
    info = "  #it{H}_{T}^{#gamma} %s 2 TeV"
    info = info%(">" if "highEMHT" in saveName else "<")
    aux.Label(sim="Closure" in saveName, info=info)
    aux.save("savedFitPredictions_"+saveName+"_fit", log=False)

    err = aux.sqrt(fitErr**2 + (1-fitScale)**2)
    fitScaleUp = fitScale + err
    fitScaleDn = fitScale - err

    fitScaleR = round(fitScale, 2)
    fitScaleUpR = round(fitScaleUp, 2)
    fitScaleDnR = round(fitScaleDn, 2)
    if fitScaleUpR > 1.2:
        print "Error: Fit scale Up = ", fitScaleUpR
        fitScaleUpR = 1.2
    if fitScaleDnR < 0.85:
        print "Error: Fit scale Dn = ", fitScaleDnR
        fitScaleDnR = 0.85
    if fitScaleR < 0.85 or fitScaleR > 1.2:
        print "Error: Fit scale    = ", fitScaleR
        fitScaleR = 1.00

    preHist = metHist(preSet, "0_0/{}/jCR/{:.6f}".format(preDirJet, fitScaleR), nBins, lowEmht)
    preHistUp = metHist(preSet, "0_0/{}/jCR/{:.6f}".format(preDirJet, fitScaleUpR), nBins, lowEmht)
    preHistDn = metHist(preSet, "0_0/{}/jCR/{:.6f}".format(preDirJet, fitScaleDnR), nBins, lowEmht)
    if subSet:
        mcPreHist = metHist(subSet, "0_0/{}/jCR/{:.6f}".format(preDirJet, fitScaleR), nBins, lowEmht)
        preHistMcUp = preHist.Clone(aux.randomName())
        preHistMcDn = preHist.Clone(aux.randomName())
        preHistMcUp.Add(mcPreHist, -0.7)
        preHistMcDn.Add(mcPreHist, -1.3)
        mcPreHistUp = metHist(subSet, "0_0/{}/jCR/{:.6f}".format(preDirJet, fitScaleUpR), nBins, lowEmht)
        mcPreHistDn = metHist(subSet, "0_0/{}/jCR/{:.6f}".format(preDirJet, fitScaleDnR), nBins, lowEmht)
        preHist.Add(mcPreHist, -1)
        preHistUp.Add(mcPreHistUp, -1)
        preHistDn.Add(mcPreHistDn, -1)

        # consider 30% uncertainy of mc background
        systFromMc = aux.getSystFromDifference(preHistMcDn, preHistMcUp)

    syst = aux.getSystFromDifference(preHistDn, preHistUp)
    if subSet:
        syst = aux.addUncertaintiesQuadratic([syst,systFromMc])


    # Scale
    preInt, preIntErr = aux.integralAndError(preHist, 1, maxBin)
    preIntUp, preIntErrUp = aux.integralAndError(preHistUp, 1, maxBin)
    preIntDn, preIntErrDn = aux.integralAndError(preHistDn, 1, maxBin)
    preIntErr2 = abs(preIntDn-preIntUp)/2.
    relScaleUncert = aux.sqrt( (dirIntErr/dirInt)**2 + (preIntErr/preInt)**2 + (preIntErr2/preInt)**2 )

    for h in preHist, syst:
        h.Scale(dirInt/preInt)
        h.SetDirectory(0)

    # Force symmetric systematic uncertainty
    for bin in range(syst.GetNbinsX()+2):
        syst.SetBinContent(bin, preHist.GetBinContent(bin))
        syst.SetBinError(bin, aux.sqrt(syst.GetBinError(bin)**2 + (relScaleUncert*syst.GetBinContent(bin))**2))
#    if style.divideByBinWidth:
#        preHist.Scale(1., "width")
#        syst.Scale(1., "width")

    # do not allow negative event yields
    # set to small value, see here: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit#Why_does_combine_have_trouble_wi
    for h in preHist, syst:
        for bin in aux.loopH(h):
            if h.GetBinContent(bin)<0: h.SetBinContent(bin,1e-7)

    info = {"scale": dirInt/preInt, "scaleErrRel": relScaleUncert, "shift": fitScale, "shiftErr": err, "shiftErrStat": fitErr}
    return preHist, syst, info

def getDatacardUncertFromHist(h,b):
    c = h.GetBinContent(b)
    return 1.+max(0,h.GetBinError(b)/c if c else 0)

def onlyPositiveContents(h):
    for b in aux.loopH(h):
        if h.GetBinContent(b)<0: h.SetBinContent(b,0)

def metHist(dataset, name, binning=None, lowEmht=True):
    h = dataset.getHist(name)
    if not h: return
    if isinstance(h, ROOT.TH2):
        if lowEmht: h = h.ProjectionX(aux.randomName(), 1, 1)
        else:       h = h.ProjectionX(aux.randomName(), 2, 2)
    if binning: h = aux.rebin(h, binning)
    aux.appendFlowBin(h, False)
    h.SetYTitle(aux.getYAxisTitle(h))
    return h

def finalDistributionSignalHist(name, lowEmht, dirSet, dirDir, preSet=None):
    name += "_lowEMHT" if lowEmht else "_highEMHT"
    style.divideByBinWidth = True

    nBins = range(0,200,10)+[200, 250, 300, 350, 450, 600, 800]
    if name == "electronClosure_highEMHT": nBins = [0, 50, 100, 150, 200, 250, 300, 350, 450, 600, 800]

    # direct stuff
    dirHist = metHist(dirSet, "0_0/{}/nominal".format(dirDir), nBins, lowEmht)
    if "electronClosure" in name:
        dirHist = metHist(dirSet, "0_0/{}/genE".format(dirDir), nBins, lowEmht)
    style.additionalPoissonUncertainty = False
    aux.drawOpt(dirHist, "data")

    if "electronClosure" not in name:
        if "qcdClosure" in name:
            gjetHist, gjetSyst, info = gjetPredictionHist(dirHist, preSet, None, nBins, lowEmht, name, preDirJet="original")
        else:
            gjetHist, gjetSyst, info = gjetPredictionHist(dirHist, preSet, zg+wg+ttg+wjets+ttjets_nlo+znunu, nBins, lowEmht, name, "original")
        print info
        gjetHist.SetLineColor(rwth.myLightBlue)
        gjetHist.GetXaxis().SetTitle("#it{p}_{T}^{miss} (GeV)")

    if "qcdClosure" not in name:
        eHist = metHist(dirSet, "0_0/{}/eControl".format(dirDir), nBins, lowEmht)
        eHist.GetXaxis().SetTitle("#it{p}_{T}^{miss} (GeV)")
        eHist.Scale( 0.0267 if dirSet == data else 0.0154 )
        eHist.SetLineColor(rwth.myYellow)
        eSyst = aux.getSysHisto(eHist, 0.3)

    if "Closure" not in name:
        zgHist = metHist(zg+znunu, "0_0/{}/nominal".format(dirDir), nBins, lowEmht)
        wgHist = metHist(wg+wjets, "0_0/{}/nominal".format(dirDir), nBins, lowEmht)
        tgHist = metHist(ttjets_ht+ttg, "0_0/{}/nominal".format(dirDir), nBins, lowEmht)

        zgHist.SetLineColor(rwth.myRed)
        wgHist.SetLineColor(rwth.myOrange)
        tgHist.SetLineColor(rwth.myBlue)

        zgPdfUnc = pdfUncertainty(zg+znunu, lowEmht, dirDir, nBins)
        wgPdfUnc = pdfUncertainty(wg+wjets, lowEmht, dirDir, nBins)
        tgPdfUnc = pdfUncertainty(ttjets_nlo+ttg, lowEmht, dirDir, nBins)

        zgPuUnc = puUncertainty(zg+znunu, lowEmht, dirDir, nBins)
        wgPuUnc = puUncertainty(wg+wjets, lowEmht, dirDir, nBins)
        tgPuUnc = puUncertainty(ttjets_nlo+ttg, lowEmht, dirDir, nBins)

        zgScaleUnc = scaleUncertainty(zg+znunu, lowEmht, dirDir, nBins)
        wgScaleUnc = scaleUncertainty(wg+wjets, lowEmht, dirDir, nBins)
        tgScaleUnc = scaleUncertainty(ttjets_nlo+ttg, lowEmht, dirDir, nBins)

        zgJesUnc = jecUncertainty(zg+znunu, lowEmht, dirDir, nBins)
        wgJesUnc = jecUncertainty(wg+wjets, lowEmht, dirDir, nBins)
        tgJesUnc = jecUncertainty(ttjets_nlo+ttg, lowEmht, dirDir, nBins)

        mcSystUncert = 0.04 # SF, lumi, trigger
        zgSyst = aux.getSysHisto(zgHist, mcSystUncert)
        wgSyst = aux.getSysHisto(wgHist, mcSystUncert)
        tgSyst = aux.getSysHisto(tgHist, mcSystUncert)

        zgSyst = aux.addUncertaintiesQuadratic([zgSyst,zgPdfUnc,zgScaleUnc,zgJesUnc,zgPuUnc])
        wgSyst = aux.addUncertaintiesQuadratic([wgSyst,wgPdfUnc,wgScaleUnc,wgJesUnc,wgPuUnc])
        tgSyst = aux.addUncertaintiesQuadratic([tgSyst,tgPdfUnc,tgScaleUnc,tgJesUnc,tgPuUnc])

        totStat = aux.addHists(gjetHist, eHist, zgHist, wgHist, tgHist)
        totSyst = aux.addHists(gjetSyst, eSyst, zgSyst, wgSyst, tgSyst)

        signal1 = metHist(t5wg, "1600_100/{}/nominal".format(dirDir), nBins, lowEmht)
        signal2 = metHist(t5wg, "1750_1650/{}/nominal".format(dirDir), nBins, lowEmht)
        for h in signal1, signal2:
            aux.drawOpt(h, "signal")
        #    h.Add(totStat)
        signal1.SetLineColor(ROOT.kMagenta+2)
        signal2.SetLineColor(ROOT.kMagenta)
        signal2.SetLineStyle(2)

        #signal1_pre = aux.createHistoFromDatasetTree(t5wg_1600_100, "met*{}".format(info["shift"]), weight, nBins, "tr_jControl/simpleTree")
        #signal1_pre.Scale(info["scale"])

    if "electronClosure" in name:
        totStat = eHist
        totSyst = eSyst
    if "qcdClosure" in name:
        totStat = gjetHist
        totSyst = gjetSyst

    totUnc = aux.addHistUncert(totStat, totSyst)
    aux.drawOpt(totUnc, "totUnc")

    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    if dirSet == data:
        if "dPhiJet_lowEMHT" in name or "dPhiJet_highEMHT" in name:
            for b in range(dirHist.GetNbinsX()-3, dirHist.GetNbinsX()+2):
                dirHist.SetBinContent(b,0)
        m.add(dirHist, "Data")
    else:
        m.add(dirHist, "Direct simulation")
    if "Closure" not in name:
#        m.add(signal1_pre, "contamination")
        m.add(signal1, "T5Wg 1600 100")
        m.add(signal2, "T6gg 1750 1650")
        #m.add(signal2, "TChiWG 700")
        m.addStack(eHist, "e#rightarrow#gamma")
        m.addStack(zgHist, "#gammaZ")
        m.addStack(tgHist, "#gammat#bar{t}")
        m.addStack(wgHist, "#gammaW")
        m.addStack(gjetHist, "Nongenuine #it{p}_{T}^{miss}")
    if "electronClosure" in name:
        m.addStack(eHist, "e#rightarrow#gamma prediction")
    if "qcdClosure" in name:
        m.addStack(gjetHist, "Nongenuine #it{p}_{T}^{miss} prediction")

    m.add(totUnc, "Total uncertainty")
    m.maximum = 2.6*m.getMaximum()
    m.minimum = m.getMinimum()
    if "qcdClosure_lowEMHT" in name: m.minimum = 8e-3
    if "qcdClosure_highEMHT" in name: m.minimum = 6e-3
    if "ee_lowEMHT" in name: m.minimum = 1e-2
    if "ee_highEMHT" in name: m.maximum = 4.0*m.getMaximum()
    if "ee_highEMHT" in name: m.minimum = 5e-4
    if "final_lowEMHT" in name: m.minimum = 4e-2
    if "final_highEMHT" in name: m.minimum = 2e-3
    if "allMC_lowEMHT" in name: m.minimum = 7e-3
    if "allMC_highEMHT" in name: m.minimum = 8e-4
    if "dPhiJet_lowEMHT" in name: m.minimum = 1e-2
    if "dPhiJet_highEMHT" in name: m.minimum = 2e-4
    legInfo = "#it{H}_{T}^{#gamma} < 2TeV" if "lowEMHT" in name else "2TeV < #it{H}_{T}^{#gamma}"
    if "ee" in name: legInfo += ", EE"
    legInfo += ", |#Delta#phi| #lower[-.1]{>} 0.3"
    m.leg.SetHeader(legInfo)
    if "Closure" in name:
        m.leg.SetX1(.51)
        m.leg.SetX2(.93)
    else:
        m.leg.SetY1(.56)
        m.leg.SetX1(.56)
        m.leg.SetX2(.99)


    m.Draw()

    # draw other labels
    if "electronClosure" not in name:
        l = ROOT.TLine()
        l.SetLineStyle(2)
        l.SetLineColor(ROOT.kGray+2)
        text = ROOT.TLatex()
        text.SetTextSize(0.8*text.GetTextSize())
        l.DrawLine(100, 0, 100, totUnc.GetBinContent(totUnc.FindBin(100)))
        text.SetTextAngle(90)
        text.DrawLatexNDC(.23,.315, "Normalization")
        text.SetTextAngle(0)
        text.DrawLatexNDC(.311,.315, "Validation")
        if "final" in name:
            l.DrawLine(350, 0, 350, totUnc.GetBinContent(totUnc.FindBin(350)))
#            if "lowEMHT" in name:  text.DrawLatexNDC(.58,.30, "#font[62]{#color[1]{Search regions}}")
#            if "highEMHT" in name: text.DrawLatexNDC(.58,.31, "#font[62]{#color[1]{Search regions}}")
            #text.DrawLatexNDC(.0,.0, "#color[0]{We are one planet}")

    if "Closure" in name:
        r = ratio.Ratio("Ratio  ", dirHist, totStat, totSyst)
        rMax = 2
        if name=="electronClosure_lowEMHT": rMax = 2.7
        r.draw(0., rMax, None, True)
    else:
        r = ratio.Ratio("#scale[.9]{#lower[.24]{#splitline{Data/Pred.}{Bkg. frac.}}}", dirHist, totStat, totSyst)
        rMax = 1.5
        if name == "ee_highEMHT": rMax = 1.8
        if name == "final_lowEMHT": rMax = 1.6
        if name == "final_highEMHT": rMax = 3.6
        r.draw(0., rMax, m.getStack(), True)

    aux.Label(sim= not dirSet==data, status="" if "allMC" not in name else "Private Work")
    aux.save(name, normal=False, changeMinMax=False)

    if name == "final_lowEMHT": dc = limitTools.MyDatacard()
    elif name == "final_highEMHT": dc = limitTools.MyDatacard("testDatacard.txt")
    else: return
    for bin in range(dirHist.GetNbinsX()-2, dirHist.GetNbinsX()+1):
        binName = "bin{}_{}".format(name.split("_")[1],bin)
        bw = dirHist.GetBinWidth(bin) if style.divideByBinWidth else 1
        dc.addBin(binName, int(round(dirHist.GetBinContent(bin)*bw)),
            {
                "signal": (signal1.GetBinContent(bin)-totStat.GetBinContent(bin))*bw,
                #"signal": (signal1.GetBinContent(bin)-totStat.GetBinContent(bin)-signal1_pre.GetBinContent(bin))*bw,
                "gqcd": gjetHist.GetBinContent(bin)*bw,
                "ele": eHist.GetBinContent(bin)*bw,
                "wg": wgHist.GetBinContent(bin)*bw,
                "zg": zgHist.GetBinContent(bin)*bw,
                "tg": tgHist.GetBinContent(bin)*bw,
#                "cont": signal1_pre.GetBinContent(bin)*bw
            }, {
                "gqcdStat_"+binName: {"gqcd": getDatacardUncertFromHist(gjetHist,bin)},
                "eleStat_"+binName: {"ele": getDatacardUncertFromHist(eHist,bin)},
                "wgStat_"+binName: {"wg": getDatacardUncertFromHist(wgHist,bin)},
                "zgStat_"+binName: {"zg": getDatacardUncertFromHist(zgHist,bin)},
                "tgStat_"+binName: {"tg": getDatacardUncertFromHist(tgHist,bin)},
                "signalStat_"+binName: {"signal": getDatacardUncertFromHist(signal1,bin)},
                "gqcdSyst": {"gqcd": getDatacardUncertFromHist(gjetSyst,bin)},
                "eleSyst": {"ele": getDatacardUncertFromHist(eSyst,bin)},
                "pdf": {
                    "wg": getDatacardUncertFromHist(wgPdfUnc,bin),
                    "zg": getDatacardUncertFromHist(zgPdfUnc,bin),
                    "tg": getDatacardUncertFromHist(tgPdfUnc,bin)},
                "scale": {
                    "wg": getDatacardUncertFromHist(wgScaleUnc,bin),
                    "zg": getDatacardUncertFromHist(zgScaleUnc,bin),
                    "tg": getDatacardUncertFromHist(tgScaleUnc,bin)},
                "lumi": {
                    "signal": 1.026,
                    "wg": 1.026,
                    "zg": 1.026,
                    "tg": 1.026},
                "pu": {
                    "wg": getDatacardUncertFromHist(wgPuUnc,bin),
                    "zg": getDatacardUncertFromHist(zgPuUnc,bin),
                    "tg": getDatacardUncertFromHist(tgPuUnc,bin)},
                "jes": {
                    "wg": getDatacardUncertFromHist(wgJesUnc,bin),
                    "zg": getDatacardUncertFromHist(zgJesUnc,bin),
                    "tg": getDatacardUncertFromHist(tgJesUnc,bin)},
                "dataMC": {
                    "signal": 1.05,
                    "wg": 1.025,
                    "zg": 1.025,
                    "tg": 1.025},
                "trigger": {
                    "signal": 1.004,
                    "wg": 1.004,
                    "zg": 1.004,
                    "tg": 1.004},
                "isr": {"signal": 1.001},
                "genMet": {"signal": 1.001},
                # jes, jer splitting
            }
        )
    dc.write("testDatacard.txt")

if __name__ == "__main__":
    ewk_highestHT = wjets1200+wjets2500+ttjets600+ttjets800+ttjets1200+ttjets2500
    gqcd_highestHT = gjets600dr+qcd2000
    gqcd_highestHT.label = "(#gamma)+jet"
    #finalDistributionSignalHist("qcdClosure", True, gqcd, "original", gqcd)
    #finalDistributionSignalHist("qcdClosure", False, gqcd_highestHT, "original", gqcd)
    #finalDistributionSignalHist("electronClosure", True, ttjets_ht+wjets, "original")
    #finalDistributionSignalHist("electronClosure", False, ewk_highestHT, "original")
    #finalDistributionSignalHist("ee", True, data, "original_ee", dataHt)
    #finalDistributionSignalHist("ee", False, data, "original_ee", dataHt)
    finalDistributionSignalHist("final", True, data, "original", dataHt)
    finalDistributionSignalHist("final", False, data, "original", dataHt)

