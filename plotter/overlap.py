#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *
import array

def getSr(met, htg):
    if htg < 2000:
        if met > 600: return 3
        elif met > 450: return 2
        else: return 1
    else:
        if met > 600: return 6
        elif met > 450: return 5
        else: return 4

def ente():
    with open(inName) as f:
        lines = f.readlines()

    texTable = "Run & RunNr & LumNr & EventNr & ptmiss & htg & pt & \\\\\n"

    srBins = [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    for i in range(0,len(lines)-1, 5):
        e = Event()
        e.evtId = lines[i].split(" ")[0]
        e.area = lines[i].split(" ")[-1].strip()
        e.met = float(lines[i+1].split(" ")[2])
        e.htg = float(lines[i+1].split(" ")[5])
        pt = float(lines[i+1].split(" ")[9])
        run, lum, nummer = e.evtId.split(":")
        texTable += " {} & {} & {} & {} & {} & {} & {} &".format(e.area, run, lum, nummer, e.met, e.htg, pt)
        e.gam = " 1 1" in lines[i+2]
        e.lep = " 1 1" in lines[i+3]
        if e.gam: texTable += "Selected by $\\gamma+\\gamma$"
        if e.lep: texTable += "Selected by $\\ell+\\gamma$"
        texTable += "\\\\\n"
        e.sr = getSr(e.met, e.htg)
        if e.gam: srBins[e.sr-1][2] += 1
        elif e.lep: srBins[e.sr-1][1] += 1
        else: srBins[e.sr-1][0] += 1
    print srBins


    c = ROOT.TCanvas()
    c.Divide(3,2)

    saves = []
    for i in range(6):
        c.cd(i+1)
        vals = array.array("d", srBins[i])
        cols = array.array("i", [ROOT.kGray, ROOT.kRed, ROOT.kBlue])
        pie = ROOT.TPie("pie%s"%i, "", 3, vals, cols)
        pie.Draw("nol")
        pie.SetLabelFormat("%val")
        saves.append(pie)
    c.cd(0)
    c.SaveAs("test.pdf")

    print texTable


def parseInputFile(inName):
    with open(inName) as f:
        lines = f.read()
    pattern = re.compile("(?P<run>\d+):(?P<lum>\d+):(?P<evt>\d+) in Run (?P<area>.)\nptmiss = (?P<ptmiss>[\d\.]+)\thtg = (?P<htg>[\d\.]+)\tphoton pt = (?P<pt>[\d\.]+)\nselected by gam\+gam: (?P<gg_1>.) (?P<gg_2>.)\nselected by gam\+lep: (?P<lg_1>.) (?P<lg_2>.)\nselected by gam\+stg: (?P<st_1>.) (?P<st_2>.)\nselected by gam\+htg: (?P<ht_1>.) (?P<ht_2>.)")
    out = [entry.groupdict() for entry in pattern.finditer(lines)]
    if len(out) != len(re.compile("in Run").findall(lines)): print "WARNING: Maybe not all entries parsed correctly"
    return out

def createPlots(infos):
    hDefault = ROOT.TH1F("", ";signal regions;", 6, 0, 6)
    searchRegionLabels = "350-450", "450-600", "600-#infty  "
    for b in range(6):
        hDefault.GetXaxis().SetBinLabel(b+1, searchRegionLabels[b%3])
    hDefault.SetXTitle("#it{p}_{T}^{miss} in #it{H}_{T}^{#gamma}>2TeV     #it{p}_{T}^{miss} in #it{H}_{T}^{#gamma}<2TeV")

    hs = dict((x,hDefault.Clone(x)) for x in ["ht", "st", "lg", "gg"])
    for info in infos:
        sr = getSr(float(info["ptmiss"]), float(info["htg"]))
        if info["gg_1"] == "1": hs["gg"].AddBinContent(sr)
        elif info["lg_1"] == "1": hs["lg"].AddBinContent(sr)
        elif info["st_1"] == "1": hs["st"].AddBinContent(sr)
        else: hs["ht"].AddBinContent(sr)
    hs["ht"].SetLineColor(ROOT.kCyan)
    hs["st"].SetLineColor(ROOT.kBlue)
    hs["lg"].SetLineColor(ROOT.kOrange)
    hs["gg"].SetLineColor(ROOT.kGreen)
    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    m.leg.SetHeader("Overlap with:")
    m.addStack(hs["gg"], "#gamma#gamma")
    m.addStack(hs["lg"], "#gammal")
    m.addStack(hs["st"], "#gammaS_{T}^{#gamma}")
    m.addStack(hs["ht"], "#gammaH_{T}^{#gamma}")

    m.Draw()
    ratio.clearXaxisCurrentPad()
    ratio.createBottomPad()
    aux.drawContributions(m.getStack())
    aux.Label(sim=False)
    aux.save("overlapInSignalRegions")



def createTexFile(infos):
    texTable = "Run & RunNr & LumNr & EventNr & ptmiss & htg & pt & \\\\\n"
    for info in infos:
        print info
        texTable += " {area} & {run} & {lum} & {evt} & {ptmiss} & {htg} & {pt} & ".format(**info)
        if info["gg_1"] == "1": texTable += "Selected by $\\gamma+\\gamma$ "
        if info["lg_1"] == "1": texTable += "Selected by $\\gamma+\\ell$ "
        if info["st_1"] == "1": texTable += "Selected by $\\gamma+S_{T}^{\\gamma}$ "
        texTable += "\\\\\n"
    print texTable

def overlapCheck(inName):
    infos = parseInputFile(inName)
    createPlots(infos)
    #createTexFile(infos)

if __name__ == "__main__":
    inName = "../histogramProducer/overlap_v27.txt"
    overlapCheck(inName)

