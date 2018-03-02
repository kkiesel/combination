#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *
import array

class Event:
    evtId = ""
    area = ""
    met = 0
    htg = 0
    sr = 0
    gam = False
    lep = False

def getSr(met, htg):
    if htg < 2000:
        if met > 600: return 3
        elif met > 450: return 2
        else: return 1
    else:
        if met > 600: return 6
        elif met > 450: return 5
        else: return 4

if __name__ == "__main__":
    inName = "../histogramProducer/overlap_v26.txt"

    with open(inName) as f:
        lines = f.readlines()

    texTable = "Run & RunNr & LumNr & EventNr & ptmiss & htg & pt & \\\\\n"

    events = []
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
        events.append(e)

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
