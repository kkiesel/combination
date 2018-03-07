#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

from finalPredictions import metHist


def compareHists(h1, h2):
    for b in aux.loopH(h1):
        c1 = h1.GetBinContent(b)
        c2 = h2.GetBinContent(b)
        if abs(c1-c2) > 1e-8:
            print b, c1, c2

def validateHist():
    filePart = "SinglePhoton_Run2016H-03Feb2017_ver3-v1"
    filePart = "WGJets_MonoPhoton_PtG-130"
    filePart = "WJetsToLNu_HT-200To400"
    filePart = "SMS-T5Wg_1600_100"
    divNgen = "Run2016" not in filePart
    lowEmht = True

    dataset1 = Dataset(filePart, 0)
    file2 = "../../phd/histogramProducer/{}_hists.root".format(filePart)

    p1 = "0_0/original/nominal"
    p2 = "signal_lowEMHT/met"

    #p1 = "0_0/original/eControl"
    #p2 = "signal_lowEMHT_eControl/met"

    #p1 = "0_0/original/genE"
    #p2 = "signal_lowEMHT_genE/met"

    dataset1 = Dataset("SMS-T5Wg", "data/xSec_SMS_Gluino_13TeV.pkl")
    p1 = "1600_100/original/nominal"
    p2 = "signal_lowEMHT/met"



#    p1 = "0_0/original/jCR/0.860000"
#    p2 = "tr_jControl/simpleTree"



    if not lowEmht: p2 = p2.replace("low", "high")

    h1 = metHist(dataset1, p1, lowEmht=lowEmht)
    if "simpleTree" in p2:
        tree = ROOT.TChain(p2)
        tree.AddFile(file2)
        scale = float(p1.split("/")[-1])
        weight = "weight*(emht>700 && emht<2000)"
        if not lowEmht:
            weight = "weight*(emht>2000)"
        h2 = aux.createHistoFromTree(tree, "met*{}".format(scale), weight, 200, 0, 2000)
    else:
        h2 = aux.getFromFile(file2, p2)
    if divNgen:
        nGen = aux.getNgen(file2)
        h2.Scale(1./nGen)
    compareHists(h1, h2)



def validateSignalScan(scanName, point):

    xsec = 0.00810078 # for m_g = 1600
    weight = xsec * aux.intLumi

    binning = [350, 450, 600]
    # get signal yield old
    fOld = "../../phd/histogramProducer/{}_signalScan.root".format(scanName)
    hOldL = aux.getFromFile(fOld, "Wg_{}/signal_lowEMHT/met".format(point))
    hOldH = aux.getFromFile(fOld, "Wg_{}/signal_highEMHT/met".format(point))
    #hOldL.Add(aux.getFromFile(fOld, "Wg_{}/signal_lowEMHT/met_contamination".format(point)), -1)
    #hOldH.Add(aux.getFromFile(fOld, "Wg_{}/signal_highEMHT/met_contamination".format(point)), -1)

    hOldL = aux.rebin(hOldL, binning)
    hOldH = aux.rebin(hOldH, binning)
    acceptancesOld = [hOldL.GetBinContent(i+1) for i in range(3)] + [hOldH.GetBinContent(i+1) for i in range(3)]
    print [a*weight for a in acceptancesOld]

    dataset1 = Dataset(scanName, 0)
    hL = metHist(dataset1, point+"/original/nominal", binning+[800], True)
    hH = metHist(dataset1, point+"/original/nominal", binning+[800], False)
    acceptances = [hL.GetBinContent(i+1) for i in range(3)] + [hH.GetBinContent(i+1) for i in range(3)]
    print [a*weight for a in acceptances]
    print





if __name__ == "__main__":
    #validateHist()
    validateSignalScan("SMS-T5Wg", "1600_100")
    validateSignalScan("SMS-T6Wg", "1600_100")
