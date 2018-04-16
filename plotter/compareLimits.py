#!/usr/bin/env python2

from include import *

def compareStrong():
    fname = "limitCalculations/SMS-T5Wg_{}/saved_graphs1d_limit.root"

    hist = aux.getFromFile(fname.format("original_test"), "obs_hist")
    #hist.SetTitle("")
    hist.GetXaxis().SetRangeUser(1400, 2800)
    hist.Reset()

    c = ROOT.TCanvas()
    hist.Draw("col")

    scenarios = collections.OrderedDict([
        ("original_test", "original"),
        ("di_cleaned_test", "#gamma#gamma veto"),
        ("lep_cleaned_test", "#gammal veto"),
        ("all_cleaned_test", "#gamma#gamma, #gammal, S^{#gamma}_{T} veto"),
        ("dilep_cleaned_test_highHtg", "#gamma#gamma, #gammal veto, H_{T}^{#gamma}>2TeV")])

    leg = ROOT.TLegend(.56,.49,.94,.915)
    leg.SetFillColor( ROOT.kWhite )
    leg.SetFillStyle(0)

    colors = [ rwth.black, rwth.red, rwth.lila, rwth.petrol50, rwth.green ]+range(100)

    for iN, (n, lab) in enumerate(scenarios.iteritems()):
        gr = aux.getFromFile(fname.format(n), "exp")
        gr.SetLineColor(colors[iN])
        gr.SetLineWidth(2)
        gr.Draw("same")
        aux.saveStuff.append(gr)
        leg.AddEntry(gr, lab, "l")
    leg.Draw()
    aux.save("limits_compare_strong", log=False)

def compareWeak():
    fname = "limitCalculations/SMS-TChiNG_BF50N50G_{}/Graphs2d.root"

    hist = ROOT.TH1F("", "TChiNG;m_{#chi} (GeV); cross section upper limit (pb)", 10, 200, 1400)

    c = ROOT.TCanvas()
    hist.Draw()

    scenarios = collections.OrderedDict([
        ("original_test", "original"),
        ("dilep_cleaned_test", "#gamma#gamma, #gammal veto"),
        ("all_cleaned_test", "#gamma#gamma, #gammal, S^{#gamma}_{T} veto"),
        ("dilep_cleaned_test_highHtg", "#gamma#gamma, #gammal veto, H_{T}^{#gamma}>2TeV")])

    leg = ROOT.TLegend(.56,.49,.94,.915)
    leg.SetFillColor( ROOT.kWhite )
    leg.SetFillStyle(0)

    colors = [ rwth.black, rwth.red, rwth.lila, rwth.petrol50, rwth.green ]+range(100)

    for iN, (n, lab) in enumerate(scenarios.iteritems()):
        gr = aux.getFromFile(fname.format(n), "exp")
        gr.SetLineColor(colors[iN])
        gr.SetLineWidth(2)
        gr.Draw("same")
        aux.saveStuff.append(gr)
        leg.AddEntry(gr, lab, "l")
    leg.Draw()
    aux.save("limits_compare_weak", log=False)

if __name__ == "__main__":
    compareStrong()
    #compareWeak()
