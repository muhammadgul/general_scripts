#!/usr/bin/env python

"""Plots a data/MC plot to visualize colour coding of simulation.

All shapes are completely made up.

source root6/build/bin/thisroot.sh
export PATH="/home/muhammad/root6/build/bin:$PATH"
export LD_LIBRARY_PATH="/home/muhammad/root6/build/lib:$LD_LIBRARY_PATH"
python3.5 dummyDataMCTTReco.py
"""

import math

import ROOT
ROOT.gROOT.SetBatch(True)


if __name__ == '__main__':
    
    # Global style
    ROOT.gStyle.SetHistMinimumZero(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetStripDecimals(False)
    ROOT.TH1.SetDefaultSumw2(True)
    ROOT.TGaxis.SetMaxDigits(3)
    
    ROOT.gStyle.SetTitleFont(42)
    ROOT.gStyle.SetTitleFontSize(0.04)
    ROOT.gStyle.SetTitleFont(42, 'XYZ')
    ROOT.gStyle.SetTitleXOffset(0.9)
    ROOT.gStyle.SetTitleYOffset(1.1)
    ROOT.gStyle.SetTitleSize(0.045, 'XYZ')
    ROOT.gStyle.SetLabelFont(42, 'XYZ')
    ROOT.gStyle.SetLabelOffset(0.007, 'XYZ')
    ROOT.gStyle.SetLabelSize(0.04, 'XYZ')
    ROOT.gStyle.SetNdivisions(508, 'XYZ')
    
    
    # Dummy histograms for simulation.  Fill colour will be specified
    # later.  Histogram titles will be used as labels in the legend.
    simHists = {}
    def build_hist(name, title):
        return ROOT.TH1D(name, title, 45, 10., 25.)
    
    hist = build_hist('ttMatched', 't#bar{t}, matched')
    for bin in range(1, hist.GetNbinsX() + 1):
        x = hist.GetBinCenter(bin)
        hist.SetBinContent(bin, 10e3 * ROOT.TMath.Gaus(x, 13., 0.8))
    simHists[hist.GetName()] = hist
    
    hist = build_hist('ttMatchable', 't#bar{t}, matchable')
    for bin in range(1, hist.GetNbinsX() + 1):
        x = hist.GetBinCenter(bin)
        hist.SetBinContent(
            bin,
            0.5e3 * ROOT.TMath.Gaus(x, 13., 0.8) + 0.5e3 * ROOT.TMath.Gaus(x, 18., 5.)
        )
    simHists[hist.GetName()] = hist
    
    hist = build_hist('ttNonMatchable', 't#bar{t}, non-matchable')
    for bin in range(1, hist.GetNbinsX() + 1):
        x = hist.GetBinCenter(bin)
        hist.SetBinContent(
            bin,
            2e3 * ROOT.TMath.Gaus(x, 13., 1.) + 2e3 * ROOT.TMath.Gaus(x, 18., 5.)
        )
    simHists[hist.GetName()] = hist
    
    hist = build_hist('ttOther', 't#bar{t}, other decays')
    for bin in range(1, hist.GetNbinsX() + 1):
        x = hist.GetBinCenter(bin)
        hist.SetBinContent(
            bin,
            600. * ROOT.TMath.Gaus(x, 13., 1.) + 500. * ROOT.TMath.Gaus(x, 18., 5.)
        )
    simHists[hist.GetName()] = hist
    
    hist = build_hist('other', 'Other')
    for bin in range(1, hist.GetNbinsX() + 1):
        x = hist.GetBinCenter(bin)
        hist.SetBinContent(
            bin,
            500. * ROOT.TMath.Gaus(x, 13., 2.) + 500. * ROOT.TMath.Gaus(x, 18., 5.)
        )
    simHists[hist.GetName()] = hist
    
    
    # Create a fake data histogram
    dataHist = build_hist('data', 'Data')
    
    sumSimHist = build_hist('sumSim', '')
    for hist in simHists.values():
        sumSimHist.Add(hist)
    
    rGen = ROOT.TRandom3(1696)
    
    for bin in range(1, dataHist.GetNbinsX() + 1):
        content = rGen.Poisson(sumSimHist.GetBinContent(bin))
        dataHist.SetBinContent(bin, content)
        dataHist.SetBinError(bin, math.sqrt(content))
    
    del rGen
    del sumSimHist
    
    
    # Predefined elements of the decoration
    dataHist.SetMarkerStyle(ROOT.kFullCircle)
    dataHist.SetMarkerColor(ROOT.kBlack)
    dataHist.SetLineColor(ROOT.kBlack)
    
    for hist in simHists.values():
        hist.SetLineWidth(1)
        hist.SetLineColor(ROOT.kBlack)
    
    
    # Custom decoration
    simHists['ttMatched'].SetFillColor(ROOT.kOrange + 9)
    simHists['ttMatchable'].SetFillColor(ROOT.kRed - 1)
    simHists['ttNonMatchable'].SetFillColor(ROOT.kMagenta - 2)
    simHists['ttOther'].SetFillColor(ROOT.kCyan - 1)
    simHists['other'].SetFillColor(ROOT.kGray)
    
    
    # Define what background histograms and in which order should be
    # drawn
    simHistOrder = ['ttMatched', 'ttMatchable', 'ttNonMatchable', 'ttOther', 'other']
    
    
    # Draw the histograms
    canvas = ROOT.TCanvas('canvas', '', 1500, 1000)
    canvas.SetTicks()
    
    histStack = ROOT.THStack('histStack', ';#lambda_{comb};Events')
    for name in reversed(simHistOrder):
        histStack.Add(simHists[name], 'hist')
    histStack.Draw()
    
    dataHist.Draw('p0 e1 x0 same')
    
    histMax = 1.1 * max(histStack.GetMaximum(), dataHist.GetMaximum())
    histStack.SetMaximum(histMax)
    dataHist.SetMaximum(histMax)
    
    
    # Draw the legend
    legend = ROOT.TLegend(0.62, 0.87 - 0.04 * (len(simHistOrder) + 1), 0.87, 0.87)
    legend.SetFillColorAlpha(ROOT.kWhite, 0.)
    legend.SetTextFont(42)
    legend.SetTextSize(0.035)
    legend.SetBorderSize(0)
    
    legend.AddEntry(dataHist, dataHist.GetTitle(), 'p')
    for name in simHistOrder:
        legend.AddEntry(simHists[name], simHists[name].GetTitle(), 'f')
    
    legend.Draw()
    
    
    # Add some labels to the plot
    cmsLabel = ROOT.TLatex(0.17, 0.91, '#scale[1.2]{#font[62]{CMS}} #font[52]{Fake}')
    cmsLabel.SetNDC()
    cmsLabel.SetTextFont(42)
    cmsLabel.SetTextSize(0.04)
    cmsLabel.SetTextAlign(11)
    cmsLabel.Draw()
    
    energyLabel = ROOT.TLatex(0.90, 0.91, '30 fb^{-1} (13 TeV)')
    energyLabel.SetNDC()
    energyLabel.SetTextFont(42)
    energyLabel.SetTextSize(0.04)
    energyLabel.SetTextAlign(31)
    energyLabel.Draw()
    
    
    # Everything is done
    canvas.Print('dummyDataMCTTReco.pdf')
