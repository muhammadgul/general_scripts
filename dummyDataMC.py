#!/usr/bin/env python

"""Plots a data/MC plot to visualize colour coding of simulation.

Most of templates are taken from simulation.
"""
"""
source root6/build/bin/thisroot.sh
export PATH="/home/muhammad/root6/build/bin:$PATH"
export LD_LIBRARY_PATH="/home/muhammad/root6/build/lib:$LD_LIBRARY_PATH"
export PATH="/home/muhammad/anaconda3/bin:$PATH"

run as: python3.5 dummyDataMC.py
"""
path = '/home/muhammad/work/sahil/scripts/root_files/'
from contextlib import ExitStack
import glob
import math
import os
import numpy as np
from matplotlib import pyplot as plt

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
    # xsec order:    {ttbar, st_tch, st_tw, wjets,    DY,    WW,    WZ,    ZZ,   ,ttW,    ttZ,    qcd pt 15-20,       qcd pt 20-30,      qcd pt 30-50,       qcd pt 50-80,        qcd pt 80-120, qcd pt 120-170, qcd pt 170-300, qcd pt 300-470, qcd pt 470-600,     qcd pt 600-800, qcd pt 800-1000, qcd pt 1000-Inf
    xsec_dic = {'ttbar': 831.76, 'st_tch': 70.70, 'st_twch_t': 19.55, 'st_twch_at':35.6, 'st_sch':11.36*0.32, 'wjets': 61526.7, 'DY': 5765.4, 'WW':  118.7, 'WZ': 47.13 ,'ZZ': 16.523, 'ttW': 0.2043, 'ttZ': 0.2529,'qcd_pt_15_20': 1273190000.*0.003,'qcd_pt_20_30':  558528000.*0.0053, 'qcd_pt_30_50': 139803000.*0.01182, 'qcd_pt_50_80': 19222500.*0.02276,'qcd_pt_80_120':  2758420.*0.03844,'qcd_pt_120_170': 469797.*0.05362,'qcd_pt_170_300': 117989.*0.07335, 'qcd_pt_300_470': 7820.25*0.10196,'qcd_pt_470_600': 645.528* 0.12242, 'qcd_pt_600_800':187.109*0.13412,'qcd_pt_800_1000': 32.3486*0.14552,'qcd_pt_1000_Inf': 10.4305*0.15544}
     #qcd pt20-30,        qcd pt 30-50,     , qcd pt 50-80,       qcd pt 80-120,  qcd pt 120-170,   qcd pt 170-300, qcd pt 300-Inf
    emqcd_xsec = [558528000.*0.0053, 139803000.*0.01182, 19222500.*0.02276, 2758420.*0.03844, 469797.*0.05362, 117989.*0.07335, 7820.25*0.10196]

    evnt_dic = {'ttbar': 75676414.188449, 'st_tch':1220447540.405972 , 'st_twch_t':5477541.757331 , 'st_twch_at':5494360.755072, 'st_sch':5477541.757331, 'wjets': 6814531160594.228516, 'DY': 1662002371648.298828, 'WW': 6660781.662860 , 'WZ':1016101.679828  ,'ZZ':1005798.714968 , 'ttW': 741220.330057, 'ttZ':508217.97641 ,'qcd_pt_15_20':4693705.289944 ,'qcd_pt_20_30':32032749.277747, 'qcd_pt_30_50':28809945.038684 , 'qcd_pt_50_80':20720147.969795 ,'qcd_pt_80_120': 14085432.511000 ,'qcd_pt_120_170': 8056061.989751,'qcd_pt_170_300':7918425.238714 , 'qcd_pt_300_470':7609923.745327 ,'qcd_pt_470_600': 3926058.285534, 'qcd_pt_600_800':4027549.343103,'qcd_pt_800_1000': 4019175.159715,'qcd_pt_1000_Inf': 3947838.760831}
    lumi = 12877.4 #muon = 12877.401701, e = 12883.8601
    Hists  = {}
    simHists = {}
    dataHists = {}
    hist_basket = []
    dataHistnames= []
    simHistnames = []

    def build_hist(name, title):
        return ROOT.TH1D(name, title, 45, 300., 1200.)
    for filenames in glob.glob(os.path.join(path, '*.root')):
      scale = 1.0
      file_ = ROOT.TFile.Open(filenames, 'read')
      hist_name = file_.Get('mujets_2btag/Pt_lep')
      hist_bins = hist_name.GetNbinsX()
      hist_basket_temp = []
      hist_basket_tuple_temp1 = []
      histname = []

      for ibin in range(hist_bins):
        bin_cont = hist_name.GetBinContent(ibin)
#        print('this is hist bins',bin_cont)
        hist_basket_temp.append(bin_cont)
#        print('this is cont: ',hist_basket)
#      canvas_test = ROOT.TCanvas('canvas', '', 500, 500)
#      plt.figure()
      if 'MC13TeV_TTJets' in filenames:
        histname = 'MC13TeV_TTJets'
        scale = xsec_dic['ttbar']*lumi/evnt_dic['ttbar']
      elif 'TTWToLNu' in filenames:
        histname = 'TTWToLNu'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'MC13TeV_TTZToLLNuNu' in filenames:
        histname = 'MC13TeV_TTZToLLNuNu'
        scale = 0.0
      elif 'MC13TeV_TTZToQQ' in filenames:
        histname = 'MC13TeV_TTZToQQ'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'MC13TeV_TTWToQQ' in filenames:
        histname = 'MC13TeV_TTWToQQ'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']


      elif 'MC13TeV_SingleT_t_' in filenames:
        histname = 'MC13TeV_SingleT_t'
        scale = xsec_dic['st_tch']*lumi/evnt_dic['st_tch']
      elif 'MC13TeV_ST_sch' in filenames:
        histname = 'MC13TeV_ST_sch'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'MC13TeV_SingleTbar_tW' in filenames:
        histname = 'MC13TeV_SingleTbar_tW'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'MC13TeV_SingleT_tW' in filenames:
        histname = 'MC13TeV_SingleT_tW'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'SingleTbar_t_' in filenames:
        histname = 'SingleTbar_t_'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'W1Jet' in filenames:
        histname = 'W1Jet'
        scale = 0.0
      elif 'W2Jet' in filenames:
        histname = 'W2Jet'
        scale = 0.0
      elif 'W3Jet' in filenames:
        histname = 'W3Jet'
        scale = 0.0
      elif 'W4Jet' in filenames:
        histname = 'W4Jet'
        scale = 0.0
      elif 'MC13TeV_ZZ' in filenames:
        histname = 'MC13TeV_ZZ'
        scale = xsec_dic['ZZ']*lumi/evnt_dic['ZZ']
      elif 'MC13TeV_WW' in filenames:
        histname = 'MC13TeV_WW'
        scale = xsec_dic['WW']*lumi/evnt_dic['WW']
      elif 'MC13TeV_WZ' in filenames:
        histname = 'MC13TeV_WZ'
        scale = xsec_dic['WZ']*lumi/evnt_dic['WZ']
      elif 'MC13TeV_WJets' in filenames:
        histname = 'MC13TeV_WJets'
        scale = xsec_dic['wjets']*lumi/evnt_dic['wjets']

      elif 'MC13TeV_DY50toInf_mlm' in filenames:
        histname = 'MC13TeV_DY50toInf_mlm'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'QCD_Pt-15to20_MuEnriched' in filenames:
        histname = 'QCD_Pt-15to20_MuEnriched'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'MC13TeV_QCD_Pt-20to30_MuEnriched' in filenames:
        histname = 'MC13TeV_QCD_Pt-20to30_MuEnriched'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'MC13TeV_QCD_Pt-30to50_MuEnriched' in filenames:
        histname = 'MC13TeV_QCD_Pt-30to50_MuEnriched'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'QCD_Pt-50to80_MuEnriched' in filenames:
        histname = 'QCD_Pt-50to80_MuEnriched'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'MC13TeV_QCD_Pt-80to120_MuEnriched' in filenames:
        histname = 'MC13TeV_QCD_Pt-80to120_MuEnriched'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'MC13TeV_QCD_Pt-120to170_MuEnriched' in filenames:
        histname = 'MC13TeV_QCD_Pt-120to170_MuEnriched'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'MC13TeV_QCD_Pt-170to300_MuEnriched' in filenames:
        histname = 'MC13TeV_QCD_Pt-170to300_MuEnriched'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'MC13TeV_QCD_Pt-300to470_MuEnriched' in filenames:
        histname = 'MC13TeV_QCD_Pt-300to470_MuEnriched'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'MC13TeV_QCD_Pt-470to600_MuEnriched' in filenames:
        histname = 'MC13TeV_QCD_Pt-470to600_MuEnriched'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'QCD_Pt-600to800_MuEnriched' in filenames:
        histname = 'QCD_Pt-600to800_MuEnriched'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'MC13TeV_QCD_Pt-800to1000_MuEnriched' in filenames:
        histname = 'MC13TeV_QCD_Pt-800to1000_MuEnriched'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'MC13TeV_QCD_Pt-1000toInf_MuEnriched' in filenames:
        histname = 'MC13TeV_QCD_Pt-1000toInf_MuEnriched'
        scale = xsec_dic['ttbar']*lumi*0.0/evnt_dic['ttbar']
      elif 'Data13TeV_SingleMuon_2016B_bk' in filenames:
        histname = 'Data13TeV_SingleMuon_2016B_bk'
      elif 'Data13TeV_SingleMuon_2016C' in filenames:
        histname = 'Data13TeV_SingleMuon_2016C'
      elif 'Data13TeV_SingleMuon_2016D' in filenames:
        histname = 'Data13TeV_SingleMuon_2016D'
      elif 'Data13TeV_SingleMuon_2016E' in filenames:
        histname = 'Data13TeV_SingleMuon_2016E'
      elif 'Data13TeV_SingleMuon_2016F' in filenames:
        histname = 'Data13TeV_SingleMuon_2016F'
      elif 'Data13TeV_SingleMuon_2016G' in filenames:
        histname = 'Data13TeV_SingleMuon_2016G'
      elif 'Data13TeV_SingleMuon_2016H_v2' in filenames:
        histname = 'Data13TeV_SingleMuon_2016H_v2'
      elif 'Data13TeV_SingleMuon_2016H_v3' in filenames:
        histname = 'Data13TeV_SingleMuon_2016H_v3'
      if "Data" not in filenames:
        simHistnames.append(histname)
      if "Data" in filenames:
        dataHistnames.append(histname)

      hist_basket_tuple_temp1.append(histname)
      hist_basket_tuple_temp1.append(histname)
      hist_basket_temp = [x * scale for x in hist_basket_temp]
      hist_basket_tuple_temp1.append(hist_basket_temp)
      hist_basket.append(tuple(hist_basket_tuple_temp1))
#      print("this is is file: ",filenames)
#      print('this is scale: ',scale)
    print('dataHistnames: ',simHistnames)
    for histName, histTitle, binContent in hist_basket:
#      print('this is histName: ',histName)
      hist = build_hist(histName, histTitle)
      for bin, content in enumerate(binContent):
        hist.SetBinContent(bin + 1, content)
      Hists[hist.GetName()] = hist
#    del hist_basket_tuple_temp1[:]
#    del hist_basket_temp[:]
#    del hist_basket[:]
    dataHist = build_hist('data', 'Data')
    for name in dataHistnames:
      dataHist.Add(Hists[name])
    Hists[dataHist.GetName()] = dataHist


    # Predefined elements of the decoration
    dataHist.SetMarkerStyle(ROOT.kFullCircle)
    dataHist.SetMarkerColor(ROOT.kBlack)
    dataHist.SetLineColor(ROOT.kBlack)

    histStack = ROOT.THStack('histStack', ';m_{t#bar{t}} [GeV];Events')
    for name in simHistnames:
      histStack.Add(Hists[name], 'hist')
#    histStack.Add(Hists['MC13TeV_TTJets'], 'hist')
      histStack.Draw()

    canvas = ROOT.TCanvas('canvas', '', 1500, 1000)
    canvas.SetRightMargin(0.15)
    canvas.SetTicks()
    histStack.Draw()
    dataHist.Draw('p0 e1 x0 same')

    canvas.Print('dummyDataMC.png')

"""
    print('this is datahistnames',dataHistnames)







    for histName, histTitle, binContent in [
        ('tt', 't#bar{t}', [307.157, 742.939, 1288.87, 1608.22, 1681.45, 1608.72, 1486.5, 1335.82, 1198.34, 1055.04, 929.448, 820.387, 721.512, 643.229, 564.359, 500.1, 439.005, 389.639, 340.815, 304.501, 270.951, 239.384, 210.131, 187.127, 167.341, 148.691, 131.85, 118.003, 104.475, 94.6695, 81.6937, 74.4348, 69.1456, 61.9767, 55.2322, 46.1541, 43.6128, 39.5578, 34.2585, 31.423, 27.5023, 26.3051, 23.3278, 19.7781, 19.2375]),
        ('t', 't', [17.8125, 31.2617, 53.9018, 64.1869, 65.5968, 64.3125, 60.4394, 55.3419, 49.3731, 45.768, 45.5789, 38.2493, 35.3382, 31.5194, 25.0101, 23.623, 22.4006, 18.8849, 17.3471, 14.6294, 13.2772, 13.4174, 9.8228, 10.8933, 8.35154, 8.61042, 6.53047, 7.20063, 6.36081, 5.19556, 5.24923, 4.26612, 4.2257, 4.21102, 3.81628, 4.01056, 2.92168, 2.58807, 1.87836, 2.97243, 1.71062, 1.46511, 1.47121, 1.13739, 1.12863]),
        ('W', 'W', [8.19584, 15.9918, 21.8019, 47.7245, 22.7086, 26.0647, 38.1045, 21.449, 18.5422, 28.5823, 34.6552, 22.1832, -0.172199, 3.17567, 11.2051, 11.6841, 2.50422, 9.26819, 10.568, 3.11688, 8.6025, 16.4938, 8.03359, 2.68415, 0.464071, 1.17749, 6.75812, 11.5448, -3.05647, 6.95907, 0.25755, 0.110415, 7.20289, 4.48571, 1.25671, 2.79944, -1.39284, 1.71902, -0.52002, 1.0918, -1.10865, 3.33987, 0.96795, 1.6982, 4.5997]),
        ('Z', 'Z/gamma*', [1.8229, 5.56976, 0.57534, 1.51224, 3.04875, 7.09425, -0.675345, 5.03736, 3.97793, 0.721789, -0.583307, -1.28495, 2.17486, 2.75874, 2.4544, 0.994734, 0.438286, 0.470855, 1.2729, 1.15163, 0.0254503, -0.256914, 1.44765, 0, 0.848283, -1.30131, 0.0417647, 2.20941, -0.238122, 0.745578, 1.67972, 1.46361, 0.899132, 0.0409446, 0.228341, 0.632722, 0, -0.714332, 0, 0, 0, -0.733536, -0.711733, 0, 0])
    ]:
        hist = build_hist(histName, histTitle)
        
        for bin, content in enumerate(binContent):
            hist.SetBinContent(bin + 1, content)
        
        simHists[hist.GetName()] = hist
    
    hist = build_hist('QCD', 'QCD')
    for bin in range(1, hist.GetNbinsX() + 1):
        hist.SetBinContent(bin, 20. * math.exp(-bin / 20.))
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
    
    
    # Merge some histograms into groups
    hist = build_hist('EW', 'EW')
    hist.Add(simHists['W'])
    hist.Add(simHists['Z'])
    # Also VV and ttV would go here, but I don't create histograms for
    # them
    simHists[hist.GetName()] = hist
    
    
    # Predefined elements of the decoration
    dataHist.SetMarkerStyle(ROOT.kFullCircle)
    dataHist.SetMarkerColor(ROOT.kBlack)
    dataHist.SetLineColor(ROOT.kBlack)
    
    for hist in simHists.values():
        hist.SetLineWidth(1)
        hist.SetLineColor(ROOT.kBlack)
    
    
    # Custom decoration
    simHists['tt'].SetFillColor(ROOT.kOrange + 1)
    simHists['t'].SetFillColor(ROOT.kMagenta)
    simHists['EW'].SetFillColor(ROOT.kGreen + 1)
    simHists['QCD'].SetFillColor(ROOT.kGray)
    
    
    # Define what background histograms and in which order should be
    # drawn
    simHistOrder = ['tt', 't', 'EW', 'QCD']
    
    
    # Draw the histograms
    canvas = ROOT.TCanvas('canvas', '', 1500, 1000)
    canvas.SetRightMargin(0.15)
    canvas.SetTicks()
    
    histStack = ROOT.THStack('histStack', ';m_{t#bar{t}} [GeV];Events')
    for name in reversed(simHistOrder):
        histStack.Add(simHists[name], 'hist')
    histStack.Draw()
    
    dataHist.Draw('p0 e1 x0 same')
    
    histMax = 1.1 * max(histStack.GetMaximum(), dataHist.GetMaximum())
    histStack.SetMaximum(histMax)
    dataHist.SetMaximum(histMax)
    
    
    # Draw the legend
    legend = ROOT.TLegend(0.86, 0.9 - 0.04 * (len(simHistOrder) + 1), 0.99, 0.9)
    legend.SetFillColor(ROOT.kWhite)
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
    
    energyLabel = ROOT.TLatex(0.85, 0.91, '30 fb^{-1} (13 TeV)')
    energyLabel.SetNDC()
    energyLabel.SetTextFont(42)
    energyLabel.SetTextSize(0.04)
    energyLabel.SetTextAlign(31)
    energyLabel.Draw()
    
    
    # Everything is done
    canvas.Print('dummyDataMC.pdf')
    canvas.Print('dummyDataMC.png')
    """
