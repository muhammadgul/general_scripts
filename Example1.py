#!/usr/bin/env python

import sys
import time
import ROOT

try:
  input = raw_input
except:
  pass

if len(sys.argv) < 2:
  print(" Usage: Example1.py input_file")
  sys.exit(1)

ROOT.gSystem.Load("libDelphes")

try:
  ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
  ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
  pass
# Create a ROOT file to save histograms
file = ROOT.TFile("output_file.root", "RECREATE")

inputFile = sys.argv[1]

# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

# Get pointers to branches used in this analysis
branchWeight   = treeReader.UseBranch("Weight")
branchEvent    = treeReader.UseBranch("Event")
branchJet = treeReader.UseBranch("Jet")
branchElectron = treeReader.UseBranch("Electron")
branchMuon = treeReader.UseBranch("Muon")

# Book histograms
histJetsSize = ROOT.TH1F("jet_Size", "jet Size", 18, 0.0, 17.0)

histbJetsSize = ROOT.TH1F("bjets_Size", "bjets Size", 17, 0.0, 17.0)
bjets_hPt_l = []
bjets_hEta_l = []
bjets_hPhi_l = []
bjets_hMass_l = []
lep_hPt_l = []
lep_hEta_l = []
lep_hPhi_l = []
for i in range (4):
    bjets_Pt_hName = f"bjet{i}_hPt"
    hist_bjet_Pt = ROOT.TH1F(bjets_Pt_hName, bjets_Pt_hName, 100, 0, 400)
    bjets_hPt_l.append(hist_bjet_Pt)

    bjets_Eta_hName = f"bjet{i}_hEta"
    hist_bjet_Eta = ROOT.TH1F(bjets_Eta_hName, bjets_Eta_hName, 100, -4.0, 4.0)
    bjets_hEta_l.append(hist_bjet_Eta)

    bjets_Phi_hName = f"bjet{i}_hPhi"
    hist_bjet_Phi = ROOT.TH1F(bjets_Phi_hName, bjets_Phi_hName, 100, 0.0, 3.20)
    bjets_hPhi_l.append(hist_bjet_Phi)

    bjets_Mass_hName = f"bjet{i}_hMass"
    hist_bjet_Mass = ROOT.TH1F(bjets_Mass_hName, bjets_Mass_hName, 100, 0.0, 200.0)
    bjets_hMass_l.append(hist_bjet_Mass)

    lep_Pt_hName = f"lep{i}_hPt"
    hist_lep_Pt = ROOT.TH1F(lep_Pt_hName, lep_Pt_hName, 100, 0, 400)
    lep_hPt_l.append(hist_lep_Pt)

    lep_Eta_hName = f"lep{i}_hEta"
    hist_lep_Eta = ROOT.TH1F(lep_Eta_hName, lep_Eta_hName, 100, -4.0, 4.0)
    lep_hEta_l.append(hist_lep_Eta)

    lep_Phi_hName = f"lep{i}_hPhi"
    hist_lep_Phi = ROOT.TH1F(lep_Phi_hName, lep_Phi_hName, 100, 0.0, 3.2)
    lep_hPhi_l.append(hist_lep_Phi)


#histMass = ROOT.TH1F("mass", "M_{inv}(e_{1}, e_{2})", 100, 40.0, 140.0)

# Loop over all events
for entry in range(0, numberOfEntries):
  if (entry+1)%100 == 0:
    print (' ... processed {} events ...'.format(entry+1))

  # Load selected branches with data from specified event
  treeReader.ReadEntry(entry)
#  print("No of entries: ",entry)

  ## main MC event weight
  w =  branchEvent[0].Weight

  # If event contains at least 1 jet
  histJetsSize.Fill(branchJet.GetEntries(), w)
  if branchJet.GetEntries() > 3:

    ## 0 - Loose , 1 - Medium, 2 - Tight
    wp = 1
    bjets_list = []
    for bj in range(0,branchJet.GetEntries()):
        #    BtagOk = ( jet1.BTag & (1 << wp) )
        BtagOk = ( branchJet.At(bj).BTag )
        pt = branchJet.At(bj).PT
        eta = abs(branchJet.At(bj).Eta)
        phi = branchJet.At(bj).Phi
#    print("btag>>>>>>>>>>>>>>>>:",BtagOk)
    # Plot bjets transverse momentum
        if (BtagOk and pt > 30. and eta < 5.):
            bjets_list.append(branchJet.At(bj))
    histbJetsSize.Fill(len(bjets_list), w)
    print("Jets no: ",branchJet.GetEntries())
    if len(bjets_list) > 3:
        for hist in bjets_hPt_l:
            idx = bjets_hPt_l.index(hist)
            hist.Fill(branchJet.At(idx).PT, w)

        for hist in bjets_hEta_l:
            idx = bjets_hEta_l.index(hist)
            hist.Fill(branchJet.At(idx).Eta, w)

        for hist in bjets_hPhi_l:
            idx = bjets_hPhi_l.index(hist)
            hist.Fill(branchJet.At(idx).Phi, w)

        for hist in bjets_hMass_l:
            idx = bjets_hMass_l.index(hist)
            hist.Fill(branchJet.At(idx).Mass, w)




  # If event contains at least 4 leptons
  leptons_list = []
  if (branchElectron.GetEntries() + branchMuon.GetEntries() ) > 3:
    # Take first two electrons
    for e in range(0,branchElectron.GetEntries()):
        leptons_list.append(branchElectron.At(e))
    for m in range(0,branchMuon.GetEntries()):
        leptons_list.append(branchMuon.At(m))
    leptons_list.sort()
    for lep in lep_hPt_l:
      idx = lep_hPt_l.index(lep)
      lep.Fill(leptons_list[idx].PT, w)

    for lep in lep_hEta_l:
      idx = lep_hEta_l.index(lep)
      lep.Fill(leptons_list[idx].Eta, w)

    for lep in lep_hPhi_l:
      idx = lep_hPhi_l.index(lep)
      lep.Fill(leptons_list[idx].Phi, w)

#      print ("leptons list", leptons_list)
#      print ("lepton index: ", idx)
#    print("Electron no: ",branchElectron.GetEntries())
#    print("Muon no: ",branchMuon.GetEntries())
#  for i in bjets_list:
#    print ("bjets pt: ", i.PT)

    # Plot their invariant mass
 #   histMass.Fill(((elec1.P4()) + (elec2.P4())).M())

# Write the histograms to the ROOT file
file.Write()

# Close the file
file.Close()

