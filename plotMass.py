#to run the script
#eosmount eos
#python plotMass.py
# import ROOT in batch mode
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
from ROOT import *
gROOT.SetBatch(True)
sys.argv = oldargv
gStyle.SetOptStat(0)
gStyle.SetOptFit(1)

# load FWLite C++ libraries
gSystem.Load("libFWCoreFWLite.so");
gSystem.Load("libDataFormatsFWLite.so");
AutoLibraryLoader.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

gen, genLabel = Handle("GenEventInfoProduct"), ("generator")

handlePruned, prunedLabel  = Handle ("std::vector<reco::GenParticle>"), ("prunedGenParticles")

bins = 500; lo = 40.0; hi = 140.0

h = TH1F("h","",bins,lo,hi); h.Sumw2()

dir = "eos/cms/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/70000/"

events = Events([
dir+"6E0A93BE-E15C-E611-8A59-002590D8C7FA.root",
dir+"084EDB90-1B5D-E611-87C9-002590D6012E.root",
dir+"AAED48BA-D95C-E611-B55E-141877637A57.root",
dir+"9415B4AD-E15C-E611-8A36-44A84225C3D0.root",
dir+"4E440A13-D45C-E611-8306-0025904CDDEC.root",
dir+"3E92BD12-D45C-E611-BA41-0025905C3D6A.root",
dir+"0AC259E3-1F5D-E611-B750-44A84225C7BB.root",
dir+"46594BF0-E95C-E611-822B-44A84225D36F.root",
dir+"262A03B5-255D-E611-8252-0025907277FE.root"
])

plotZ = False;

for i,event in enumerate(events):
  event.getByLabel(prunedLabel, handlePruned)
  event.getByLabel(genLabel, gen)
  
  if (gen.product().weight()>0.0): w = 1.0
  else: w = -1.0

  pruned = handlePruned.product()

  nz=0; nl=0; nt=0;

  ll=TLorentzVector(0,0,0,0)

  for p in pruned :

    if ( (abs(p.pdgId())==11 or abs(p.pdgId())==13) and p.status()==23):
      nl+=1
      l=TLorentzVector(p.px(),p.py(),p.pz(),p.energy())
      ll+=l;

    if (abs(p.pdgId())==15): 
      nt+=1

    if p.pdgId()==23 and p.status()==22:
      nz+=1
      if (plotZ): h.Fill(p.mass(),w)

  print 'nz',nz,'ntau',nt,'nl',nl,'ll.M',ll.M()
      
  if nt==0 and nl==2: 
    if (not plotZ): h.Fill(ll.M(),w)

print h.GetEntries()
h.SaveAs("genmz.root")

c1 = TCanvas("c1","c1",800,800)
c1.SetTopMargin(0.08)
c1.SetRightMargin(0.05)
c1.SetBottomMargin(0.12)
c1.SetLeftMargin(0.12)

c1.cd()

h.GetXaxis().SetTitle("m(Z)")
h.GetYaxis().SetTitle("d#sigma/dm")
h.GetYaxis().SetTitleOffset(1.5)

h.SetLineColor(1)
h.SetLineWidth(2)
h.Draw("ehist")

f1 = TF1("f1","[2]*TMath::BreitWigner(x,[0],[1])+[3]*TMath::Exp([4]*x)",50.0,140.0)
f1.SetParameter(0,91.2)
f1.SetParameter(1,2.5)
f1.SetNpx(1000)

h.Fit("f1","R")
f1.Draw("same")

latex2 = TLatex()
latex2.SetNDC()
latex2.SetTextSize(0.5*c1.GetTopMargin())
latex2.SetTextFont(42)
latex2.SetTextAlign(31) # align right                                                                                         
latex2.DrawLatex(0.92, 0.94,"13 TeV")
latex2.SetTextSize(0.6*c1.GetTopMargin())
latex2.SetTextFont(62)
latex2.SetTextAlign(11) # align right                                                                                         
latex2.DrawLatex(0.12, 0.94, "CMS")
latex2.SetTextSize(0.5*c1.GetTopMargin())
latex2.SetTextFont(52)
latex2.SetTextAlign(11)
latex2.DrawLatex(0.235, 0.94, "Simulation")

c1.SaveAs("genmz.pdf")
