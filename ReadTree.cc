#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include "UserCode/TopAnalysis/interface/MiniEvent.h"
#include "UserCode/TopAnalysis/interface/ReadTree.h"
#include "UserCode/TopAnalysis/interface/NeutrinoSolver.h"
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <vector>
#include <TLorentzVector.h>
#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidateFwd.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "UserCode/TopAnalysis/interface/BtagUncertaintyComputer.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include <math.h>
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <TFileMerger.h>

//void  NeutrinoSolver::NeutrinoSolver( TLorentzVector * aa, TLorentzVector *bb, double MW , double MT );
/*
std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt)
{
  std::vector<float> res(3,1.0);
  float ptSF(1.0), ptSF_err(0.0);
  if(TMath::Abs(eta)<0.8)       { ptSF=1.061; ptSF_err = 0.023; }
  else if(TMath::Abs(eta)<1.3)  { ptSF=1.088; ptSF_err = 0.029; }
  else if(TMath::Abs(eta)<1.9)  { ptSF=1.106; ptSF_err = 0.030; }
  else if(TMath::Abs(eta)<2.5)  { ptSF=1.126; ptSF_err = 0.094; }
  else if(TMath::Abs(eta)<3.0)  { ptSF=1.343; ptSF_err = 0.123; }
  else if(TMath::Abs(eta)<3.2)  { ptSF=1.303; ptSF_err = 0.111; }
  else if(TMath::Abs(eta)<5.0)  { ptSF=1.320; ptSF_err = 0.286; }
  res[0] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF)*(pt-genjpt)))/pt;
  res[1] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF-ptSF_err)*(pt-genjpt)))/pt;
  res[2] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF+ptSF_err)*(pt-genjpt)))/pt;
  return res;
}
*/

std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt)                                          
{                                                                                                                     
  std::vector<float> res(3,1.0);                                                                                      
  float ptSF(1.0), ptSF_err(0.0);                                                                                     
  if(TMath::Abs(eta)<0.5)       { ptSF=1.095; ptSF_err = 0.018; }                                                     
  else if(TMath::Abs(eta)<0.8)  { ptSF=1.120; ptSF_err = 0.028; }                                                     
  else if(TMath::Abs(eta)<1.1)  { ptSF=1.097; ptSF_err = 0.017; }                                                     
  else if(TMath::Abs(eta)<1.3)  { ptSF=1.103; ptSF_err = 0.033; }                                                     
  else if(TMath::Abs(eta)<1.7)  { ptSF=1.118; ptSF_err = 0.014; }                                                     
  else if(TMath::Abs(eta)<1.9)  { ptSF=1.100; ptSF_err = 0.033; }                                                     
  else if(TMath::Abs(eta)<2.1)  { ptSF=1.162; ptSF_err = 0.044; }                                                     
  else if(TMath::Abs(eta)<2.3)  { ptSF=1.160; ptSF_err = 0.048; }                                                     
  else if(TMath::Abs(eta)<2.5)  { ptSF=1.161; ptSF_err = 0.060; }                                                     
  else if(TMath::Abs(eta)<2.8)  { ptSF=1.209; ptSF_err = 0.059; }                                                     
  else if(TMath::Abs(eta)<3.0)  { ptSF=1.564; ptSF_err = 0.321; }                                                     
  else if(TMath::Abs(eta)<3.2)  { ptSF=1.384; ptSF_err = 0.033; }                                                     
  else if(TMath::Abs(eta)<5.0)  { ptSF=1.216; ptSF_err = 0.050; }                                                     
  res[0] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF)*(pt-genjpt)))/pt;                                           
  res[1] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF-ptSF_err)*(pt-genjpt)))/pt;                                  
  res[2] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF+ptSF_err)*(pt-genjpt)))/pt;                                  
  return res;                                                                                                         
}                                                                                                         

Float_t computeMT(TLorentzVector &a, TLorentzVector &b)
{
  return TMath::Sqrt(2*a.Pt()*b.Pt()*(1-TMath::Cos(a.DeltaPhi(b))));
}

bool sortBySignificance(std::pair<int,std::pair<float,float> > a,
			std::pair<int,std::pair<float,float> > b)
{
  if( a.second.first>0 || b.second.first>0 ) return (a.second.first>b.second.first);
  return (a.second.second>b.second.second);
}

void ReadTree(TString filename,  TString outname)
//void ReadTree(TString filename)
{
  gROOT->Reset();
  std::cout<< "filename inside read tree: " << filename <<std::endl;
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TTree *t = (TTree*)f->Get("demo/AnaTree");
  attachToMiniEventTree(t,ev);
  Int_t nentries(t->GetEntriesFast());
  t->GetEntry(0);
  TFile *file1 = new TFile(outname,"recreate");
  TTree *t1_ = new TTree("t1_", "title");
  MiniEvent_t evA_;
  createAnaTree(t1_,evA_);

  using namespace std;
  
  bool isDEl  ( filename.Contains("DoubleEG") );
  bool isDMu  ( filename.Contains("DoubleMu") );
  bool isMuEG ( filename.Contains("MuonEG") );
  bool isSEl  ( filename.Contains("SingleEl") );
  bool isSMu  ( filename.Contains("SingleMu") );
  
  //PILEUP WEIGHTING
  std::vector<TGraph *>puWgtGr;
  if(!ev.isData)
    {
      TString puWgtUrl("${CMSSW_BASE}/src/UserCode/TopAnalysis/data/pileupWgts.root");
      gSystem->ExpandPathName(puWgtUrl);
      TFile *fIn=TFile::Open(puWgtUrl);
      std::cout<< " file opened : " << std::endl;
      if(fIn)
	{
	  puWgtGr.push_back( (TGraph *)fIn->Get("puwgts_nom") );
	  puWgtGr.push_back( (TGraph *)fIn->Get("puwgts_down") );
	  puWgtGr.push_back( (TGraph *)fIn->Get("puwgts_up") );
	  fIn->Close();
	}
    }

  //BTAG SFs
  TString btagUncUrl("${CMSSW_BASE}/src/UserCode/TopAnalysis/data/CSVv2.csv");
  gSystem->ExpandPathName(btagUncUrl);
  std::vector<BTagCalibrationReader *> sfbReaders, sflReaders;
  TString btagEffExpUrl("${CMSSW_BASE}/src/UserCode/TopAnalysis/data/expTageff.root");
  gSystem->ExpandPathName(btagEffExpUrl);
  std::map<TString, TGraphAsymmErrors *> expBtagEff, expBtagEffPy8;
  BTagSFUtil myBTagSFUtil;
  if(!ev.isData)
    {
      BTagCalibration btvcalib("csvv2", btagUncUrl.Data());
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "central") );
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "down") ); 
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "up") );
      
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "central") );
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "down") ); 
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "up") );

      TFile *beffIn=TFile::Open(btagEffExpUrl);
      expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
      expBtagEff["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
      expBtagEff["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
      beffIn->Close();
    }

  //MUON EFFICIENCIES
  TString MuEffUrl("${CMSSW_BASE}/src/UserCode/TopAnalysis/data/leptonEfficiencies.root");
  gSystem->ExpandPathName(MuEffUrl);
  std::map<TString,TH2 *> mu_eff, el_eff;
  if(!ev.isData)
    {
      TFile *fIn=TFile::Open(MuEffUrl);
      mu_eff["m_sel"]=(TH2 *)fIn->Get("m_sel");
      for(auto& it : mu_eff) it.second->SetDirectory(0);
      fIn->Close();
    }

  //ELECTRON EFFICIENCIES
  TString ElEffUrl="${CMSSW_BASE}/src/UserCode/TopAnalysis/data/CutBasedID_TightWP_76X_18Feb.txt_SF2D.root";
  gSystem->ExpandPathName(ElEffUrl);
  if(!ev.isData)
    {
      TFile *fIn=TFile::Open(ElEffUrl);
      el_eff["e_sel"]=(TH2 *)fIn->Get("EGamma_SF2D");
      for(auto& it : el_eff) it.second->SetDirectory(0);
      fIn->Close();
    }
  /*
  //JET ENERGY SCALE: https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Summer15_uncertainties
  TString jecUncUrl("${CMSSW_BASE}/src/UserCode/TopAnalysis/data/Fall15_25nsV2_DATA_UncertaintySources_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  //FactorizedJetCorrector *jetCorr=getFactorizedJetEnergyCorrector("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/jectxt",!ev.isData);
  std::vector<TString> jecUncSrcs;
  std::vector<JetCorrectionUncertainty*> jecUncs;

  jecUncSrcs.push_back("AbsoluteStat");        
  jecUncSrcs.push_back("AbsoluteScale");        
  jecUncSrcs.push_back("AbsoluteFlavMap");        
  jecUncSrcs.push_back("AbsoluteMPFBias");        
  jecUncSrcs.push_back("Fragmentation");        
  jecUncSrcs.push_back("SinglePionECAL");  
  jecUncSrcs.push_back("SinglePionHCAL"); 
  jecUncSrcs.push_back("TimeEta");
  jecUncSrcs.push_back("TimePt");
  jecUncSrcs.push_back("RelativeJEREC1");  
  jecUncSrcs.push_back("RelativeJEREC2"); 
  jecUncSrcs.push_back("RelativeJERHF");
  jecUncSrcs.push_back("RelativePtBB");    
  jecUncSrcs.push_back("RelativePtEC1");  
  jecUncSrcs.push_back("RelativePtEC2");   
  jecUncSrcs.push_back("RelativePtHF");  
  jecUncSrcs.push_back("RelativeFSR");
  jecUncSrcs.push_back("RelativeStatEC"); 
  jecUncSrcs.push_back("RelativeStatHF");
  jecUncSrcs.push_back("PileUpDataMC");    
  jecUncSrcs.push_back("PileUpPtRef");    
  jecUncSrcs.push_back("PileUpPtBB");     
  jecUncSrcs.push_back("PileUpPtEC1");      
  jecUncSrcs.push_back("PileUpPtEC2");      
  jecUncSrcs.push_back("PileUpPtHF");    
  jecUncSrcs.push_back("FlavorPureGluon"); 
  jecUncSrcs.push_back("FlavorPureQuark");
  jecUncSrcs.push_back("FlavorPureCharm"); 
  jecUncSrcs.push_back("FlavorPureBottom");
  for(size_t i=0; i<jecUncSrcs.size(); i++)
    {
      JetCorrectorParameters *p = new JetCorrectorParameters(jecUncUrl.Data(), jecUncSrcs[i].Data());
      jecUncs.push_back( new JetCorrectionUncertainty(*p) );
    }
*/
  //LIST OF SYSTEMATICS
  std::vector<TString> expSysts;
  //  expSysts=jecUncSrcs;
  
  expSysts.push_back("JER");
  expSysts.push_back("Pileup");
  expSysts.push_back("Trigger");
  expSysts.push_back("MuEfficiency");
  expSysts.push_back("EleEfficiency");
  
  expSysts.push_back("BtagEff");
  expSysts.push_back("CtagEff");
  expSysts.push_back("LtagEff");
  
  Int_t nExpSysts=expSysts.size();
  cout << "\t..." << "/" << expSysts.size() << " generator level/experimental systematics will be considered" << endl;

  std::map<Int_t,Float_t> lumiMap;
  if(!ev.isData) lumiMap = lumiPerRun();

  // Book Histograms 
  std::map<TString, TH1 *> plots;
  std::map<TString, TH2 *> plots2d;

  plots["counter_"]         = new TH1F("counter_",";  ;Events" ,2, 0.,2.);   
  plots["weights_"]         = new TH1F("weights_",";  ;Events" ,2, 0.,2.);   
  plots["cutflow"]         = new TH1F("cutflow",";  ;Events" ,13, 0.,13.);   
  plots["cutflow"]->GetXaxis()->SetBinLabel(1,"All Events");
  plots["cutflow"]->GetXaxis()->SetBinLabel(2,"#trigger cut");
  plots["cutflow"]->GetXaxis()->SetBinLabel(3,"#pv cut");
  plots["cutflow"]->GetXaxis()->SetBinLabel(4,"#3 tight leptons");
  plots["cutflow"]->GetXaxis()->SetBinLabel(5,"#pre-sel(uuu)");
  plots["cutflow"]->GetXaxis()->SetBinLabel(6,"#pre-sel(eee)");
  plots["cutflow"]->GetXaxis()->SetBinLabel(7,"#pre-sel(uue)");
  plots["cutflow"]->GetXaxis()->SetBinLabel(8,"#pre-sel(eeu)");
  plots["cutflow"]->GetXaxis()->SetBinLabel(9,"#final sel(uuu)");
  plots["cutflow"]->GetXaxis()->SetBinLabel(10,"#final sel(eee)");
  plots["cutflow"]->GetXaxis()->SetBinLabel(11,"#final sel(uue)");
  plots["cutflow"]->GetXaxis()->SetBinLabel(12,"#final sel(eeu)");
  plots["cutflow"]->GetXaxis()->SetBinLabel(13,"#final sel(eeu)");

  plots["finalCutflow"]         = new TH1F("finalCutflow",";  ;Events" ,5, 0.,5.);   
  plots["finalCutflow"]->GetXaxis()->SetBinLabel(1,"total");
  plots["finalCutflow"]->GetXaxis()->SetBinLabel(2,"#mu#mu#mu");
  plots["finalCutflow"]->GetXaxis()->SetBinLabel(3,"eee");
  plots["finalCutflow"]->GetXaxis()->SetBinLabel(4,"#mu#mue");
  plots["finalCutflow"]->GetXaxis()->SetBinLabel(5,"ee#mu");

  plots["bjD_cutflow"] = new TH1F("bjD_cutflow",";  ;Events" ,5, 0.,5.);
  plots["bjU_cutflow"] = new TH1F("bjU_cutflow",";  ;Events" ,5, 0.,5.);
  plots["bjN_cutflow"] = new TH1F("bjN_cutflow",";  ;Events" ,5, 0.,5.);

  plots["JERD_cutflow"] = new TH1F("JERD_cutflow",";  ;Events" ,5, 0.,5.);
  plots["JERU_cutflow"] = new TH1F("JERU_cutflow",";  ;Events" ,5, 0.,5.);

  plots["presel_cutflow"]         = new TH1F("presel_cutflow",";  ;Events" ,5, 0.,5.);
  plots["presel_cutflow"]->GetXaxis()->SetBinLabel(1,"total");
  plots["presel_cutflow"]->GetXaxis()->SetBinLabel(2,"#mu#mu#mu");
  plots["presel_cutflow"]->GetXaxis()->SetBinLabel(3,"eee");
  plots["presel_cutflow"]->GetXaxis()->SetBinLabel(4,"#mu#mue");
  plots["presel_cutflow"]->GetXaxis()->SetBinLabel(5,"ee#mu");

  plots["cr0_cutflow"]         = new TH1F("cr0_cutflow",";  ;Events" ,5, 0.,5.);
  plots["cr0_cutflow"]->GetXaxis()->SetBinLabel(1,"total");
  plots["cr0_cutflow"]->GetXaxis()->SetBinLabel(2,"#mu#mu#mu");
  plots["cr0_cutflow"]->GetXaxis()->SetBinLabel(3,"eee");
  plots["cr0_cutflow"]->GetXaxis()->SetBinLabel(4,"#mu#mue");
  plots["cr0_cutflow"]->GetXaxis()->SetBinLabel(5,"ee#mu");

  if(lumiMap.size()) plots["ratevsrun_"] = new TH1F("ratevsrun_",";Run number; Events/pb",lumiMap.size(),0,lumiMap.size());
  Int_t runCtr(0);
  for(std::map<Int_t,Float_t>::iterator it=lumiMap.begin(); it!=lumiMap.end(); it++,runCtr++)
    plots["ratevsrun_"]->GetXaxis()->SetBinLabel(runCtr+1,Form("%d",it->first));


  plots["Nvtx_before_bc"]   = new TH1F("Nvtx_before_bc", "; Nvtx; Events" ,40, 0.0,40.0);   
  plots["Nvtx_after_bc"]    = new TH1F("Nvtx_after_bc","; Nvtx; Events" ,40, 0.0,40.0);   
  plots["Nvtx_before"]   = new TH1F("Nvtx_before", "; Nvtx; Events" ,40, 0.0,40.0);   
  plots["Nvtx_after"]    = new TH1F("Nvtx_after","; Nvtx; Events" ,40, 0.0,40.0);   
  plots["csv_bc"]    = new TH1F("csv_bc","; csv; Events" ,100, 0.,1.);   
  plots["csv_"]    = new TH1F("csv_","; csv; Events" ,100, 0.,1.);   
  plots["electron_iso"]  = new TH1F("electron_iso","; electron_iso; Events" ,100, 0.,0.2);   
  plots["muon_iso"]      = new TH1F("muon_iso","; muon_iso; Events" ,100, 0.,1.);   
  plots["gen_wgt_"]      = new TH1F("gen_wgt_","; generator weight; Events" ,20, -2.0,2.0);   
  plots["z_mass_uu"]     = new TH1F("z_mass_uu","; Z Mass ;Events" ,60, 60.,120.);   
  plots["z_pt_uu"]       = new TH1F("z_pt_uu","; Z Mass ;Events" ,100, 0.,200.);   
  plots["z_mass_ee"]     = new TH1F("z_mass_ee","; Z Mass ;Events" ,60, 60.,120.);
  plots["z_pt_ee"]       = new TH1F("z_pt_ee","; Z Mass ;Events" ,100, 0.,200.);
  plots["w_tm_munu"]     = new TH1F("w_tm_munu","; W tMass ;Events" ,150, 0.,150.);
  plots["w_tm_elnu"]     = new TH1F("w_tm_elnu","; W tMass ;Events" ,150, 0.,150.);
  plots["Wb_dphi"]       = new TH1F("Wb_dphi", "; #Delta #Phi (W,b)", 100, -5.0,5.0);
  plots["top_pt_mu"]     = new TH1F("top_pt_mu","; top pt ;Events" ,400, 0.,400.);   
  plots["top_tm_mu"]     = new TH1F("top_tm_mu","; top tMass ;Events" ,400, 0.,400.);   
  plots["top_invM_mu"]   = new TH1F("top_invM_mu","; top invMass ;Events" ,400, 0.,400.);   
  plots["top_pt_el"]     = new TH1F("top_pt_el","; top pt ;Events" ,400, 0.,400.);
  plots["top_tm_el"]     = new TH1F("top_tm_el","; top tMass ;Events" ,400, 0.,400.);
  plots["top_invM_el"]   = new TH1F("top_invM_el","; top invMass ;Events" ,400, 0.,400.);
  plots["f_top_tm_el"]   = new TH1F("f_top_tm_el","; top tMass ;Events" ,400, 0.,400.);
  plots["f_top_tm_mu"]   = new TH1F("f_top_tm_mu","; top tMass ;Events" ,400, 0.,400.);   
  plots["vetoLep_multi_"] = new TH1F("vetoLep_multi_","; Electron Multiplicity ;Events" ,6, 1.,7.);
  plots["bjets_multi"]   = new TH1F("bjets_multi","; bJets Multiplicity ;Events" ,2, 1.,3.);
  plots["bjets_multi"]->GetXaxis()->SetBinLabel(1,"bjet=1");
  plots["bjets_multi"]->GetXaxis()->SetBinLabel(2,"bjet=2");
  plots["ljets_multi"]   = new TH1F("ljets_multi","; lJets Multiplicity ;Events" ,5, 1.,6.);
  plots["ljets_multi"]->GetXaxis()->SetBinLabel(1,"lightjet=1");
  plots["ljets_multi"]->GetXaxis()->SetBinLabel(2,"lightjet=2");
  plots["ljets_multi"]->GetXaxis()->SetBinLabel(3,"lightjet=3");
  plots["ljets_multi"]->GetXaxis()->SetBinLabel(4,"lightjet=4");
  plots["ljets_multi"]->GetXaxis()->SetBinLabel(5,"lightjet=5");
  plots["jets_multi"]   = new TH1F("jets_multi","; lJets Multiplicity ;Events" ,5, 1.,6.);
  plots["jets_multi"]->GetXaxis()->SetBinLabel(1,"jet=1");
  plots["jets_multi"]->GetXaxis()->SetBinLabel(2,"jet=2");
  plots["jets_multi"]->GetXaxis()->SetBinLabel(3,"jet=3");
  plots["jets_multi"]->GetXaxis()->SetBinLabel(4,"jet=4");
  plots["jets_multi"]->GetXaxis()->SetBinLabel(5,"jet=5");
  plots["bjets_multi_bc"]   = new TH1F("bjets_multi_bc","; bJets Multiplicity ;Events" ,3, 1.,4.);
  plots["ljets_multi_bc"]   = new TH1F("ljets_multi_bc","; lJets Multiplicity ;Events" ,5, 1.,6.);
  plots["jets_multi_bc"]   = new TH1F("jets_multi_bc","; lJets Multiplicity ;Events" ,5, 1.,6.);
  plots["elec_multi_bc"]   = new TH1F("elec_multi_bc","; Electron Multiplicity ;Events" ,4, 1.,5.);
  plots["elec_multi_bc"]->GetXaxis()->SetBinLabel(1,"e");
  plots["elec_multi_bc"]->GetXaxis()->SetBinLabel(2,"ee");
  plots["elec_multi_bc"]->GetXaxis()->SetBinLabel(3,"eee");
  plots["elec_multi_bc"]->GetXaxis()->SetBinLabel(4,"eeee");

  plots["muon_multi_bc"]   = new TH1F("muon_multi_bc","; Muon Multiplicity ;Events" ,4, 1.,5.);
  plots["muon_multi_bc"]->GetXaxis()->SetBinLabel(1,"#mu");
  plots["muon_multi_bc"]->GetXaxis()->SetBinLabel(2,"#mu#mu");
  plots["muon_multi_bc"]->GetXaxis()->SetBinLabel(3,"#mu#mu#mu");
  plots["muon_multi_bc"]->GetXaxis()->SetBinLabel(4,"#mu#mu#mu#mu");

  plots["lepton_multi_bc"] = new TH1F("lepton_multi_bc","; Muon pt ;Events" ,4, 1.,5.);
  plots["lepton_multi_bc"]->GetXaxis()->SetBinLabel(1,"l");
  plots["lepton_multi_bc"]->GetXaxis()->SetBinLabel(2,"ll");
  plots["lepton_multi_bc"]->GetXaxis()->SetBinLabel(3,"lll");
  plots["lepton_multi_bc"]->GetXaxis()->SetBinLabel(4,"llll");
  plots["elec_multi_"]   = new TH1F("elec_multi_","; Electron Multiplicity ;Events" ,4, 1.,5.);
  plots["elec_multi_"]->GetXaxis()->SetBinLabel(1,"e");
  plots["elec_multi_"]->GetXaxis()->SetBinLabel(2,"ee");
  plots["elec_multi_"]->GetXaxis()->SetBinLabel(3,"eee");
  plots["elec_multi_"]->GetXaxis()->SetBinLabel(4,"eeee");
  plots["muon_multi_"]   = new TH1F("muon_multi_","; Muon Multiplicity ;Events" ,4, 1.,5.);
  plots["muon_multi_"]->GetXaxis()->SetBinLabel(1,"#mu");
  plots["muon_multi_"]->GetXaxis()->SetBinLabel(2,"#mu#mu");
  plots["muon_multi_"]->GetXaxis()->SetBinLabel(3,"#mu#mu#mu");
  plots["muon_multi_"]->GetXaxis()->SetBinLabel(4,"#mu#mu#mu#mu");

  plots["dR_emu_"] = new TH1F("dR_emu_","; #Delta R (#mu, e) ;Events" ,100, 0.,5.);
  plots["dR_emu_bc"] = new TH1F("dR_emu_bc","; #Delta R (#mu, e) ;Events" ,100, 0.,5.);
  plots["muon_pt_bc"]      = new TH1F("muon_pt_bc","; Muon pt ;Events" ,100, 0.,300.);
  plots["electron_pt_bc"]  = new TH1F("electron_pt_bc","; Electron pt ;Events" ,100, 0.,300.);
  plots["lmuon_pt_bc"]      = new TH1F("lmuon_pt_bc","; Muon pt ;Events" ,100, 0.,300.);
  plots["lmuon_eta_bc"]     = new TH1F("lmuon_eta_bc","; Muon eta ;Events" ,100, -3.0,3.0);
  plots["sublmuon_pt_bc"]      = new TH1F("sublmuon_pt_bc","; Muon pt ;Events" ,100, 0.,300.);
  plots["sublmuon_eta_bc"]     = new TH1F("sublmuon_eta_bc","; Muon eta ;Events" ,100, -3.0,3.0);
  plots["tmuon_pt_bc"]      = new TH1F("tmuon_pt_bc","; Muon pt ;Events" ,100, 0.,300.);
  plots["tmuon_eta_bc"]     = new TH1F("tmuon_eta_bc","; Muon eta ;Events" ,100, -3.0,3.0);
  plots["lelectron_pt_bc"]  = new TH1F("lelectron_pt_bc","; Electron pt ;Events" ,100, 0.,300.);
  plots["lelectron_eta_bc"] = new TH1F("lelectron_eta_bc","; Muon eta ;Events" ,100, -3.0,3.0);
  plots["sublelectron_pt_bc"]  = new TH1F("sublelectron_pt_bc","; Electron pt ;Events" ,100, 0.,300.);
  plots["sublelectron_eta_bc"] = new TH1F("sublelectron_eta_bc","; Muon eta ;Events" ,100, -3.0,3.0);
  plots["telectron_pt_bc"]  = new TH1F("telectron_pt_bc","; Electron pt ;Events" ,100, 0.,300.);
  plots["telectron_eta_bc"] = new TH1F("telectron_eta_bc","; Muon eta ;Events" ,100, -3.0,3.0);
  plots["jet_pt_"]   = new TH1F("jet_pt_","; Leading jet pt ;Events" ,100, 0.,400.);
  plots["jet_eta_"]  = new TH1F("jet_eta_","; Leading jet eta ;Events" ,100, -5.,5.);

  TString tag ="";
  for (int i=0; i<=4; ++i)
    {
      if(i==0) tag= ""; if(i==1) tag= "uuu"; if(i==2) tag= "uue";if(i==3) tag= "eeu";if(i==4) tag= "eee";

      plots["lmuon_pt_"+tag]      = new TH1F("lmuon_pt_"+tag,"; Muon pt ;Events" ,100, 0.,300.);   
      plots["lelectron_pt_"+tag]  = new TH1F("lelectron_pt_"+tag,"; Electron pt ;Events" ,100, 0.,300.);   
      plots["lmuon_eta_"+tag]     = new TH1F("lmuon_eta_"+tag,"; Muon eta ;Events" ,100, -3.0,3.0);   
      plots["lelectron_eta_"+tag] = new TH1F("lelectron_eta_"+tag,"; Muon eta ;Events" ,100, -3.0,3.0);   
      plots["lbjet_pt_"+tag]  = new TH1F("lbjet_pt_"+tag,"; Leading bjet pt ;Events" ,100, 0.,400.);
      plots["lbjet_eta_"+tag] = new TH1F("lbjet_eta_"+tag,"; Leading bjet #eta ;Events" ,100, -5.0,5.);
      plots["ljet_pt_"+tag]   = new TH1F("ljet_pt_"+tag,"; Leading jet pt ;Events" ,100, 0.,400.);
      plots["ljet_eta_"+tag]  = new TH1F("ljet_eta_"+tag,"; Leading jet eta ;Events" ,100, -5.,5.);
      plots["z_mass_"+tag]    = new TH1F("z_mass_"+tag,"; Z Mass ;Events" ,40, 70.,110.);   
      plots["z_pt_"+tag]      = new TH1F("z_pt_"+tag,";   Z pt ;Events" ,100, 0.,200.);   
      plots["w_tm_"+tag]      = new TH1F("w_tm_"+tag,";   W tMass ;Events" ,150, 0.,150.);   
      plots["met_"+tag]          = new TH1F("met_"+tag,"; MET ;Events" ,100, 0.,200.);   

      plots["bjets_multi_"+tag]   = new TH1F("bjets_multi_"+tag,"; bJets Multiplicity ;Events" ,2, 1.,3);
      plots["ljets_multi_"+tag]   = new TH1F("ljets_multi_"+tag,"; bJets Multiplicity ;Events" ,5, 1.,6);
      plots["f_lmuon_pt_"+tag]      = new TH1F("f_lmuon_pt_"+tag,"; Muon pt ;Events" ,100, 0.,300.);   
      plots["f_lelectron_pt_"+tag]  = new TH1F("f_lelectron_pt_"+tag,"; Electron pt ;Events" ,100, 0.,300.);   
      plots["f_lmuon_eta_"+tag]     = new TH1F("f_lmuon_eta_"+tag,"; Muon eta ;Events" ,100, -3.0,3.0);   
      plots["f_lelectron_eta_"+tag] = new TH1F("f_lelectron_eta_"+tag,"; Muon eta ;Events" ,100, -3.0,3.0);   
      plots["f_lbjet_pt_"+tag]   = new TH1F("f_lbjet_pt_"+tag,"; Leading bjet pt ;Events" ,100, 0.,400.);
      plots["f_lbjet_eta_"+tag]  = new TH1F("f_lbjet_eta_"+tag,"; Leading bjet #eta ;Events" ,100, -5.0,5.);
      plots["f_ljet_pt_"+tag]    = new TH1F("f_ljet_pt_"+tag,"; Leading jet pt ;Events" ,100, 0.,400.);
      plots["f_ljet_eta_"+tag]   = new TH1F("f_ljet_eta_"+tag,"; Leading jet eta ;Events" ,100, -5.,5.);
      plots["f_z_mass_"+tag]     = new TH1F("f_z_mass_"+tag,"; Z Mass ;Events" ,40, 70.,110.);   
      plots["f_z_pt_"+tag]       = new TH1F("f_z_pt_"+tag,";   Z pt ;Events" ,100, 0.,200.);   
      plots["f_w_tm_"+tag]       = new TH1F("f_w_tm_"+tag,";   W tMass ;Events" ,150, 0.,150.);   
      plots["f_top_pt_"+tag]     = new TH1F("f_top_pt_"+tag,"; top pt ;Events" ,100, 0.,400.);
      plots["f_top_tm_"+tag]     = new TH1F("f_top_tm_"+tag,"; top tMass ;Events" ,100, 0.,400.);
      plots["f_top_invM_"+tag]   = new TH1F("f_top_invM_"+tag,"; top invMass ;Events" ,100, 0.,400.);
      plots["f_met_"+tag]        = new TH1F("f_met_"+tag,"; MET ;Events" ,100, 0.,200.);

      plots["z_mass_"+tag+"_cr0"]    = new TH1F("z_mass_"+tag+"_cr0","; Z Mass ;Events" ,40, 70.,110.);
      plots["w_tm_"+tag+"_cr0"]      = new TH1F("w_tm_"+tag+"_cr0",";   W tMass ;Events" ,150, 0.,150.);   
      plots["met_"+tag+"_cr0"]          = new TH1F("met_"+tag+"_cr0","; MET ;Events" ,100, 0.,200.);   
      plots["jets_multi_"+tag+"_cr0"]   = new TH1F("jets_multi_"+tag+"_cr0","; lJets Multiplicity ;Events" ,2, 0.,2.);
    }
 
  TString itag ="";
  for (int i=0; i<=4; ++i)
    {
      if(i==0) itag= "3m0e"; if(i==1) itag= "2m1e";if(i==2) itag= "1m2e";if(i==3) itag= "0m3e"; if(i==4) itag = "";
      plots2d["metptshapes_"+itag+"_exp"]                  
	= new TH2F("metptshapes_"+itag+"_exp", ";Missing transverse energy [GeV];Events" , 100,0.,200., 2*nExpSysts,0,2*nExpSysts);

      plots2d["nbjetshapes_"+itag+"_exp"]                  
	= new TH2F("nbjetshapes_"+itag+"_exp", "; bJets Multiplicity ;Events" , 2,1,3, 2*nExpSysts,0,2*nExpSysts);

      plots2d["nljetshapes_"+itag+"_exp"]                  
	= new TH2F("nljetshapes_"+itag+"_exp", "; light Jets Multiplicity ;Events" , 5,1,6, 2*nExpSysts,0,2*nExpSysts);

      //      plots2d["sys_cutflow_"]                  
      //	= new TH2F("sys_cutflow", "; sys cutflow ;Events" , 5,0,5, 2*nExpSysts,0,2*nExpSysts);

      //  plots["bjets_multi"]   = new TH1F("bjets_multi","; bJets Multiplicity ;Events" ,2, 1.,3.);

      for(Int_t isyst=0; isyst<nExpSysts; isyst++)
	{
	  for(Int_t ivar=0; ivar<2; ivar++)
	    {
	      TString label(expSysts[isyst] + (ivar==0 ? "Down" : "Up"));
	      plots2d["metptshapes_"+itag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
	      plots2d["nbjetshapes_"+itag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
	      plots2d["nljetshapes_"+itag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
	      //      plots2d["sys_cutflow"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
	    }
	}
    }

  for (auto& it : plots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  for (auto& it : plots2d) { it.second->Sumw2(); it.second->SetDirectory(0); }
  
  //get the original number of events in the dataset
  TH1F *TuplesCutFlow=(TH1F *)f->Get("demo/counter");
  if(TuplesCutFlow) {
    Float_t TuplesEventsBin1=TuplesCutFlow->GetBinContent(1);
    Float_t TuplesEventsBin2=TuplesCutFlow->GetBinContent(2);
    plots["counter_"]->SetBinContent(1,TuplesEventsBin1);
    plots["counter_"]->SetBinContent(2,TuplesEventsBin2); 
  }

  //get the original number of events in the dataset
  TH1F *tempCutFlow=(TH1F *)f->Get("demo/weights");
  if(tempCutFlow) {
    Float_t Bin1=tempCutFlow->GetBinContent(1);
    Float_t Bin2=tempCutFlow->GetBinContent(2);
    plots["weights_"]->SetBinContent(1,Bin1);
    plots["weights_"]->SetBinContent(2,Bin2); 
  }

  double nlo=1.;
  /*
  if (filename.Contains("tZq_ll_") )    nlo = 0.2657343;
  //  if (filename.Contains("WZJets_TuneCUETP8M1_13TeV") )    nlo = 0.65763;
  if (filename.Contains("TTJets_TuneCUETP8M1") )    nlo = 0.3316263;
  if (filename.Contains("TTZToLLNuNu_") )        nlo = 0.464799;
  if (filename.Contains("TTWJetsToLNu_") )        nlo = 0.513428;
  if (filename.Contains("DYJetsToLL_M-10to50_") ) nlo = 0.7316367;
  if (filename.Contains("DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX") ) nlo = 0.6702222;
  std::cout << "NLO weight : " << nlo <<std::endl;
*/

  //loop over all entries
  for (Int_t i=0;i<nentries;i++)
    {
      t->GetEntry(i);
      plots["cutflow"]->Fill(0);

      //update pileup weights, if found
      float wgt=1.0;
      float genwgt=1.0;
      std::vector<float> puWeight(3,1.0);
      if(!ev.isData)
	{
	  if(puWgtGr.size())
	    {
	      puWeight[0]=puWgtGr[0]->Eval(ev.putrue);  
	      puWeight[1]=puWgtGr[1]->Eval(ev.putrue); 
	      puWeight[2]=puWgtGr[2]->Eval(ev.putrue);
	    }
	  wgt = puWeight[0];
	  genwgt=ev.generator_weight;
	  wgt*=genwgt;
	  wgt*=nlo;
	  evA_.wgt1=wgt;
	}
      //      std::cout<< "------Gen Weight --------"<< genwgt<<std::endl;
      plots["cutflow"]->Fill(1,genwgt);
      if(ev.tzq_nw>0)   plots["gen_wgt_"]->Fill(genwgt);

      if (!ev.isData && ev.muTrigger == 0 && ev.elTrigger == 0 && ev.muegTrigger == 0 && ev.singElTrigger == 0 && ev.singMuTrigger == 0) continue;


      bool DMu=false, DEl=false, MuEG= false, SMu=false, SEl=false;

      if (ev.muegTrigger !=0) MuEG=true;
      if (ev.muTrigger != 0 && ev.muegTrigger == 0) DMu=true;
      if (ev.elTrigger != 0 && ev.muegTrigger ==0 && ev.muTrigger == 0) DEl=true;
      if (ev.singElTrigger != 0 && ev.elTrigger == 0 && ev.muegTrigger ==0 && ev.muTrigger == 0) SEl=true;
      if (ev.singMuTrigger != 0 && ev.singElTrigger == 0 &&ev.elTrigger == 0 && ev.muegTrigger ==0 && ev.muTrigger == 0) SMu=true;


      std::cout<<" DMu: "<<DMu<<" DEl: "<<DEl<<" MuEG: "<<MuEG<<" SEl: "
	       <<SEl<<" SMu: "<<SMu<<std::endl;
      
      if(ev.isData && isDMu  && !DMu)  continue;
      if(ev.isData && isDEl  && !DEl)  continue;
      if(ev.isData && isMuEG && !MuEG) continue;
      if(ev.isData && isSMu  && !SMu ) continue;
      if(ev.isData && isSEl  && !SEl)  continue;
      
      std::cout<< "-------------after trigger ----------------"<<std::endl;

      MEzCalculator mu_neutrinoPz, el_neutrinoPz;
      std::vector<TLorentzVector> selectedMuons, selectedElectrons, selectedLeptons;
      std::vector<TLorentzVector> vetoMuons, vetoElectrons;
      std::vector<int> mu_charge, el_charge, mu_iso, el_iso; 
            
      //Selected Muons
      for(int m=0; m<ev.SelMu;m++)
	{
	  TLorentzVector tmp;
	  evA_.muoniso=ev.SelMu_iso[m];
	  //	  if(ev.SelMu_iso[m] > 0.12) continue;
	  tmp.SetPtEtaPhiE(ev.SelMu_pt[m], ev.SelMu_eta[m], ev.SelMu_phi[m], ev.SelMu_en[m]);
	  selectedMuons.push_back(tmp);
	  mu_charge.push_back(ev.SelMu_charge[m]);
	  mu_iso.push_back(ev.SelMu_iso[m]);
	  plots["muon_pt_bc"]->Fill(ev.SelMu_pt[m],wgt);
	  plots["muon_iso"]->Fill(ev.SelMu_iso[m]);	      
	}
      //Selected Electrons
      for(int e=0; e<ev.SelEl;e++)
	{
	  TLorentzVector tmp1;
	  evA_.electrononiso=ev.SelEl_iso[e];
	  //	  if(ev.SelEl_iso[e] > 0.0646) continue;
	  tmp1.SetPtEtaPhiE(ev.SelEl_pt[e], ev.SelEl_eta[e], ev.SelEl_phi[e], ev.SelEl_en[e]);
	  selectedElectrons.push_back(tmp1);
	  el_charge.push_back(ev.SelEl_charge[e]);
	  el_iso.push_back(ev.SelEl_iso[e]);
	  plots["electron_pt_bc"]->Fill(ev.SelEl_pt[e],wgt);
	  plots["electron_iso"]->Fill(ev.SelEl_iso[e]);
	}
      
      /*
      //Selected Non Muons
      for(int m=0; m<ev.selNonIsoMu;m++)
	{
	  TLorentzVector tmp;
	  tmp.SetPtEtaPhiE(ev.selNonIsoMu_pt[m], ev.selNonIsoMu_eta[m], ev.selNonIsoMu_phi[m], ev.selNonIsoMu_en[m]);
	  selectedMuons.push_back(tmp);
	  mu_charge.push_back(ev.selNonIsoMu_charge[m]);
	  plots["muon_pt_bc"]->Fill(ev.selNonIsoMu_pt[m],wgt);
	}

      //Selected Non Iso Electrons
      for(int e=0; e<ev.SelNonIsoEl;e++)
	{
	  TLorentzVector tmp1;
	  tmp1.SetPtEtaPhiE(ev.SelNonIsoEl_pt[e], ev.SelNonIsoEl_eta[e], ev.SelNonIsoEl_phi[e], ev.SelNonIsoEl_en[e]);
	  selectedElectrons.push_back(tmp1);
	  el_charge.push_back(ev.SelNonIsoEl_charge[e]);
	  plots["electron_pt_bc"]->Fill(ev.SelNonIsoEl_pt[e],wgt);
	}
*/

      //Veto Muons
      for(int m=0; m<ev.VetoMu;m++)
	{
	  TLorentzVector tmp2;
	  //	  if(ev.VetoMu_iso[m] > 0.25) continue;
	  tmp2.SetPtEtaPhiE(ev.VetoMu_pt[m], ev.VetoMu_eta[m], ev.VetoMu_phi[m], ev.VetoMu_en[m]);
	  vetoMuons.push_back(tmp2);
	  mu_charge.push_back(ev.VetoMu_charge[m]);
	}

      //Veto Electrons
      for(int e=0; e<ev.VetoEl;e++)
	{
	  TLorentzVector tmp3;
	  tmp3.SetPtEtaPhiE(ev.VetoEl_pt[e], ev.VetoEl_eta[e], ev.VetoEl_phi[e], ev.VetoEl_en[e]);
	  vetoElectrons.push_back(tmp3);
	  el_charge.push_back(ev.VetoEl_charge[e]);
	}

      int nSelMuons = selectedMuons.size();
      int nSelElectrons = selectedElectrons.size();
      int nSelLeptons = nSelElectrons + nSelMuons;
      int nvetoLeptons = vetoMuons.size()+vetoElectrons.size();
      plots["lepton_multi_bc"]->Fill(nSelLeptons,wgt);
      plots["vetoLep_multi_"]->Fill(nvetoLeptons,wgt);
      plots["cutflow"]->Fill(2,wgt);

      double lmuonpt(0), lmuoneta(0), lelectronpt(0), lelectroneta(0);
      double tmuonpt(0), tmuoneta(0), telectronpt(0), telectroneta(0);
      double sublmuonpt(0), sublmuoneta(0), sublelectronpt(0), sublelectroneta(0);

      if(nSelMuons >0.){
        lmuonpt   = selectedMuons[0].Pt();
        lmuoneta  = selectedMuons[0].Eta();
        plots["lmuon_pt_bc"]->Fill(lmuonpt,wgt);
        plots["lmuon_eta_bc"]->Fill(lmuoneta,wgt);
	plots["muon_multi_bc"]->Fill(nSelMuons,wgt);
      }

      if(nSelElectrons >0.){
        lelectronpt = selectedElectrons[0].Pt();
        lelectroneta = selectedElectrons[0].Eta();
        plots["lelectron_pt_bc"]->Fill(lelectronpt,wgt);
        plots["lelectron_eta_bc"]->Fill(lelectroneta,wgt);
	plots["elec_multi_bc"]->Fill(nSelElectrons,wgt);
      }

      if(nSelMuons >1.){
        sublmuonpt   = selectedMuons[1].Pt();
        sublmuoneta  = selectedMuons[1].Eta();
        plots["sublmuon_pt_bc"]->Fill(sublmuonpt,wgt);
        plots["sublmuon_eta_bc"]->Fill(sublmuoneta,wgt);
      }

      if(nSelElectrons >1.){
        sublelectronpt = selectedElectrons[1].Pt();
        sublelectroneta = selectedElectrons[1].Eta();
        plots["sublelectron_pt_bc"]->Fill(sublelectronpt,wgt);
        plots["sublelectron_eta_bc"]->Fill(sublelectroneta,wgt);
      }

      if(nSelMuons >2.){
        tmuonpt   = selectedMuons[2].Pt();
        tmuoneta  = selectedMuons[2].Eta();
        plots["tmuon_pt_bc"]->Fill(tmuonpt,wgt);
        plots["tmuon_eta_bc"]->Fill(tmuoneta,wgt);
      }

      if(nSelElectrons >2.){
        telectronpt = selectedElectrons[2].Pt();
        telectroneta = selectedElectrons[2].Eta();
        plots["telectron_pt_bc"]->Fill(telectronpt,wgt);
        plots["telectron_eta_bc"]->Fill(telectroneta,wgt);
      }

      plots["Nvtx_before_bc"]->Fill(ev.nvtx,genwgt);      
      plots["Nvtx_after_bc"]->Fill(ev.nvtx,wgt);      

      //THREE LEPTON SELECTION CUT
      if (nSelLeptons!=3) continue;
      if (nvetoLeptons != 3) continue;

      bool eee=false, mmm=false, mme=false, eem=false; 
      if(nSelMuons == 3 && nSelElectrons == 3) mmm=true; 
      if(nSelMuons == 2 && nSelElectrons == 1) mme=true; 
      if(nSelMuons == 1 && nSelElectrons == 2) eem=true; 
      if(nSelMuons == 0 && nSelElectrons == 3) eee=true; 
      std::cout<< "eee:"<<eee<<"mmm"<<mmm<<"eem:"<<eem<<"mme:"<<mme<<std::endl;

      //LEPTON SCALE FACTORS
      std::vector<float> lepTriggerSF(3,1.0),muSelSF(3,1.0), elSelSF(3,1.0);
      if(!ev.isData)
	{
	  if (nSelElectrons > 0.)
	    {
	      for(int el=0; el<nSelElectrons; el++)
		{
		  float minEtaForEff( el_eff["e_sel"]->GetXaxis()->GetXmin() );
		  float maxEtaForEff( el_eff["e_sel"]->GetXaxis()->GetXmax()-0.01 );
		  float etaForEff=TMath::Max(TMath::Min(float(fabs(selectedElectrons[el].Eta())),maxEtaForEff),minEtaForEff);
		  Int_t etaBinForEff=el_eff["e_sel"]->GetXaxis()->FindBin(etaForEff);

		  float minPtForEff( el_eff["e_sel"]->GetYaxis()->GetXmin() );
		  float maxPtForEff( el_eff["e_sel"]->GetYaxis()->GetXmax()-0.01 );
		  float ptForEff=TMath::Max(TMath::Min(float(selectedElectrons[el].Pt()),maxPtForEff),minPtForEff);
		  Int_t ptBinForEff=el_eff["e_sel"]->GetYaxis()->FindBin(ptForEff);

		  float elSF(el_eff["e_sel"]->GetBinContent(etaBinForEff,ptBinForEff));
		  float elSFUnc(el_eff["e_sel"]->GetBinError(etaBinForEff,ptBinForEff));
		  elSelSF[0]*=elSF; elSelSF[1]*=(elSF-elSFUnc); elSelSF[2]*=(elSF+elSFUnc);
		}
	    }

	  if (nSelMuons > 0.)
	    {
	      for(int mu=0; mu<nSelMuons; mu++)
		{
		  float minEtaForEff( mu_eff["m_sel"]->GetXaxis()->GetXmin() );
		  float maxEtaForEff( mu_eff["m_sel"]->GetXaxis()->GetXmax()-0.01 );
		  float etaForEff=TMath::Max(TMath::Min(float(fabs(selectedMuons[mu].Eta())),maxEtaForEff),minEtaForEff);
		  Int_t etaBinForEff=mu_eff["m_sel"]->GetXaxis()->FindBin(etaForEff);

		  float minPtForEff( mu_eff["m_sel"]->GetYaxis()->GetXmin() );
		  float maxPtForEff( mu_eff["m_sel"]->GetYaxis()->GetXmax()-0.01 );
		  float ptForEff=TMath::Max(TMath::Min(float(selectedMuons[mu].Pt()),maxPtForEff),minPtForEff);
		  Int_t ptBinForEff=mu_eff["m_sel"]->GetYaxis()->FindBin(ptForEff);

		  float muSF(mu_eff["m_sel"]->GetBinContent(etaBinForEff,ptBinForEff));
		  float muSFUnc(mu_eff["m_sel"]->GetBinError(etaBinForEff,ptBinForEff));
		  muSelSF[0]*=muSF; muSelSF[1]*=(muSF-muSFUnc); muSelSF[2]*=(muSF+muSFUnc);
		}
	    }

	  float trigSF(1.0), trigSFUnc(0.09);
	  lepTriggerSF[0]*=trigSF; lepTriggerSF[1]*=(trigSF-trigSFUnc); lepTriggerSF[2]*=(trigSF+trigSFUnc);
	  std::cout<< " elSelSF : " << elSelSF[0] << " muSelSF : " << muSelSF[0] << std::endl;
	  wgt*=lepTriggerSF[0]*muSelSF[0]*elSelSF[0];
	}

      //Jet Selection   
      Float_t htsum(0);
      std::vector<TLorentzVector> bJets,lightJets, allJets;
      TLorentzVector jetDiff(0,0,0,0);
      int nlightjets(0), nbjets(0), ncjets(0),nljets(0),nalljets(0), leadingJetIdx(-1);
      std::vector<int> resolvedJetIdx;
      std::vector<TLorentzVector> resolvedJetP4;
      for (int k=0; k<ev.nj;k++)
	{
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
	  resolvedJetIdx.push_back(k);
	  //	  jetDiff -= jp4;
	  float genJet_pt(ev.genj_pt[k]); 
	  if(!ev.isData && genJet_pt>0) 
	    {
	      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt)[0];
	      jp4 *= jerSmear;
	    }
	  //	  jetDiff += jp4;
	  resolvedJetIdx.push_back(k);
	  resolvedJetP4.push_back(jp4);
	  if(jp4.Pt()<=30) continue;
	  if(fabs(jp4.Eta()) > 4.5) continue;
	  if(leadingJetIdx<0) leadingJetIdx=k;
	  /*
	  float dR_jet_mu=999., dR_jet_el=999.;
	  for(unsigned int mu=0; mu <selectedMuons.size(); ++mu)
	    {
	      dR_jet_mu = jp4.DeltaR(selectedMuons[mu]); 
	      if(dR_jet_mu <= 0.4) break; 
	    }
	  if (dR_jet_mu <= 0.4) continue;
	  for(unsigned int el=0; el <selectedElectrons.size(); ++el)
	    {
	      dR_jet_el = jp4.DeltaR(selectedElectrons[el]); 
	      if(dR_jet_el <= 0.4) break; 
	    }
	  if (dR_jet_el <= 0.4) continue;
	  */
	  htsum     += jp4.Pt();

	  //b-tag
	  float csv = ev.j_csv[k];  
	  plots["csv_bc"]->Fill(csv,wgt);
	  //bool isBTagged(csv>0.935);         // Tight WP
	  bool isBTagged(csv>0.800); // Medium WP
	  //	  bool isBTagged(csv>0.460); // Loose WP
	  if(!ev.isData)
	    {
	      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
	      float expEff(1.0), jetBtagSF(1.0);
	      if(abs(ev.j_hadflav[k])==4) 
		{ 
		  ncjets++;
		  expEff    = expBtagEff["c"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff>0 ? 1 : 0. ; 
		}
	      else if(abs(ev.j_hadflav[k])==5) 
		{ 
		  nbjets++;
		  expEff    = expBtagEff["b"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff>0 ? 1 : 0. ; 
		}
	      else
		{
		  nljets++;
		  expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
                  jetBtagSF = sflReaders[0]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff>0 ? 1 : 0. ;
		}
	            
	      //updated b-tagging decision with the data/MC scale factor
	      myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
	    }

	  //save jet
	  allJets.push_back(jp4);
	  if(isBTagged) bJets.push_back(jp4);
	  else lightJets.push_back(jp4);
	}
      double lbjetpt(0), lbjeteta(0), ljetpt(0), ljeteta(0), jetpt(0), jeteta(0);
      double bjets_size=bJets.size();
      nlightjets = lightJets.size();
      nalljets = allJets.size();

      plots["Nvtx_before"]->Fill(ev.nvtx,genwgt);      
      plots["Nvtx_after"]->Fill(ev.nvtx,wgt);      
      plots["cutflow"]->Fill(3,wgt);
      //plots["csv_"]->Fill(csv,wgt);

      if(nSelMuons >0.){
	lmuonpt   = selectedMuons[0].Pt();
	lmuoneta  = selectedMuons[0].Eta();
	evA_.leadmuonpt=lmuonpt;
	evA_.leadmuoneta=lmuoneta;
	plots["lmuon_pt_"]->Fill(lmuonpt,wgt);      
	plots["lmuon_eta_"]->Fill(lmuoneta,wgt);
	plots["muon_multi_"]->Fill(nSelMuons,wgt);

      }      
      if(nSelElectrons >0.){
	lelectronpt = selectedElectrons[0].Pt();
	lelectroneta = selectedElectrons[0].Eta();
	evA_.leadelectronpt=lelectronpt;
	evA_.leadelectroneta=lelectroneta;
	plots["lelectron_pt_"]->Fill(lelectronpt,wgt);
	plots["lelectron_eta_"]->Fill(lelectroneta,wgt);
	plots["elec_multi_"]->Fill(nSelElectrons,wgt);
      }
      if(allJets.size()>0.){
	jetpt=allJets[0].Pt();
	jeteta=allJets[0].Eta();

	plots["jet_pt_"]->Fill(jetpt,wgt);
	plots["jet_eta_"]->Fill(jeteta,wgt);
	plots["jets_multi"]->Fill(nalljets,wgt);
      }
      if(bJets.size()>0.){
	lbjetpt=bJets[0].Pt();
	lbjeteta=bJets[0].Eta();
	evA_.leadbjetpt=lbjetpt;
	evA_.leadbjeteta=lbjeteta;
	plots["lbjet_pt_"]->Fill(lbjetpt,wgt);
	plots["lbjet_eta_"]->Fill(lbjeteta,wgt);
	plots["bjets_multi"]->Fill(bjets_size,wgt);
	if(mmm)	plots["bjets_multi_uuu"]->Fill(bjets_size,wgt);
	if(mme)	plots["bjets_multi_uue"]->Fill(bjets_size,wgt);
	if(eem)	plots["bjets_multi_eeu"]->Fill(bjets_size,wgt);
	if(eee)	plots["bjets_multi_eee"]->Fill(bjets_size,wgt);
      }

      if(lightJets.size()>0.){
	ljetpt=lightJets[0].Pt();
	ljeteta=lightJets[0].Eta();
	evA_.lightjetpt=ljetpt;
	evA_.lightjeteta=ljeteta;
	plots["ljet_pt_"]->Fill(ljetpt,wgt);
	plots["ljet_eta_"]->Fill(ljeteta,wgt);
	plots["ljets_multi"]->Fill(nlightjets,wgt);
	if(mmm)	plots["ljets_multi_uuu"]->Fill(nlightjets,wgt);
	if(mme)	plots["ljets_multi_uue"]->Fill(nlightjets,wgt);
	if(eem)	plots["ljets_multi_eeu"]->Fill(nlightjets,wgt);
	if(eee)	plots["ljets_multi_eee"]->Fill(nlightjets,wgt);
      }      

      float deta_emu=0., dphi_emu=0., dR_emu=999.;
      for(unsigned int mu=0; mu <selectedMuons.size(); ++mu)
	{
	  for(unsigned int el=0; el <selectedElectrons.size(); ++el)
	    {
	      deta_emu = selectedMuons[mu].Eta() - selectedElectrons[el].Eta();;
	      dphi_emu = deltaPhi(selectedMuons[mu].Phi(),  selectedElectrons[el].Phi());
	      dR_emu = sqrt( (deta_emu * deta_emu) + (dphi_emu * dphi_emu) );
	      if(dR_emu < 0.1) break;
	    }
	}
      plots["dR_emu_bc"]->Fill(dR_emu,wgt);
      if(dR_emu < 0.1) continue;
      plots["dR_emu_"]->Fill(dR_emu,wgt);

      std::map<Int_t,Float_t>::iterator rIt=lumiMap.find(ev.run);  
      if(rIt!=lumiMap.end()) 
	{
	  Int_t runCtr=std::distance(lumiMap.begin(),rIt);
	  plots["ratevsrun_"]->Fill(runCtr,1.e+6/rIt->second);
	}

      //Z Mass from Muons
      double min_zmass = 76., max_zmass= 106.;
      double min_zmass_wz = 81., max_zmass_wz= 101.;
      double z_mass = 91.1876, min_m =1000.0, ehard_pt =0., m_tmp, mmass_tmp;
      double min_e =1000.0, mhard_pt =0.,e_tmp, emass_tmp, zmass=0. , zpt =0;
      int j=0, jj=0, m_i =999, m_ii=999, e_i = 999, e_ii = 999;
      for(unsigned int mu1=0; mu1 <selectedMuons.size(); ++mu1){
	jj=0;
	for(unsigned int mu2=0; mu2 <selectedMuons.size(); ++mu2) {
	  if( (mu2 > mu1) &&  (mu_charge[mu1]*mu_charge[mu2]) < 0. ) {
	    //	    double dphi = (selectedMuons)[mu1].Phi() - (selectedMuons)[mu2].Phi();
	    //	    if( fabs(dphi)< 0.5) continue;
	    mmass_tmp  =  (selectedMuons[mu1] + selectedMuons[mu2]).M() ;
	    if( selectedMuons[mu1].Pt() > selectedMuons[mu2].Pt() ) {
	      mhard_pt =  selectedMuons[mu1].Pt();
	    }
	    else {
	      mhard_pt =  selectedMuons[mu2].Pt();
	    }
	    m_tmp =  z_mass - mmass_tmp;
	    if(fabs(m_tmp) <min_m  && mhard_pt > 20.)
	      {
		min_m = fabs(m_tmp); m_i = j; m_ii = jj;
	      }
	  }//end of sfos if loop
	  jj++;
	} //end of nested mu 2 loop
	j++;
      } //end of nested mu 1 loop
      
      std::pair < int,int> ZmuonPair;
      ZmuonPair= std::make_pair(m_i, m_ii);
      bool is2muon = false; j=0; jj=0;
      if(m_i != 999 && m_ii != 999) is2muon =true;
      if (is2muon )
	//      if (is2muon && mu_iso[m_i] < 0.15 && mu_iso[m_ii] < 0.15)
	{
	zmass  =  (selectedMuons[m_i] + selectedMuons[m_ii]).M() ;
	zpt    =  (selectedMuons[m_i] + selectedMuons[m_ii]).Pt() ;
	evA_.zmass_uu=zmass;
	evA_.zpt_uu=zpt;
	plots["z_mass_uu"]->Fill(zmass,wgt);
	plots["z_pt_uu"]->Fill(zpt,wgt);
	}

      ///For Electrons
      double e_zmass=0., e_zpt=0.;
      for(unsigned int el1=0; el1 <selectedElectrons.size(); ++el1){
	jj=0;
	for(unsigned int el2=0; el2 <selectedElectrons.size(); ++el2){
	  if( (el2 > el1) &&  (el_charge[el1]*el_charge[el2]) < 0. ){
	    //	    double dphi = (selectedElectrons)[el1].Phi() - (selectedElectrons)[el2].Phi();
	    //	    if( fabs(dphi)< 0.5) continue; 
	    emass_tmp  =  (selectedElectrons[el1] + selectedElectrons[el2]).M() ;
	    if( selectedElectrons[el1].Pt() > selectedElectrons[el2].Pt() ){
	      ehard_pt =  selectedElectrons[el1].Pt();
	    }
	    else {
	      ehard_pt =  selectedElectrons[el2].Pt();
	    }
	    e_tmp =  z_mass - emass_tmp;
	    if(fabs(e_tmp) < min_e  && ehard_pt > 20.) 
	      {
		min_e = fabs(e_tmp); e_i = j; e_ii = jj;
	      }
	  }//end of sfos if loop
	  jj++;
	} //end of nested mu 2 loop 
	j++;
      } //end of nested mu 1 loop 

      std::pair < int,int> ZelectronPair;
      ZelectronPair= std::make_pair(e_i, e_ii);
      bool is2electron = false;
      if(e_i != 999 && e_ii != 999) is2electron =true;
      if (is2electron)
	{
	  e_zmass  =  (selectedElectrons[e_i] + selectedElectrons[e_ii]).M() ;
	  e_zpt    =  (selectedElectrons[e_i] + selectedElectrons[e_ii]).Pt() ;
	  evA_.zmass_ee=e_zmass;
	  evA_.zpt_ee=e_zpt;
	  plots["z_mass_ee"]->Fill(e_zmass,wgt);
	  plots["z_pt_ee"]->Fill(e_zpt,wgt);
	}
      TLorentzVector met(0,0,0,0),  mup4_w (0,0,0,0), elp4_w (0,0,0,0); 
      bool isWmuon =false, isWelectron =false;
      if(is2electron || is2muon)
      	{
	  for( int mu=0; mu < (int) selectedMuons.size(); ++mu)
	    {
	      if( mu == ZmuonPair.first || mu == ZmuonPair.second) continue;
	      mup4_w = TLorentzVector((selectedMuons)[mu].Px(), (selectedMuons)[mu].Py(),
				      (selectedMuons)[mu].Pz(), (selectedMuons)[mu].Energy() );
	      isWmuon =true;
	    }
	  for( int el=0; el < (int) selectedElectrons.size(); ++el)
	    {
	      if( el == ZelectronPair.first || el == ZelectronPair.second) continue;
	      elp4_w = TLorentzVector((selectedElectrons)[el].Px(), (selectedElectrons)[el].Py(),
				      (selectedElectrons)[el].Pz(), (selectedElectrons)[el].Energy() );
	      isWelectron =true;
	    }
	}
           
      //MET and transverse mass
      double met_pt=0, mt_mu=0., mt_el=0; TLorentzVector mu_neutrino(0,0,0,0), el_neutrino(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt,0,ev.met_phi,0.);
      met.SetPz(0.); met.SetE(met.Pt());
      met_pt = ev.met_pt;
      plots["met_"]->Fill(met_pt,wgt);
      if(isWmuon)
	{
	  mt_mu = computeMT(mup4_w,met);
	  evA_.tm_mu=mt_mu;
	  plots["w_tm_munu"]->Fill(mt_mu,wgt);
	  mu_neutrinoPz.SetMET(met);
	  mu_neutrinoPz.SetLepton(mup4_w);
	  float mu_nupz=mu_neutrinoPz.Calculate();
	  //mu_nupz=0.;
	  float mu_nuE = sqrt( met.Px()*met.Px() + met.Py()*met.Py() + mu_nupz*mu_nupz);
	  mu_neutrino =  TLorentzVector(met.Px(), met.Py(), mu_nupz, mu_nuE );
	}

      if(isWelectron)
	{
	  mt_el = computeMT(elp4_w,met);
	  evA_.tm_el=mt_el;
	  plots["w_tm_elnu"]->Fill(mt_el,wgt);
	  el_neutrinoPz.SetMET(met);
	  el_neutrinoPz.SetLepton(elp4_w);
	  float el_nupz=el_neutrinoPz.Calculate();
	  //el_nupz=0.;
	  float el_nuE = sqrt(met.Px()*met.Px()+met.Py()*met.Py()+ el_nupz*el_nupz);
	  el_neutrino =  TLorentzVector(met.Px(), met.Py(), el_nupz, el_nuE);
	}

      /// Top reconstruction for Muon      
      TLorentzVector  Wmunu(0,0,0,0), top_Wmu (0,0,0,0);
      double top_mu=0, top_el=0;
      if(isWmuon == true && bjets_size == 1 )
	{
	  Wmunu = mup4_w + mu_neutrino;
	  //	  for(unsigned i=0; i< bJets.size();i++)
	  //	    {
	  top_Wmu =Wmunu + bJets.at(0);
	  top_mu = top_Wmu.M();
	  evA_.topmass_mu=top_mu;
	  evA_.toppt_mu=top_Wmu.Pt();
	  if(top_Wmu.Mt() > 0.)  plots["top_tm_mu"]->Fill(top_Wmu.Mt(),wgt);
	  if(top_mu       > 0.)  plots["top_invM_mu"]->Fill(top_mu,wgt);
	  if(top_Wmu.Pt() > 0.)  plots["top_pt_mu"]->Fill( top_Wmu.Pt(),wgt);
	      //	    }
	}
      //      TLorentzVector *aa, *bb;     
      //      TLorentzVector ee, metT_;     
      //      double metpx=met.Px(), metpy=met.Py(), Xerr=1.,Yerr=1., rhomet=0., te=1 ;
      /// Top reconstruction for Electron
      TLorentzVector  Welnu(0,0,0,0), top_Wel (0,0,0,0);
      if(isWelectron == true && bjets_size==1 )
        {
	  Welnu = elp4_w + el_neutrino;
	  //	  aa =  &elp4_w;
	  //	  bb=   &bJets.at(0);
	  //	  std::cout<<" aa : " << aa->Pt()<< " bb : " << bb->Pt()<<std::endl;
	  //	  NeutrinoSolver::NeutrinoSolver(elp4_w, bJets.at(0))
	  //	    NeutrinoSolver::Neutrino(const &elp4_w, const &bJets.at(0));
	  //	  double cc=78,dd=173;
	  //	  NeutrinoSolver NS( aa, bb,  cc ,  dd);
	  //	  metT_ = TLorentzVector(NS.GetBest(metpx, metpy, Xerr, Yerr, rhomet, te));
	  //	  Welnu = elp4_w + el_neutrino;
	  //	  std::cout<< "mine pz : " << el_neutrino.Pz() << "lulla Pz : " << metT_.Pz() << std::endl;
	  //	  plots["lulaa_pz"]->Fill(metT_.Pz(),wgt);
	  //	  plots["mine_pz"]->Fill(mu_neutrino.Pz(),wgt);
	  //	  for(unsigned i=0; i< bJets.size();i++)
	    //	    {
	  top_Wel =Welnu + bJets.at(0);
	  top_el = top_Wel.M();
	  evA_.topmass_el=top_el;
	  evA_.toppt_mu=top_Wel.Pt();
	  if(top_Wel.Mt() > 0.)  plots["top_tm_el"]->Fill(top_Wel.Mt(),wgt);
	  if(top_el       > 0.)  plots["top_invM_el"]->Fill(top_el,wgt);
	  if(top_Wel.Pt() > 0.)  plots["top_pt_el"]->Fill(top_Wel.Pt(),wgt);
	  //	    }
        }

      //--------------------------- Pre-Selection Plots-------------------------------
      bool isEEE=false, isMMM=false, isMME=false, isEEM=false; 
      if(selectedMuons.size() == 3 && is2muon && isWmuon && selectedElectrons.size() == 0 
	 && zmass > 0. && zmass>min_zmass && zmass < max_zmass) isMMM=true; 

      if(selectedMuons.size() == 0 && is2electron && isWelectron && selectedElectrons.size() == 3
	 && e_zmass > 0 && e_zmass >min_zmass && e_zmass < max_zmass )isEEE=true;
	
      if(selectedMuons.size() == 2 && is2muon && isWelectron && selectedElectrons.size() == 1
	 && zmass > 0 && zmass >min_zmass && zmass < max_zmass) isMME=true;

      if(selectedMuons.size() == 1 && is2electron && isWmuon && selectedElectrons.size() == 2
	 && e_zmass > 0 && e_zmass >min_zmass && e_zmass < max_zmass ) isEEM=true;

      if(isMMM || isEEE || isMME || isEEM)   plots["presel_cutflow"]->Fill(0.,wgt);

      if(isMMM)
	{
	  plots["presel_cutflow"]->Fill(1,wgt);
	  plots["cutflow"]->Fill(4,wgt);
	  plots["met_uuu"]->Fill(met_pt,wgt);
	  plots["lmuon_pt_uuu"]->Fill(lmuonpt,wgt);
	  plots["lmuon_eta_uuu"]->Fill(lmuoneta,wgt);
	  if(bJets.size() >0.)
	    {
	      plots["lbjet_pt_uuu"]->Fill(lbjetpt,wgt);
	      plots["lbjet_eta_uuu"]->Fill(lbjeteta,wgt);
	    }
	  if(lightJets.size() >0.)
	    {
	      plots["ljet_pt_uuu"]->Fill(ljetpt,wgt);
	      plots["ljet_eta_uuu"]->Fill(ljeteta,wgt);
	    }
	  plots["z_mass_uuu"]->Fill(zmass,wgt);
	  plots["z_pt_uuu"]->Fill(zpt,wgt);
	  plots["w_tm_uuu"]->Fill(mt_mu,wgt);
	}

      // eee
      if(isEEE)
	{
	  plots["presel_cutflow"]->Fill(2,wgt);
	  plots["cutflow"]->Fill(5,wgt);
	  plots["met_eee"]->Fill(met_pt,wgt);
	  plots["lelectron_pt_eee"]->Fill(lelectronpt,wgt);
	  plots["lelectron_eta_eee"]->Fill(lelectroneta,wgt);
	  if(bJets.size() >0.)
	    {
	      plots["lbjet_pt_eee"]->Fill(lbjetpt,wgt);
	      plots["lbjet_eta_eee"]->Fill(lbjeteta,wgt);
	    }
	  if(lightJets.size() >0.)
	    {
	      plots["ljet_pt_eee"]->Fill(ljetpt,wgt);
	      plots["ljet_eta_eee"]->Fill(ljeteta,wgt);
	    }
	  plots["z_mass_eee"]->Fill(e_zmass,wgt);
	  plots["z_pt_eee"]->Fill(e_zpt,wgt);
	  plots["w_tm_eee"]->Fill(mt_el,wgt);
	}

      // uue
      if(isMME)
	{
	  plots["presel_cutflow"]->Fill(3,wgt);
	  plots["cutflow"]->Fill(6,wgt);
	  plots["met_uue"]->Fill(met_pt,wgt);
	  plots["lmuon_pt_uue"]->Fill(lmuonpt,wgt);
	  plots["lmuon_eta_uue"]->Fill(lmuoneta,wgt);
	  plots["lelectron_pt_uue"]->Fill(lelectronpt,wgt);
	  plots["lelectron_eta_uue"]->Fill(lelectroneta,wgt);
	  if(bJets.size() >0.)
	    {
	      plots["lbjet_pt_uue"]->Fill(lbjetpt,wgt);
	      plots["lbjet_eta_uue"]->Fill(lbjeteta,wgt);
	    }
	  if(lightJets.size() >0.)
	    {
	      plots["ljet_pt_uue"]->Fill(ljetpt,wgt);
	      plots["ljet_eta_uue"]->Fill(ljeteta,wgt);
	    }
	  plots["z_mass_uue"]->Fill(zmass,wgt);
	  plots["z_pt_uue"]->Fill(zpt,wgt);
	  plots["w_tm_uue"]->Fill(mt_el,wgt);
	}

      // eeu
      if(isEEM)
	{
	  plots["presel_cutflow"]->Fill(4,wgt);
	  plots["cutflow"]->Fill(7,wgt);
	  plots["met_eeu"]->Fill(met_pt,wgt);
	  plots["lmuon_pt_eeu"]->Fill(lmuonpt,wgt);
	  plots["lmuon_eta_eeu"]->Fill(lmuoneta,wgt);
	  plots["lelectron_pt_eeu"]->Fill(lelectronpt,wgt);
	  plots["lelectron_eta_eeu"]->Fill(lelectroneta,wgt);
	  if(bJets.size() >0.)
	    {
	    plots["lbjet_pt_eeu"]->Fill(lbjetpt,wgt);
	    plots["lbjet_eta_eeu"]->Fill(lbjeteta,wgt);
	    }
	  if(lightJets.size() >0.)
	    {
	    plots["ljet_pt_eeu"]->Fill(ljetpt,wgt);
	    plots["ljet_eta_eeu"]->Fill(ljeteta,wgt);
	    }
	  plots["z_mass_eeu"]->Fill(e_zmass,wgt);
	  plots["z_pt_eeu"]->Fill(e_zpt,wgt);
	  plots["w_tm_eeu"]->Fill(mt_mu,wgt); 
	}


      //--------------------------- WZ Control Region -------------------------------

      bool isEEE_cr=false, isMMM_cr=false, isMME_cr=false, isEEM_cr=false; 
      if(selectedMuons.size() == 3 && is2muon && isWmuon && selectedElectrons.size() == 0 
	 && zmass > 0. && zmass>min_zmass_wz && zmass < max_zmass_wz) isMMM_cr=true; 

      if(selectedMuons.size() == 0 && is2electron && isWelectron && selectedElectrons.size() == 3
	 && e_zmass > 0 && e_zmass >min_zmass_wz && e_zmass < max_zmass_wz )isEEE_cr=true;
	
      if(selectedMuons.size() == 2 && is2muon && isWelectron && selectedElectrons.size() == 1
	 && zmass > 0 && zmass >min_zmass_wz && zmass < max_zmass_wz) isMME_cr=true;

      if(selectedMuons.size() == 1 && is2electron && isWmuon && selectedElectrons.size() == 2
	 && e_zmass > 0 && e_zmass >min_zmass_wz && e_zmass < max_zmass_wz ) isEEM_cr=true;


      if (nalljets <= 1 && bjets_size == 0 && met_pt > 30.)
	{ 
	  if( (mt_mu > 50. || mt_el > 50.) && (isMMM_cr || isEEE_cr || isMME_cr || isEEM_cr) )   plots["cr0_cutflow"]->Fill(0.,wgt);

	  if(isMMM_cr && mt_mu > 50.)
	    {
	      plots["jets_multi_uuu_cr0"]->Fill(nalljets,wgt);
	      plots["cr0_cutflow"]->Fill(1,wgt);
	      plots["met_uuu_cr0"]->Fill(met_pt,wgt);
	      plots["z_mass_uuu_cr0"]->Fill(zmass,wgt);
	      plots["w_tm_uuu_cr0"]->Fill(mt_mu,wgt);
	    }
	  // eee
	  if(isEEE_cr && mt_el > 50.)
	    {
	      plots["jets_multi_eee_cr0"]->Fill(nalljets,wgt);
	      plots["cr0_cutflow"]->Fill(2,wgt);
	      plots["met_eee_cr0"]->Fill(met_pt,wgt);
	      plots["z_mass_eee_cr0"]->Fill(e_zmass,wgt);
	      plots["w_tm_eee_cr0"]->Fill(mt_el,wgt);

	    }
	  // uue
	  if(isMME_cr && mt_el > 50.)
	    {
	      plots["jets_multi_uue_cr0"]->Fill(nalljets,wgt);
	      plots["cr0_cutflow"]->Fill(3,wgt);
	      plots["met_uue_cr0"]->Fill(met_pt,wgt);
	      plots["z_mass_uue_cr0"]->Fill(zmass,wgt);
	      plots["w_tm_uue_cr0"]->Fill(mt_el,wgt);

	    }
	  // eeu
	  if(isEEM_cr && mt_mu > 50.)
	    {
	      plots["jets_multi_eeu_cr0"]->Fill(nalljets,wgt);
	      plots["cr0_cutflow"]->Fill(4,wgt);
	      plots["met_eeu_cr0"]->Fill(met_pt,wgt);
	      plots["z_mass_eeu_cr0"]->Fill(e_zmass,wgt);
	      plots["w_tm_eeu_cr0"]->Fill(mt_mu,wgt);
	    }
	}

      //experimental systematics
      int MuonCat=TMath::Min((int)(nSelMuons),(int)3);
      int ElectronCat=TMath::Min((int)(nSelElectrons),(int)3);
      std::vector<TString> catsToFill (2,Form("%dm", MuonCat));
      catsToFill[1]+= Form("%de",ElectronCat);
      TString tag=catsToFill[1];
      //      std::cout<< "  tag: " << tag << std::endl; 

      for(size_t ivar=0; ivar<expSysts.size(); ivar++)
	{
	  TString varName=expSysts[ivar];
	  //	  std::cout << " ----------------varName :------------------- " << varName << std::endl;
	  //	  bool updateJEn(ivar<=jecUncSrcs.size());
	  bool updateBtag(varName.Contains("tagEff"));
	  bool updateJER(varName.Contains("JER"));
	  for(int isign=0; isign<2; isign++)
	    {
	      float newWgt(wgt);
	      //	      std::cout<<"weight : "<< wgt<< std::endl;
	      if(varName=="Pileup" && puWeight[0]!=0) 
		newWgt*=(isign==0 ? puWeight[1]/puWeight[0] : puWeight[2]/puWeight[0]);
	      if(varName=="Trigger") 
		newWgt *= (isign==0 ? lepTriggerSF[1]/lepTriggerSF[0] : lepTriggerSF[2]/lepTriggerSF[0]);
	      if(varName=="MuEfficiency") 
		newWgt *= (isign==0 ? muSelSF[1]/muSelSF[0] : muSelSF[2]/muSelSF[0]);
	      if(varName=="EleEfficiency") 
		newWgt *= (isign==0 ? elSelSF[1]/elSelSF[0] : elSelSF[2]/elSelSF[0]);

	      //jets
	      std::vector<TLorentzVector> varBJets,varLightJets;
	      TLorentzVector jetDiff(0,0,0,0),jetSum(0,0,0,0);
	      if( !(updateBtag || updateJER) )
		{ 
		  varBJets=bJets; 
		  varLightJets=lightJets; 
		}
	      else
		{
		  for (size_t ij=0; ij<resolvedJetIdx.size();ij++)
		    {
		      int k(resolvedJetIdx[ij]);
		      int jflav( abs(ev.j_hadflav[k]) );

		      //check kinematics
		      TLorentzVector jp4;
		      jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
		      jetDiff -= jp4;
		      /*
		      float dR_jet_mu_sys=999., dR_jet_el_sys=999.;
		      for(unsigned int mu=0; mu <selectedMuons.size(); ++mu)
			{
			  dR_jet_mu_sys = jp4.DeltaR(selectedMuons[mu]);
			  if(dR_jet_mu_sys <= 0.4) break;
			}
		      if (dR_jet_mu_sys <= 0.4) continue;
		      for(unsigned int el=0; el <selectedElectrons.size(); ++el)
			{
			  dR_jet_el_sys = jp4.DeltaR(selectedElectrons[el]);
			  if(dR_jet_el_sys <= 0.4) break;
			}
		      if (dR_jet_el_sys <= 0.4) continue;
		      */
		      //smear jet energy resolution for MC
		      float genJet_pt = ev.genj_pt[k];
		      if(!ev.isData && genJet_pt>0) 
			{
			  std::vector<float> jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt);
			  if(varName=="JER") jp4 *= jerSmear[isign+1];
			  else		     jp4 *= jerSmear[0];
			}
		      jetDiff += jp4;
		      jetSum += jp4;

		      //re-apply kinematics
		      if(jp4.Pt()<=30) continue;
		      if(fabs(jp4.Eta()) > 4.5) continue;

		      //b-tag
		      float csv = ev.j_csv[k];  
		      //	      bool isBTagged(csv>0.935); // Tight WP
		      bool isBTagged(csv>0.800); // Medium WP
		      //bool isBTagged(csv>0.460); // Loose WP
		      if(!ev.isData)
			{
			  float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
			  float expEff(1.0), jetBtagSF(1.0);
			  if(jflav==4) 
			    { 
			      expEff        = expBtagEff["c"]->Eval(jptForBtag); 
			      int idx(0);
			      if(varName=="CtagEff")   idx=(isign==0 ? 1 : 2);
			      jetBtagSF  = sfbReaders[idx]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
			      jetBtagSF *= expEff>0 ? 1 : 0. ;
			    }
			  else if(jflav==5)
			    { 
			      expEff=expBtagEff["b"]->Eval(jptForBtag); 
			      int idx(0);
			      if(varName=="BtagEff")   idx=(isign==0 ? 1 : 2);
			      jetBtagSF=sfbReaders[idx]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
			      jetBtagSF *= expEff>0 ? 1 : 0. ;
			    }
			  else
			    {
			      expEff=expBtagEff["udsg"]->Eval(jptForBtag);
			      int idx(0);
			      if(varName=="LtagEff")   idx=(isign==0 ? 1 : 2);
			      jetBtagSF=sflReaders[idx]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
			      jetBtagSF *= expEff>0 ? 1 : 0. ;
			    }
			        
			  //updated b-tagging decision with the data/MC scale factor
			  myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
			}

		      //save jet
		      if(isBTagged) varBJets.push_back(jp4);
		      else          varLightJets.push_back(jp4);
		    }
		}
	      //	      std::cout<< "new weight :" << newWgt <<std::endl;
	      plots2d["metptshapes_"+tag+"_exp"]->Fill(met_pt,2*ivar+isign,newWgt);
	      if (varBJets.size() > 0.)plots2d["nbjetshapes_"+tag+"_exp"]->Fill(varBJets.size(),2*ivar+isign,newWgt);
	      if (varLightJets.size() > 0.)plots2d["nljetshapes_"+tag+"_exp"]->Fill(varLightJets.size(),2*ivar+isign,newWgt);
	      if(varName=="BtagEff")
		{
		  if (isign==0 && varLightJets.size()>=1 && varBJets.size()==1)
		    {
		      if  (isMMM  || isEEE  || isMME  || isEEM )    plots["bjD_cutflow"]->Fill(0., newWgt);
		      if  (isMMM )    plots["bjD_cutflow"]->Fill(1., newWgt);
		      if  (isEEE )    plots["bjD_cutflow"]->Fill(2., newWgt);
		      if  (isMME )    plots["bjD_cutflow"]->Fill(3., newWgt);
		      if  (isEEM )    plots["bjD_cutflow"]->Fill(4., newWgt);
		    }

		  if (isign==1 && varLightJets.size()>=1 && varBJets.size()==1)
		    {
		      if  (isMMM  || isEEE  || isMME  || isEEM )    plots["bjU_cutflow"]->Fill(0., newWgt);
		      if  (isMMM )    plots["bjU_cutflow"]->Fill(1., newWgt);
		      if  (isEEE )    plots["bjU_cutflow"]->Fill(2., newWgt);
		      if  (isMME )    plots["bjU_cutflow"]->Fill(3., newWgt);
		      if  (isEEM )    plots["bjU_cutflow"]->Fill(4., newWgt);
		    }

		  if (  isign ==0 && bjets_size == 1 && nlightjets >= 1)
		    {
		      if  (isMMM  || isEEE  || isMME  || isEEM )    plots["bjN_cutflow"]->Fill(0., newWgt);
		      if  (isMMM )    plots["bjN_cutflow"]->Fill(1., newWgt);
		      if  (isEEE )    plots["bjN_cutflow"]->Fill(2., newWgt);
		      if  (isMME )    plots["bjN_cutflow"]->Fill(3., newWgt);
		      if  (isEEM )    plots["bjN_cutflow"]->Fill(4., newWgt);
		    }


		}
	      if(varName=="JER")
		{
		  if (isign==0 && varLightJets.size()>=1 && varBJets.size()==1)
		    {
		      if  (isMMM  || isEEE  || isMME  || isEEM )    plots["JERD_cutflow"]->Fill(0., newWgt);
		      if  (isMMM )    plots["JERD_cutflow"]->Fill(1., newWgt);
		      if  (isEEE )    plots["JERD_cutflow"]->Fill(2., newWgt);
		      if  (isMME )    plots["JERD_cutflow"]->Fill(3., newWgt);
		      if  (isEEM )    plots["JERD_cutflow"]->Fill(4., newWgt);
		    }

		  if (isign==1 && varLightJets.size()>=1 && varBJets.size()==1)
		    {
		      if  (isMMM  || isEEE  || isMME  || isEEM )    plots["JERU_cutflow"]->Fill(0., newWgt);
		      if  (isMMM )    plots["JERU_cutflow"]->Fill(1., newWgt);
		      if  (isEEE )    plots["JERU_cutflow"]->Fill(2., newWgt);
		      if  (isMME )    plots["JERU_cutflow"]->Fill(3., newWgt);
		      if  (isEEM )    plots["JERU_cutflow"]->Fill(4., newWgt);
		    }

		}

	    }
	}
      
      /// light jet step
      //      if (nalljets < 1) continue;
      if (nlightjets == 0) continue;

      //// bjet step
      if ( bjets_size !=1 ) continue;
      /// mwT step
      //      if ( !(mt_mu > 20.|| mt_el > 20.) ) continue;
      
      if(isMMM || isEEE || isMME || isEEM) plots["finalCutflow"]->Fill(0.,wgt);

      //--------------------------- FinalSelection Plots-------------------------------
      // uuu
      if(bjets_size == 1 && nlightjets >= 1 && isMMM)
	{
	  plots["finalCutflow"]->Fill(1,wgt);
	  plots["cutflow"]->Fill(8,wgt);
	  plots["f_met_uuu"]->Fill(met_pt,wgt);
	  plots["f_lmuon_pt_uuu"]->Fill(lmuonpt,wgt);
	  plots["f_lmuon_eta_uuu"]->Fill(lmuoneta,wgt);
	  plots["f_lbjet_pt_uuu"]->Fill(lbjetpt,wgt);
	  plots["f_lbjet_eta_uuu"]->Fill(lbjeteta,wgt);
	  plots["f_ljet_pt_uuu"]->Fill(ljetpt,wgt);
	  plots["f_ljet_eta_uuu"]->Fill(ljeteta,wgt);
	  plots["f_z_mass_uuu"]->Fill(zmass,wgt);
	  plots["f_z_pt_uuu"]->Fill(zpt,wgt);
	  plots["f_w_tm_uuu"]->Fill(mt_mu,wgt);
	  plots["f_top_tm_mu"]->Fill(top_Wel.Mt(),wgt);
	  plots["f_top_tm_uuu"]->Fill(top_Wmu.Mt(),wgt);
	  plots["f_top_invM_uuu"]->Fill(top_Wmu.M(),wgt);
	  plots["f_top_pt_uuu"]->Fill( top_Wmu.Pt(),wgt);
	}

      // eee
      if(bjets_size ==1 && nlightjets >= 1 && isEEE)
	{
	  plots["finalCutflow"]->Fill(2,wgt);
	  plots["cutflow"]->Fill(9);
	  plots["f_met_eee"]->Fill(met_pt,wgt);
	  plots["f_lelectron_pt_eee"]->Fill(lelectronpt,wgt);
	  plots["f_lelectron_eta_eee"]->Fill(lelectroneta,wgt);
	  plots["f_lbjet_pt_eee"]->Fill(lbjetpt,wgt);
	  plots["f_lbjet_eta_eee"]->Fill(lbjeteta,wgt);
	  plots["f_ljet_pt_eee"]->Fill(ljetpt,wgt);
	  plots["f_ljet_eta_eee"]->Fill(ljeteta,wgt);
	  plots["f_z_mass_eee"]->Fill(e_zmass,wgt);
	  plots["f_z_pt_eee"]->Fill(e_zpt,wgt);
	  plots["f_w_tm_eee"]->Fill(mt_el,wgt);
	  plots["f_top_tm_el"]->Fill(top_Wel.Mt(),wgt);
	  plots["f_top_tm_eee"]->Fill(top_Wel.Mt(),wgt);
	  plots["f_top_invM_eee"]->Fill(top_Wel.M(),wgt);
	  plots["f_top_pt_eee"]->Fill( top_Wel.Pt(),wgt);
	}

      // uue
      if(bjets_size ==1 && nlightjets >= 1 && isMME)
	{
	  plots["finalCutflow"]->Fill(3,wgt);
	  plots["cutflow"]->Fill(10);
	  plots["f_met_uue"]->Fill(met_pt,wgt);
	  plots["f_lmuon_pt_uue"]->Fill(lmuonpt,wgt);
	  plots["f_lmuon_eta_uue"]->Fill(lmuoneta,wgt);
	  plots["f_lelectron_pt_uue"]->Fill(lelectronpt,wgt);
	  plots["f_lelectron_eta_uue"]->Fill(lelectroneta,wgt);
	  plots["f_lbjet_pt_uue"]->Fill(lbjetpt,wgt);
	  plots["f_lbjet_eta_uue"]->Fill(lbjeteta,wgt);
	  plots["f_ljet_pt_uue"]->Fill(ljetpt,wgt);
	  plots["f_ljet_eta_uue"]->Fill(ljeteta,wgt);
	  plots["f_z_mass_uue"]->Fill(zmass,wgt);
	  plots["f_z_pt_uue"]->Fill(zpt,wgt);
	  plots["f_w_tm_uue"]->Fill(mt_el,wgt);
          plots["f_top_tm_uue"]->Fill(top_Wel.Mt(),wgt);
          plots["f_top_invM_uue"]->Fill(top_Wel.M(),wgt);
          plots["f_top_pt_uue"]->Fill( top_Wel.Pt(),wgt);
	}

      // eeu
      if(bjets_size ==1 && nlightjets >= 1 && isEEM)
	{
	  plots["cutflow"]->Fill(11,wgt);
	  plots["finalCutflow"]->Fill(4,wgt);
	  plots["f_met_eeu"]->Fill(met_pt,wgt);
	  plots["f_lmuon_pt_eeu"]->Fill(lmuonpt,wgt);
	  plots["f_lmuon_eta_eeu"]->Fill(lmuoneta,wgt);
	  plots["f_lelectron_pt_eeu"]->Fill(lelectronpt,wgt);
	  plots["f_lelectron_eta_eeu"]->Fill(lelectroneta,wgt);
	  plots["f_lbjet_pt_eeu"]->Fill(lbjetpt,wgt);
	  plots["f_lbjet_eta_eeu"]->Fill(lbjeteta,wgt);
	  plots["f_ljet_pt_eeu"]->Fill(ljetpt,wgt);
	  plots["f_ljet_eta_eeu"]->Fill(ljeteta,wgt);
	  plots["f_z_mass_eeu"]->Fill(e_zmass,wgt);
	  plots["f_z_pt_eeu"]->Fill(e_zpt,wgt);
	  plots["f_w_tm_eeu"]->Fill(mt_mu,wgt); 
          plots["f_top_tm_eeu"]->Fill(top_Wmu.Mt(),wgt);
          plots["f_top_invM_eeu"]->Fill(top_Wmu.M(),wgt);
          plots["f_top_pt_eeu"]->Fill( top_Wmu.Pt(),wgt);
	}

      t1_->Fill();
      cout<<" -----------------------------------------------------------------------" <<endl;
    }


  //file1->Close();

  f->Close();

  /*
  
  //save histos to file 
  //  TString outname ("given"); 
  TString selPrefix("out_");  
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  //  TFile *fOut=TFile::Open(dirName+"/"+selPrefix+baseName,"RECREATE");
  TFile *fOut=TFile::Open(dirName+"/"+selPrefix+baseName,"RECREATE");
  //  TFile *fOut = new TFile(selPrefix+outname,"recreate");
  //  TFile *fOut = new TFile("fOut.root","recreate");
  //  TDirectory *cdtree = fOut->mkdir("tree");
  //  cdtree->cd() ;
  fOut = t1_->GetCurrentFile(); //to get the pointer to the current file  
  fOut->Write();


  fOut->cd();
  for (auto& it : plots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : plots2d)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  //  fOut = t1_->GetCurrentFile();
  //  fOut->Add(t1_);
  fOut->Close();
  //
  */
  
  file1->cd();
  for (auto& it : plots)  {
    it.second->SetDirectory(file1); it.second->Write();
  }
  for (auto& it : plots2d)  {
    it.second->SetDirectory(file1); it.second->Write();
  }
  file1 = t1_->GetCurrentFile(); //to get the pointer to the current file
  file1->Write();

  file1->Close();
}

  std::map<Int_t,Float_t> lumiPerRun()
  {
    std::map<Int_t,Float_t> lumiMap;
    lumiMap [254231] =  33138.485    ;
    lumiMap [254232] =  110184.436   ;
    lumiMap[254790] =  11605943.026 ;
    lumiMap[254852] =  959960.760   ;
    lumiMap[254879] =  1895968.466  ;
    lumiMap[254906] =  1597591.232  ;
    lumiMap[254907] =  1091487.189  ;
    lumiMap[254914] =  993829.670   ;
    lumiMap[256630] =  1092761.593  ;
    lumiMap[256673] =  99832.950    ;
    lumiMap[256674] =  99302.405    ;
    lumiMap[256675] =  7849524.162  ;
    lumiMap[256676] =  9793747.718  ;
    lumiMap[256677] =  18367246.351 ;
    lumiMap[256801] =  9618463.129  ;
    lumiMap[256842] =  93685.940    ;
    lumiMap[256843] =  40362166.080 ;
    lumiMap[256866] =  757843.030   ;
    lumiMap[256867] =  4913613.818  ;
    lumiMap[256868] =  24546225.747 ;
    lumiMap[256869] =  1651719.692  ;
    lumiMap[256926] =  1689563.183  ;
    lumiMap[256941] =  9411864.875  ;
    lumiMap[257461] =  3427246.292  ;
    lumiMap[257531] =  9394696.392  ;
    lumiMap[257599] =  5564981.768  ;
    lumiMap[257613] =  82731544.165 ;
    lumiMap[257614] =  919074.018   ;
    lumiMap[257645] =  68160020.189 ;
    lumiMap[257682] =  14569734.275 ;
    lumiMap[257722] =  956606.613   ;
    lumiMap[257723] =  7019337.030  ;
    lumiMap[257735] =  587615.257   ;
    lumiMap[257751] =  29680404.147 ;
    lumiMap[257804] =  1317383.953  ;
    lumiMap[257805] =  18741121.322 ;
    lumiMap[257816] =  26637823.645 ;
    lumiMap[257819] =  16590859.087 ;
    lumiMap[257968] =  18564612.656 ;
    lumiMap[257969] =  42819344.139 ;
    lumiMap[258129] =  6426942.093  ;
    lumiMap[258136] =  3981543.789  ;
    lumiMap[258157] =  4246566.735  ;
    lumiMap[258158] =  115154237.210;
    lumiMap[258159] =  27649017.003 ;
    lumiMap[258177] =  115285544.224;
    lumiMap[258211] =  7181306.111  ;
    lumiMap[258213] =  12774821.409 ;
    lumiMap[258214] =  16794100.825 ;
    lumiMap[258215] =  454929.891   ;
    lumiMap[258287] =  14850896.885 ;
    lumiMap[258403] =  16927279.430 ;
    lumiMap[258425] =  11442275.566 ;
    lumiMap[258426] =  830918.449   ;
    lumiMap[258427] =  8809409.741  ;
    lumiMap[258428] =  12769721.390 ;
    lumiMap[258432] =  306090.456   ;
    lumiMap[258434] =  33448313.925 ;
    lumiMap[258440] =  48914089.775 ;
    lumiMap[258444] =  2256161.978  ;
    lumiMap[258445] =  17746898.812 ;
    lumiMap[258446] =  8067709.968  ;
    lumiMap[258448] =  38407595.657 ;
    lumiMap[258655] =  799765.962   ;
    lumiMap[258656] =  28462159.130 ;
    lumiMap[258694] =  17462408.905 ;
    lumiMap[258702] =  32814046.039 ;
    lumiMap[258703] =  34738929.010 ;
    lumiMap[258705] =  8450632.729  ;
    lumiMap[258706] =  57585416.124 ;
    lumiMap[258712] =  37898277.818 ;
    lumiMap[258713] =  11128389.088 ;
    lumiMap[258714] =  4570661.608  ;
    lumiMap[258741] =  5308753.892  ;
    lumiMap[258742] =  67529699.325 ;
    lumiMap[258745] =  23504003.502 ;
    lumiMap[258749] =  49178502.971 ;
    lumiMap[258750] =  15670309.877 ;
    lumiMap[259626] =  11677798.109 ;
    lumiMap[259637] =  16238521.437 ;
    lumiMap[259681] =  3522075.932  ;
    lumiMap[259683] =  7916719.783  ;
    lumiMap[259685] =  57452782.772 ;
    lumiMap[259686] =  27567315.296 ;
    lumiMap[259721] =  12766726.230 ;
    lumiMap[259809] =  14951741.295 ;
    lumiMap[259810] =  10165998.080 ;
    lumiMap[259811] =  7791667.150  ;
    lumiMap[259813] =  865658.417   ;
    lumiMap[259817] =  430314.760   ;
    lumiMap[259818] =  13593974.896 ;
    lumiMap[259820] =  13023747.118 ;
    lumiMap[259821] =  16688292.395 ;
    lumiMap[259822] =  33464085.107 ;
    lumiMap[259861] =  6786586.562  ;
    lumiMap[259862] =  47268844.797 ;
    lumiMap[259884] =  7143126.771  ;
    lumiMap[259890] =  10096941.110 ;
    lumiMap[259891] =  9999933.427  ;
    lumiMap[260373] =  11151316.640 ;
    lumiMap[260424] =  69025255.475 ;
    lumiMap[260425] =  24359916.200 ;
    lumiMap[260426] =  45105840.468 ;
    lumiMap[260427] =  16455592.554 ;
    lumiMap[260431] =  36097791.380 ;
    lumiMap[260532] =  71167380.483 ;
    lumiMap[260533] =  1271083.770  ;
    lumiMap[260534] =  33011855.207 ;
    lumiMap[260536] =  14853812.840 ;
    lumiMap[260538] =  23008065.482 ;
    lumiMap[260541] =  1902645.781  ;
    lumiMap[260575] =  1921606.212  ;
    lumiMap[260576] =  17923769.029 ;
    lumiMap[260577] =  8592324.703  ;
    lumiMap[260593] =  38029477.673 ;
    lumiMap[260627] =  186195320.533;
    return lumiMap;
  }
