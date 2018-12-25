/*
#ifndef _readtree_h_
#define _readtree_h_

#include <TString.h>

void ReadTree(TString filename="chargediso_QCD_1000_MuEnriched_PU20bx25.root",TString output="plots/");
void RunOverSamples(TString output="plots/");

#endif
*/

#ifndef _readtree_h_
#define _readtree_h_

#include <TString.h>
#include <TGraph.h>
#include <TH1F.h>
#include "TLorentzVector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include <map>
#include <vector>

//enum FlavourSplitting {NOFLAVOURSPLITTING=0, UDSGSPLITTING=1, CSPLITTING=4, BSPLITTING=5 };

FactorizedJetCorrector *getFactorizedJetEnergyCorrector(TString,bool);
float getLeptonEnergyScaleUncertainty(int l_id,float l_pt,float l_eta);
Float_t computeMT(TLorentzVector &a, TLorentzVector &b);
void ReadTree(TString filename,
	      TString outname);
	      //	      Int_t channelSelection, 
	      //	      Int_t chargeSelection, 
	      //	      FlavourSplitting flavourSplitting,
	      //	      TH1F *normH, 
	      //	      Bool_t debug_);

std::map<Int_t,Float_t> lumiPerRun();
std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt);
float getLeptonEnergyScaleUncertainty(int l_id,float l_pt,float l_eta);
#endif
