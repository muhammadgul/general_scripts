#include <iostream>
#include <fstream>
#include <string>
TString path= "/user/mgul/work/iso_gr15_Mtgr50_root_files/";
//TFile *ttfile= new TFile(path+"output_MC13TeV_TTJets_powheg.root");
TFile *ttfile= new TFile(path+"output_new_SampleMC13TeV_TTJets_powheg.root");
TFile *ttWfile= new TFile(path+"output_MC13TeV_TTWJetsToLNu.root");
TFile *ttZfile= new TFile(path+"output_MC13TeV_TTZToLLNuNu_M-10.root");
//TFile *wjfile= new TFile(path+"output_newMadGraph_WJetsToLNu_TuneCUETP8M1_pythia8.root");
TFile *wjfile= new TFile(path+"output_MC13TeV_WJets.root");
//TFile *DYfile= new TFile(path+"output_newMadGraph_DYJetsToLL_M_50_TuneCUETP8M1_pythia8.root");
TFile *DYfile= new TFile(path+"output_DYJetsToLL_M_50.root");

TFile *file_tch_t= new TFile(path+"output_MC13TeV_SingleT_t.root");
TFile *file_tw_at= new TFile(path+"output_MC13TeV_SingleTbar_tW_nohad.root");
TFile *file_tw_t= new TFile(path+"output_MC13TeV_SingleT_tW_nohad.root");
TFile *file_sch= new TFile(path+"output_MC13TeV_SingleTbar_tW_nohad.root");

TFile *fileww= new TFile(path+"output_MC13TeV_WWToLNuQQ.root");
TFile *filezz= new TFile(path+"output_MC13TeV_ZZ.root");
TFile *filewz= new TFile(path+"output_MC13TeV_WZ.root");

TFile *fileqcd15= new TFile(path+"output_QCD_Pt_15to20.root");
TFile *fileqcd20= new TFile(path+"output_QCD_Pt_20to30.root");
TFile *fileqcd30= new TFile(path+"output_QCD_Pt_30to50.root");
TFile *fileqcd50= new TFile(path+"output_QCD_Pt_50to80.root");
TFile *fileqcd80= new TFile(path+"output_QCD_Pt_80to120.root");
TFile *fileqcd120= new TFile(path+"output_QCD_Pt_120to170.root");
TFile *fileqcd170= new TFile(path+"output_QCD_Pt_170to300.root");
TFile *fileqcd300= new TFile(path+"output_QCD_Pt_300to470.root");
TFile *fileqcd470= new TFile(path+"output_QCD_Pt_470to600.root");
TFile *fileqcd600= new TFile(path+"output_QCD_Pt_600to800.root");
TFile *fileqcd800= new TFile(path+"output_QCD_Pt_800to1000.root");
TFile *fileqcd1000= new TFile(path+"output_QCD_Pt_1000toInf.root");

TFile *datafB= new TFile(path+"output_Data13TeV_SingleMuon_2016B.root");
TFile *datafC= new TFile(path+"output_Data13TeV_SingleMuon_2016C.root");
TFile *datafD= new TFile(path+"output_Data13TeV_SingleMuon_2016D.root");
//void plot_shape(TString name,TString dir, TString xtitle, TString ytitle, int binning = 1 )
void format_h(TH1F* h, int fillcolor){
     h->SetLineWidth(1.5);
     h->SetLineColor(fillcolor);
     h->SetMarkerStyle(20);
     h->SetFillColor(fillcolor);
     }
void format_hSignal(TH1F* h, int fillcolor){
     h->SetLineWidth(1.);
     h->SetLineColor(fillcolor);
 //    h->SetMarkerStyle(20);
 //    h->SetFillColor(fillcolor);
     }

void xy_axis(TH1F* h, TString xtitle, TString ytitle){
      h->GetXaxis()->SetTitle(xtitle);
      h->GetXaxis()->SetTitleSize(0.028);
      h->GetXaxis()->SetLabelSize(0.03);
      h->GetXaxis()->SetTitleOffset(1.35);
      h->GetYaxis()->SetTitle(ytitle);
      h->GetYaxis()->SetTitleSize(0.028);
      h->GetYaxis()->SetLabelSize(0.03);
      h->GetYaxis()->SetTitleOffset(1.65);
      }
void state (TPaveStats* p, float x1, float x2, float y1, int color){
      p->SetX1NDC(x1);
      p->SetX2NDC(x2);
      p->SetY1NDC(y1);
      p->SetTextColor(color);
      }
Int_t cut_events(TH1F* h){
      TString nam=h->GetName();
      Int_t n_bins=h->GetSize();
      for (Int_t i=0; i<=n_bins; i++){
      if (nam=="events_eachCut"){
      Int_t entr1=h->GetBinContent(i+1);
      std::cout<<entr1[i+1]<<"   ";
       }
      }}  
      
void bin_name(TH1F* h2){TString naming[]=
      {"HLT","elec=0","mu=1","M_{T}<50 GeV","Jets >= 4","bjets >=2"};
      TString nam=h2->GetName();
cout<<"Hi I am here: "<<endl;
      for (int i=0;i<sizeof(naming)/sizeof(TString);i++){
      h2->GetXaxis()->SetBinLabel(i+1,naming[i]);
      }}
      
      void overflow(TH1 *h){
      Int_t nx    = h->GetNbinsX()+1;
      h->SetBinContent(nx-1,h->GetBinContent(nx)+h->GetBinContent(nx-1));
      h->SetBinContent(nx,0);
      h->SetBinContent(1,h->GetBinContent(0)+h->GetBinContent(1));
      h->SetBinContent(0,0);
      }

void plot_shape(TString name, TString dir, TString xtitle, TString ytitle)//, int binning = 1 )
{

//  gROOT->ForceStyle();
  TCanvas *c1 = new TCanvas("c1", "H to ttbar",22,100,623,500);  // 3rd option expansion on x-axis, 4th option on y axis
  THStack *hs = new THStack("hs"," ");
  TH1::SetDefaultSumw2(true);
   gStyle->SetOptStat(0);
   c1->Range(-33.02326,-282.867,254.8837,2508.087);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.1147011);
   c1->SetRightMargin(0.19063);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
//   c1->SetLogy(kTRUE);
//   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();

      TH1F *htt = (TH1F*)ttfile->Get(dir+name);
      format_h(htt,kOrange+1);
      TH1F *httW = (TH1F*)ttWfile->Get(dir+name);
      format_h(httW,kOrange+5);
      TH1F *httZ = (TH1F*)ttZfile->Get(dir+name);
      format_h(httZ,kOrange+5);

      TH1F *hwj = (TH1F*)wjfile->Get(dir+name);
      format_h(hwj,kOrange);
      TH1F *hDY = (TH1F*)DYfile->Get(dir+name);
      format_h(hDY,kBlue+1);

      TH1F *hst_tch_t = (TH1F*)file_tch_t->Get(dir+name);
      format_h(hst_tch_t,kMagenta);
      TH1F *hst_twch_at = (TH1F*)file_tw_at->Get(dir+name);
      format_h(hst_twch_at,kMagenta);
      TH1F *hst_twch_t = (TH1F*)file_tw_t->Get(dir+name);
      format_h(hst_twch_t,kMagenta);
      TH1F *hst_sch = (TH1F*)file_sch->Get(dir+name);
      format_h(hst_sch,kMagenta);


      TH1F *hww = (TH1F*)fileww->Get(dir+name);
      format_h(hww,kGreen+1);
      TH1F *hzz = (TH1F*)filezz->Get(dir+name);
      format_h(hzz,kGreen+1);
      TH1F *hwz = (TH1F*)filewz->Get(dir+name);
      format_h(hwz,kGreen+1);
/*
for (unsigned int i=1;i<12; i++){
  TH1F *hqcd[i] = (TH1F*)fileqcd[i]->Get(dir+name);
  format_h(hqcd[i],0);
  if (nam=="events_counter")cout<<("this is events  : hqcd%i",i)<<hqcd[i]->GetBinContent(1)<<endl;
   }
*/    

      TH1F *hqcd15 = (TH1F*)fileqcd15->Get(dir+name);
      TH1F *hqcd20 = (TH1F*)fileqcd20->Get(dir+name);
      TH1F *hqcd30 = (TH1F*)fileqcd30->Get(dir+name);
      TH1F *hqcd50 = (TH1F*)fileqcd50->Get(dir+name);
      TH1F *hqcd80 = (TH1F*)fileqcd80->Get(dir+name);
      TH1F *hqcd120 = (TH1F*)fileqcd120->Get(dir+name);
      TH1F *hqcd170 = (TH1F*)fileqcd170->Get(dir+name);
      TH1F *hqcd300 = (TH1F*)fileqcd300->Get(dir+name);
      TH1F *hqcd470 = (TH1F*)fileqcd470->Get(dir+name);
      TH1F *hqcd600 = (TH1F*)fileqcd600->Get(dir+name);
      TH1F *hqcd800 = (TH1F*)fileqcd800->Get(dir+name);
      TH1F *hqcd1000 = (TH1F*)fileqcd1000->Get(dir+name);
      format_h(hqcd15,kGray);  format_h(hqcd20,kGray);
      format_h(hqcd30,kGray);  format_h(hqcd50,kGray);
      format_h(hqcd80,kGray);  format_h(hqcd120,kGray);
      format_h(hqcd170,kGray); format_h(hqcd300,kGray);
      format_h(hqcd470,kGray); format_h(hqcd600,kGray);
      format_h(hqcd800,kGray); format_h(hqcd1000,kGray);

      TH1F *hdataB = (TH1F*)datafB->Get(dir+name);
      TH1F *hdataC = (TH1F*)datafC->Get(dir+name);
      TH1F *hdataD = (TH1F*)datafD->Get(dir+name);
      TH1F* hdata  = (TH1F*)hdataB->Clone("hdata");
      hdata->Add(hdataC) ;
      hdata->Add(hdataD) ;

    /*  overflow(htt);overflow(hwj);
      overflow(hst_twch_at);overflow(hst_twch_t);
      overflow(hst_sch);overflow(hqcd);
      overflow(hww);overflow(hzz);
      overflow(hdata);

      hst_tch_t->Rebin(2.5);hst_twch_at->Rebin(2.5);hst_twch_t->Rebin(2.5);
      hst_sch->Rebin(2.5);hww->Rebin(2.5);hzz->Rebin(2.5);
      hdata->Rebin(2.5);
    */
      hdata->SetMarkerSize(0.78);
      hdata->SetLineColor(kBlack);
      hdata->SetLineWidth(1.5);
      hdata->SetMarkerColor(kBlack);
      hdata->SetMarkerStyle(20);
      hdata->Sumw2();

      Int_t n_bins=hdata->GetSize();
      TString nam=htt->GetName();
      cout << fixed;

      if (nam=="events_counter")cout<<"float tt_evnt = "<<htt->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float ttW_evnt = "<<httW->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float ttZ_evnt = "<<httZ->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float wj_evnt = "<<hwj->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float DY_evnt = "<<hDY->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float st_tch_top_evnt="<<hst_tch_t->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float st_twch_top_evnt="<<hst_twch_at->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float st_twch_atop_evnt="<<hst_twch_t->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float st_sch_evnt="<<hst_sch->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float ww_evnt = "<<hww->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float zz_evnt = "<<hzz->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float wz_evnt = "<<hwz->GetBinContent(1)<<";"<<endl;

      if (nam=="events_counter")cout<<"float qcd15_e  =  "<<hqcd15->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float qcd20_e  =  "<<hqcd20->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float qcd30_e  =  "<<hqcd30->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float qcd50_e  =  "<<hqcd50->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float qcd80_e  =  "<<hqcd80->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float qcd120_e =  "<<hqcd120->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float qcd170_e =  "<<hqcd170->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float qcd300_e =  "<<hqcd300->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float qcd470_e =  "<<hqcd470->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float qcd600_e =  "<<hqcd600->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float qcd800_e =  "<<hqcd800->GetBinContent(1)<<";"<<endl;
      if (nam=="events_counter")cout<<"float qcd1000_e=  "<<hqcd1000->GetBinContent(1)<<";"<<endl;

//      for (Int_t i=1; i<=n_bins-2; i++){if (nam=="events_counter"){cout<<"this is events: "<<htt->GetBinContent(i)<<endl;}}


      float lumi=12877.4; //muon = 12877.401701, e = 12883.8601 
      float tt_evnt = 85676414.188449;
      float ttW_evnt = 741220.330057;
      float ttZ_evnt = 508217.976413;
      float wj_evnt = 6814531160594.228516;
      float DY_evnt = 1662002371648.298828;
      float st_tch_top_evnt=1220447540.405972;
      float st_twch_top_evnt=5477541.757331;
      float st_twch_atop_evnt=5494360.755072;
      float st_sch_evnt=5477541.757331;
      float ww_evnt = 6660781.662860;
      float zz_evnt = 1005798.714968;
      float wz_evnt = 1016101.679828;
      float qcd15_e  =  4693705.289944;
      float qcd20_e  =  32032749.277747;
      float qcd30_e  =  28809945.038684;
      float qcd50_e  =  20720147.969795;
      float qcd80_e  =  14085432.511000;
      float qcd120_e =  8056061.989751;
      float qcd170_e =  7918425.238714;
      float qcd300_e =  7609923.745327;
      float qcd470_e =  3926058.285534;
      float qcd600_e =  4027549.343103;
      float qcd800_e =  4019175.159715;
      float qcd1000_e=  3947838.760831;
       
      htt->Scale((831.76*lumi)/tt_evnt); 
      httW->Scale((0.2043*lumi)/ttW_evnt); 
      httZ->Scale((0.2529*lumi)/ttZ_evnt); 
      hwj->Scale((61526.7/wj_evnt )*lumi);    
      hDY->Scale((5765.4/DY_evnt )*lumi );    
      hst_tch_t->Scale((70.70*lumi)/st_tch_top_evnt);
      hst_twch_at->Scale((35.6*lumi*0.)/st_twch_atop_evnt);
      hst_twch_t->Scale((35.6*lumi*0.)/st_twch_top_evnt);
      hst_sch->Scale((11.36*lumi*0.32*0.0*0.)/st_sch_evnt);

      hqcd15->Scale((1273190000.*lumi*0.003)/qcd15_e); hqcd20->Scale((558528000.*lumi*0.0053)/qcd20_e);
      hqcd30->Scale((139803000.*lumi*0.01182)/qcd30_e);hqcd50->Scale((19222500.*lumi*0.02276)/qcd50_e);
      hqcd80->Scale((2758420.*lumi*0.03844)/qcd80_e);  hqcd120->Scale((469797.*lumi*0.05362)/qcd120_e);
      hqcd170->Scale((117989.*lumi*0.07335)/qcd170_e); hqcd300->Scale((7820.25*lumi*0.10196)/qcd300_e);
      hqcd470->Scale((645.528*lumi*0.12242)/qcd470_e); hqcd600->Scale((187.109*lumi*0.13412)/qcd600_e);
      hqcd800->Scale((32.3486*lumi*0.14552)/qcd800_e); hqcd1000->Scale((10.4305*lumi*0.15544)/qcd1000_e);
      TH1D *hqcd = (TH1D*)hqcd15->Clone("hqcd");
      hqcd->Add(hqcd20);hqcd->Add(hqcd30);
      hqcd->Add(hqcd50);hqcd->Add(hqcd80);
      hqcd->Add(hqcd120);hqcd->Add(hqcd170);
      hqcd->Add(hqcd300);hqcd->Add(hqcd470);
      hqcd->Add(hqcd600);
      hqcd->Add(hqcd800);hqcd->Add(hqcd1000);
      hww->Scale((118.7*lumi*0.)/ww_evnt);
      hzz->Scale((16.523*lumi*0.)/zz_evnt);
      hwz->Scale((47.13*lumi*0.)/wz_evnt);
//      hqcd->Scale(0.);
//--------------------------

TH1D * dummy = (TH1D*)hww->Clone("dummy");

  float errLumi = 0.03; // lumi uncertainty
  float errQCD = 0.15;  // uncertainty
  float errVV = 0.15;   // uncertainty on VV background norm
  float errW = 0.15;    // uncertainty on W normalization 
  float errDY = 0.15;    // uncertainty on DY normalization check for it??
  float errTT = 0.1;    // uncertainty on TTbar normalization
  float errST_tch = 0.20;    // uncertainty on ST t-channel
  float errST_tWch = 0.15;    // uncertainty on ST tW-channel
  float errMuonEff = 0.02; // muon Id/Iso/Trigger uncertainty
  float errMuonTrigEff = 0.02;// Muon Trigger uncertainty
  float errBtagEff = 0.04;// Btagging uncertainty
  float errCtagEff = 0.01;// Ctagging uncertainty
  float errPU = 0.017;// Pile Up uncertainty
  float errJER = 0.017;// JER uncertainty
  float errJES = 0.017;// JES uncertainty
  for (int iB=1; iB<=n_bins; ++iB) {
    float eQCD=errQCD*hqcd->GetBinContent(iB);
    float eVV = errVV*hww->GetBinContent(iB);
    float eW = errW*hwj->GetBinContent(iB);
    float eDY = errDY*hDY->GetBinContent(iB);
    float eTT = errTT*htt->GetBinContent(iB);
    float eST_tch = errST_tch*hst_tch_t->GetBinContent(iB);
    float eST_tWch = errST_tWch*hst_twch_t->GetBinContent(iB);
    float err2 = eQCD*eQCD + eVV*eVV + eW*eW + eDY*eDY + eTT*eTT + eST_tch*eST_tch + eST_tWch*eST_tWch;
    float errTot = TMath::Sqrt(err2);
    dummy->SetBinError(iB,errTot);
  }
     TH1* totBkg= (TH1D*)htt->Clone("totBkg");
     totBkg->Add(httW);
     totBkg->Add(httZ);
     totBkg->Add(hwj);
     totBkg->Add(hDY);
     totBkg->Add(hst_tch_t);
     totBkg->Add(hst_twch_at);
     totBkg->Add(hst_twch_t);
//     totBkg->Add(hst_sch);
     totBkg->Add(hqcd);
     totBkg->Add(hww);
     totBkg->Add(hzz);
     totBkg->Add(hwz);
  TH1D * bkgdErr = (TH1D*)totBkg->Clone("bkgdErr");
  bkgdErr->SetFillStyle(3013);
  bkgdErr->SetFillColor(1);
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);  
  
  for (int iB=1; iB<=n_bins; ++iB) {
    hqcd->SetBinError(iB,0);
    hww->SetBinError(iB,0);
    htt->SetBinError(iB,0);
    hwj->SetBinError(iB,0);
    hDY->SetBinError(iB,0);
    hzz->SetBinError(iB,0);
    hww->SetBinError(iB,0);
    float eStat =  bkgdErr->GetBinError(iB);
    float X = bkgdErr->GetBinContent(iB);
    float eLumi = errLumi * X;
    float eMuonEff=errMuonEff *X;
    float eMuonTrigEff=errMuonTrigEff *X;
    float eBtagEff = errBtagEff * X;
    float eCtagEff = errCtagEff * X;
    float ePU      = errPU      * X;
    float eJER     = errJER     * X;
    float eJES     = errJES     * X;
    float eBkg = dummy->GetBinError(iB);
    float Err = TMath::Sqrt(eStat*eStat+eLumi*eLumi+eMuonEff*eMuonEff+eMuonTrigEff*eMuonTrigEff+
                            +eBtagEff*eBtagEff+eCtagEff*eCtagEff+ePU*ePU+eJER*eJER+eJES*eJES+eBkg*eBkg);
    bkgdErr->SetBinError(iB,Err);
  }

/*     double norm_data=hdata->Integral()/totBkg->Integral();
     //std::cout<<norm_data<<endl;
      htt->Scale(norm_data);
      hwj->Scale(norm_data);
      hst_tch_t->Scale(norm_data);
      hst_twch_at->Scale(norm_data);
      hst_twch_t->Scale(norm_data);
      hst_sch->Scale(norm_data);
//      hqcd->Scale(norm_data);
      hww->Scale(norm_data);
      hzz->Scale(norm_data);
*/
 //     hs->Add(hqcd);
      hs->Add(hww);
      hs->Add(hzz);
      hs->Add(hwz);
      hs->Add(hst_twch_t);
      hs->Add(hst_twch_at);
      hs->Add(hst_tch_t);
//      hs->Add(hst_sch);
      hs->Add(hDY);
      hs->Add(hwj);
      hs->Add(httZ);
      hs->Add(httW);
      hs->Add(htt);
      hs->Add(hqcd);
     
     c1->cd();
     TPad*    upperPad = new TPad("upperPad", "upperPad",0.043, 0.32, 0.90, 0.95);
     upperPad->SetLogy();
     upperPad->SetTopMargin(0);
     upperPad->SetBottomMargin(0);
     upperPad->SetTickx(1);
     upperPad->SetTicky(1);
     upperPad->Draw();
     upperPad->cd();
     hs->Draw("HIST");
//     hM400->Draw("P,SAMES");
//     hM600->Draw("P,SAMES");
//     hM800->Draw("P,SAME,L");

     hdata->Draw("P,SAMES,E");
//     bkgdErr->Draw("e2same");
     upperPad->Modified();
     upperPad->cd();
      hs->SetMaximum(hdata->GetMaximum()*2.);
      hs->GetXaxis()->SetTitle("");
      hs->GetXaxis()->SetTitleSize(0);
      hs->GetXaxis()->SetLabelSize(0);
      hs->GetXaxis()->SetTitleOffset(2.035);
	    hs->GetXaxis()->SetTitleOffset(2.80);
      Int_t n_bin=hdata->GetSize()-2;
      float maxx=hdata->GetXaxis()->GetXmax();
      float minn=hdata->GetXaxis()->GetXmin();
      Int_t n_bin=hdata->GetNbinsX();
      float bin_size=(maxx-minn)/n_bin;
      cout<<"this is no of bins "<<n_bin<<"     "<<maxx<<"     "<<minn<<"      "<<maxx-minn<<"      "<<bin_size<<endl;
      ostringstream temp;
      temp<<bin_size;
      string str=temp.str();
      string evnt="Events/"+str;
      hs->GetYaxis()->SetTitle(evnt.c_str());
      hs->GetYaxis()->SetTitleSize(0.07);
      hs->GetYaxis()->SetLabelSize(0.065);
      hs->GetYaxis()->SetTitleOffset(0.80);

     // increase arg 1 move to right, arg2 start the lower edge, increase arg3 move to left, arg 4 is the upper limit,
     c1->cd();
     TPad*    lowerPad = new TPad("lowerPad", "lowerPad",0.043, 0.0985, 0.90,0.318);
 //    dividend->GetYaxis()->SetRangeUser(2.0,0.0);
     lowerPad->SetTopMargin(0);
     lowerPad->SetBottomMargin(0);
     lowerPad->SetGridy();
     lowerPad->SetTickx(1);
     lowerPad->SetTicky(1);
     lowerPad->Draw();
     lowerPad->cd();
     
     TH1F* dividend=new TH1F(*hdata);
     TH1F* totalMC=(TH1F*)htt->Clone();
     totalMC->Add(httW);
     totalMC->Add(httZ);
     totalMC->Add(hwj);
     totalMC->Add(hDY);
     totalMC->Add(hst_tch_t);
     totalMC->Add(hst_twch_at);
     totalMC->Add(hst_twch_t);
     totalMC->Add(hst_sch);
     totalMC->Add(hwz);
     totalMC->Add(hzz);
     totalMC->Add(hww);
     totalMC->Add(hqcd);

     TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
    for (int iB=1; iB<=n_bins; ++iB) {
    float x1 = hdata->GetBinContent(iB);
    float x2 = totalMC->GetBinContent(iB);
    ratioErrH->SetBinContent(iB,1.0);
    ratioErrH->SetBinError(iB,0.0);
    float xBkg = bkgdErr->GetBinContent(iB);
    float errBkg = bkgdErr->GetBinError(iB);
    if (xBkg>0) {
      float relErr = errBkg/xBkg;
      ratioErrH->SetBinError(iB,relErr);
    }
    }

     dividend->Divide(totalMC);
//     hdata->GetXaxis()->SetLabelSize(0.0);

     bin_name(dividend);
     dividend->GetYaxis()->SetTitle("Data/MC");
     dividend->GetYaxis()->CenterTitle(true);
//     dividend->GetXaxis()->SetTitle(xtitle);
     dividend->GetXaxis()->SetLabelSize(0.17);
// dividend->GetXaxis()->SetLabelOffset(4.02);
     dividend->GetXaxis()->SetTitleSize(0.17);
     dividend->GetXaxis()->SetTitleOffset(0.88);
     dividend->GetYaxis()->SetLabelSize(0.14);
     dividend->GetYaxis()->SetTitleSize(0.20);
     dividend->GetYaxis()->SetTitleOffset(0.27);
     dividend->GetYaxis()->SetNdivisions(5);
   //  dividend->GetYaxis()->SetRange(0.0,2.0);
//     float max=dividend->SetMaximum(dividend->GetMaximum()*1.1);   
    // float min=dividend->SetMinimum(dividend->GetMinimum()*1.3);
    dividend->SetMinimum(0.1);
   dividend->SetMaximum(1.9);


     dividend->SetStats(0);
     format_h(dividend,kBlack);
     gStyle->SetOptStat(0);
     dividend->Draw("e1");
     dividend->SetTitle("");
//     ratioErrH->Draw("e2same");

//TPaveText *t = gPad->FindObject("title");
//t->AddText("something");
        upperPad->cd();

    TLatex *tex1 = new TLatex(0.10,1.02,"#bf{CMS} #it{Preliminary}");
    tex1->SetNDC();
    tex1->SetTextAngle(0);
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.06);
    tex1->SetTextAlign(11);
    tex1->Draw();
    TLatex *tex2 = new TLatex(0.65,1.02,"12.9 fb^{-1} (13 TeV)");
    tex2->SetNDC();
    tex2->SetTextAngle(0);
    tex2->SetTextFont(42);
    tex2->SetTextAlign(11);
    tex2->SetTextSize(0.06);
    tex2->Draw();
    TLatex *tex3 = new TLatex(0.4,1.02,"#mu + jets");
    tex3->SetNDC();
    tex3->SetTextAngle(0);
    tex3->SetTextFont(42);
    tex3->SetTextAlign(11);
    tex3->SetTextSize(0.06);
    tex3->Draw();


//        TLegend *leg = new TLegend(0.3,0.98,0.85,0.77,NULL,"brNDC");//1)move to left and expand 2)move down and vertically seperate 3)move left and shrink 4) shrink vertically
        TLegend *leg = new TLegend(0.91,1.02,1.099,0.4,NULL,"brNDC");//1)move to left and expand 2)move down and vertically seperate 3)move left and shrink 4) shrink vertically
//        leg-> SetNColumns(3);
        leg->SetFillStyle ( 0);
        leg->SetFillColor ( 0);
        leg->SetBorderSize( 0);
        leg->SetLineWidth(1);
        leg->SetTextSize(0.05);
        leg->SetFillColor(0);
        leg->Draw();
        leg->AddEntry(hdata,"Data","LP");
        leg->AddEntry(hqcd,"QCD","F");
        leg->AddEntry(htt,"t#bar{t}","f");
        leg->AddEntry(httW,"t#bar{t}V","f");
        leg->AddEntry(hwj,"wjets","f");
        leg->AddEntry(hDY,"Z/#gamma","F");
        leg->AddEntry(hst_tch_t,"s-top","F");
//        leg->AddEntry(hzz,"VV","F");
//        leg->AddEntry(hqcd,"QCD","F");

 //       leg->AddEntry(hM400,"HM-400(*1600)","LP");
  //      leg->AddEntry(hM600,"HM-600(*800)","LP");
  //      leg->AddEntry(hM800,"HM-800(*800)","LP");

//        leg->AddEntry(hqcd,"multijets","F");
	
    }



void FinalHistos()
{
 plot_shape("events_counter","mujets_2btag/","events_counter","No. of Events");
 plot_shape("events_eachCut","mujets_2btag/","events_eachCut","No. of Events");
        c1->SaveAs("results/events_eachCut.png");

 plot_shape("jets_Pt_afterCut","mujets_2btag/","jets_Pt_afterCut","No. of Events");
        c1->SaveAs("results/jets_Pt_afterCut.png");
 plot_shape("jets_mass_afterCut","mujets_2btag/","jets_mass_afterCut","No. of Events");
        c1->SaveAs("results/jets_mass_afterCut.png");
 plot_shape("jets_energy_afterCut","mujets_2btag/","jets_energy_afterCut","No. of Events");
        c1->SaveAs("results/jets_energy_afterCut.png");
 plot_shape("jets_lep_Pt_afterCut","mujets_2btag/","jets_lep_Pt_afterCut","No. of Events");
        c1->SaveAs("results/jets_lep_Pt_afterCut.png");

plot_shape("MET_Px","mujets_2btag/","MET_Px","No. of Events");
        c1->SaveAs("results/MET_Px.png");
plot_shape("MET_Py","mujets_2btag/","MET_Py","No. of Events");
        c1->SaveAs("results/MET_Py.png");
plot_shape("MET_Pz","mujets_2btag/","MET_Pz","No. of Events");
        c1->SaveAs("results/MET_Pz.png");
plot_shape("MET_E","mujets_2btag/","MET [GeV]","No. of Events");
        c1->SaveAs("results/MET_E.png");
plot_shape("patMuonPfIsoDbeta_after","mujets_2btag/","IsoDbeta_","No. of Events");
        c1->SaveAs("results/IsoDbeta_.png");
plot_shape("patMuonPfIsoDbeta","mujets_2btag/","IsoDbeta","No. of Events");
        c1->SaveAs("results/IsoDbeta.png");
plot_shape("Pt_lep","mujets_2btag/","P_{T} Lep [GeV]","No. of Events");
        c1->SaveAs("results/Pt_lep.png");
plot_shape("Eta_lep","mujets_2btag/","Eta_lep","No. of Events");
        c1->SaveAs("results/Eta_lep.png");
plot_shape("Phi_lep","mujets_2btag/","Phi_lep","No. of Events");
        c1->SaveAs("results/Phi_lep.png");
plot_shape("E_lep","mujets_2btag/","E_lep","No. of Events");
        c1->SaveAs("results/E_lep.png");
plot_shape("Pt_LepW","mujets_2btag/","Pt_LepW","No. of Events");
        c1->SaveAs("results/Pt_LepW.png");
plot_shape("M_LepW","mujets_2btag/","Mass_LepW","No. of Events");
        c1->SaveAs("results/mass_LepW.png");
plot_shape("M_T","mujets_2btag/","M_{T} [GeV]","No. of Events");
        c1->SaveAs("results/M_T.png");
plot_shape("M_T_aft","mujets_2btag/","M_{T} [GeV]","No. of Events");
        c1->SaveAs("results/M_T_aft.png");

plot_shape("n_Jets","mujets_2btag/","n Jets Pt>30","No. of Events");
        c1->SaveAs("results/n_Jets30.png");
plot_shape("n_Jets20","mujets_2btag/","n Jets Pt>20","No. of Events");
        c1->SaveAs("results/n_Jets20.png");
plot_shape("PtJet1","mujets_2btag/","P_{T} Jet1 [GeV]","No. of Events");
        c1->SaveAs("results/PtJet1.png");
plot_shape("EtaJet1","mujets_2btag/","#eta Jet1","No. of Events");
        c1->SaveAs("results/EtaJet1.png");
plot_shape("PhiJet1","mujets_2btag/","#phi Jet1","No. of Events");
        c1->SaveAs("results/PhiJet1.png");
plot_shape("EJet1","mujets_2btag/","EJet1","No. of Events");
        c1->SaveAs("results/EJet1.png");
plot_shape("MJet1","mujets_2btag/","MJet1","No. of Events");
        c1->SaveAs("results/MJet1.png");

plot_shape("PtJet2","mujets_2btag/","P_{T} Jet2 [GeV]","No. of Events");
        c1->SaveAs("results/PtJet2.png");
plot_shape("EtaJet2","mujets_2btag/","EtaJet2","No. of Events");
        c1->SaveAs("results/EtaJet2.png");
plot_shape("PhiJet2","mujets_2btag/","PhiJet2","No. of Events");
        c1->SaveAs("results/PhiJet2.png");
plot_shape("EJet2","mujets_2btag/","EJet2","No. of Events");
        c1->SaveAs("results/EJet2.png");
plot_shape("MJet2","mujets_2btag/","MJet2","No. of Events");
        c1->SaveAs("results/MJet2.png");

plot_shape("PtJet3","mujets_2btag/","P_{T} Jet3 [GeV]","No. of Events");
        c1->SaveAs("results/PtJet3.png");
plot_shape("EtaJet3","mujets_2btag/","EtaJet3","No. of Events");
        c1->SaveAs("results/EtaJet3.png");
plot_shape("PhiJet3","mujets_2btag/","PhiJet3","No. of Events");
        c1->SaveAs("results/PhiJet3.png");
plot_shape("EJet3","mujets_2btag/","EJet3","No. of Events");
        c1->SaveAs("results/EJet3.png");
plot_shape("MJet3","mujets_2btag/","MJet3","No. of Events");
        c1->SaveAs("results/MJet3.png");

plot_shape("PtJet4","mujets_2btag/","P_{T} Jet4 [GeV]","No. of Events");
        c1->SaveAs("results/PtJet4.png");
plot_shape("EtaJet4","mujets_2btag/","EtaJet4","No. of Events");
        c1->SaveAs("results/EtaJet4.png");
plot_shape("PhiJet4","mujets_2btag/","PhiJet4","No. of Events");
        c1->SaveAs("results/PhiJet4.png");
plot_shape("EJet4","mujets_2btag/","EJet4","No. of Events");
        c1->SaveAs("results/EJet4.png");
plot_shape("MJet4","mujets_2btag/","MJet4","No. of Events");
        c1->SaveAs("results/MJet4.png");


plot_shape("n_bjets","mujets_2btag/","n_bjets","No. of Events");
        c1->SaveAs("results/n_bJets.png");
plot_shape("n_bjets_corr","mujets_2btag/","n_bjets_corr","No. of Events");
        c1->SaveAs("results/n_bJets_corr.png");
plot_shape("n_bjets_aft","mujets_2btag/","n_bjets_aft","No. of Events");
        c1->SaveAs("results/n_bJets_aft.png");
plot_shape("Pt_bJet1","mujets_2btag/","P_{T} bJet1 [GeV]","No. of Events");
        c1->SaveAs("results/Pt_bJet1.png");
plot_shape("Eta_bJet1","mujets_2btag/","Eta_bJet1","No. of Events");
        c1->SaveAs("results/Eta_bJet1.png");
plot_shape("Phi_bJet1","mujets_2btag/","Phi_bJet1","No. of Events");
        c1->SaveAs("results/Phi_bJet1.png");
plot_shape("E_bJet1","mujets_2btag/","E_bJet1","No. of Events");
        c1->SaveAs("results/E_bJet1.png");
plot_shape("M_bJet1","mujets_2btag/","M_bJet1","No. of Events");
        c1->SaveAs("results/M_bJet1.png");

plot_shape("Pt_bJet2","mujets_2btag/","P_{T} bJet2 [GeV]","No. of Events");
        c1->SaveAs("results/Pt_bJet2.png");
plot_shape("Eta_bJet2","mujets_2btag/","Eta_bJet2","No. of Events");
        c1->SaveAs("results/Eta_bJet2.png");
plot_shape("Phi_bJet2","mujets_2btag/","Phi_bJet2","No. of Events");
        c1->SaveAs("results/Phi_bJet2.png");
plot_shape("E_bJet2","mujets_2btag/","E_bJet2","No. of Events");
        c1->SaveAs("results/E_bJet2.png");
plot_shape("M_bJet2","mujets_2btag/","M_bJet2","No. of Events");
        c1->SaveAs("results/M_bJet2.png");
plot_shape("DPhi_bj12","mujets_2btag/","DPhi_bjet12","No. of Events");
        c1->SaveAs("results/DPhi_bjet12.png");
plot_shape("DR_bj12","mujets_2btag/","DR_bjet12","No. of Events");
        c1->SaveAs("results/DR_bjet12.png");

plot_shape("n_ljets","mujets_2btag/","n_ljets","No. of Events");
        c1->SaveAs("results/n_lJets.png");
plot_shape("n_ljets_aft","mujets_2btag/","n_ljets_aft","No. of Events");
        c1->SaveAs("results/n_lJets_aft.png");
plot_shape("Pt_lJet1","mujets_2btag/","P_{T} lJet1 [GeV]","No. of Events");
        c1->SaveAs("results/Pt_lJet1.png");
plot_shape("Eta_lJet1","mujets_2btag/","Eta_lJet1","No. of Events");
        c1->SaveAs("results/Eta_lJet1.png");
plot_shape("Phi_lJet1","mujets_2btag/","Phi_lJet1","No. of Events");
        c1->SaveAs("results/Phi_lJet1.png");
plot_shape("E_lJet1","mujets_2btag/","E_lJet1","No. of Events");
        c1->SaveAs("results/E_lJet1.png");
plot_shape("M_lJet1","mujets_2btag/","M_lJet1","No. of Events");
        c1->SaveAs("results/M_lJet1.png");

plot_shape("Pt_lJet2","mujets_2btag/","P_{T} lJet2 [GeV]","No. of Events");
        c1->SaveAs("results/Pt_lJet2.png");
plot_shape("Eta_lJet2","mujets_2btag/","Eta_lJet2","No. of Events");
        c1->SaveAs("results/Eta_lJet2.png");
plot_shape("Phi_lJet2","mujets_2btag/","Phi_lJet2","No. of Events");
        c1->SaveAs("results/Phi_lJet2.png");
plot_shape("E_lJet2","mujets_2btag/","E_lJet2","No. of Events");
        c1->SaveAs("results/E_lJet2.png");
plot_shape("M_lJet2","mujets_2btag/","M_lJet2","No. of Events");
        c1->SaveAs("results/M_lJet2.png");

plot_shape("Pt_hadW","mujets_2btag/","P_{T} Had W [GeV]","No. of Events");
        c1->SaveAs("results/Pt_hadW.png");
plot_shape("M_hadW","mujets_2btag/","Mass Had W [GeV]","No. of Events");
        c1->SaveAs("results/M_hadW.png");

plot_shape("M_lepTop","mujets_2btag/","Mass Lep Top [GeV]","No. of Events");
        c1->SaveAs("results/M_lepTop.png");
plot_shape("Pt_lepTop","mujets_2btag/","P_{T} Lep Top [GeV]","No. of Events");
        c1->SaveAs("results/Pt_lepTop.png");

plot_shape("M_hadTop","mujets_2btag/","Mass Had Top [GeV]","No. of Events");
        c1->SaveAs("results/M_hadTop.png");
plot_shape("Pt_hadTop","mujets_2btag/","P_{T} Had Top [GeV]","No. of Events");
        c1->SaveAs("results/Pt_hadTop.png");

plot_shape("Mass_H_binned","mujets_2btag/","Mass t#bar{t} [GeV]","No. of Events");
        c1->SaveAs("results/Mass_H_binned.png");
plot_shape("Mass_H","mujets_2btag/","Mass t#bar{t} [GeV]","No. of Events");
        c1->SaveAs("results/Mass_H.png");
plot_shape("Pt_H","mujets_2btag/","P_{T} t#bar{t} [GeV]","No. of Events");
        c1->SaveAs("results/Pt_ttbar.png");
plot_shape("cos_theta","mujets_2btag/","cos#theta^{*} [GeV]","No. of Events");
        c1->SaveAs("results/cos_theta.png");

plot_shape("EvtInfo_NumVtx","mujets_2btag/","EvtInfo_NumVtx","No. of Events");
        c1->SaveAs("results/EvtInfo_NumVtx.png");
plot_shape("EvtInfo_NumVtx_w","mujets_2btag/","EvtInfo_NumVtx_w","No. of Events");
        c1->SaveAs("results/EvtInfo_NumVtx_w.png");
plot_shape("PU_npT","mujets_2btag/","PU_npT","No. of Events");
        c1->SaveAs("results/PU_npT.png");
plot_shape("bDiscMVAv2","mujets_2btag/","bDiscMVAv2","No. of Events");
        c1->SaveAs("results/bDiscMVAv2.png");
plot_shape("cMVA_v2","mujets_2btag/","cMVA_v2","No. of Events");
        c1->SaveAs("results/cMVAafter_v2.png");
plot_shape("Jet_cMult","mujets_2btag/","Jet_cMult","No. of Events");
        c1->SaveAs("results/Jet_cMult.png");


}



