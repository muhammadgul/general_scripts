#include <iostream>
#include <fstream>
/*TString paths1 = "/afs/cern.ch/work/m/mgul/public/Hto_ttbar/CMSSW_7_2_0/src/Tupel/Tupel/tupelanalyzer/muon_analyzer/muonHM400_allevents/";
TFile *file1 = new TFile(paths1+"HM400to_ttbar_2.root");
TString paths2= "/afs/cern.ch/work/m/mgul/public/Hto_ttbar/CMSSW_7_2_0/src/Tupel/Tupel/tupelanalyzer/muon_analyzer/muonHM600_allevents/";
TFile *file2= new TFile(paths2+"HM600to_ttbar_2.root");
TString paths3= "/afs/cern.ch/work/m/mgul/public/Hto_ttbar/CMSSW_7_2_0/src/Tupel/Tupel/tupelanalyzer/muon_analyzer/muonHM800_allevents/";
TFile *file3= new TFile(paths3+"HM800to_ttbar_2.root");*/
TString pathtt= "/home/muhammad/root_files/";
TFile *ttfile= new TFile(pathtt+"Ttbar_bkg.root");

TString pathwj= "/home/muhammad/root_files/";
TFile *wjfile= new TFile(pathwj+"wjets_bkg.root");

TString path6= "/home/muhammad/root_files/";
TFile *file6= new TFile(path6+"single_top1.root");TFile *file7= new TFile(path6+"single_top2.root");TFile *file8= new TFile(path6+"single_top3.root");
TFile *file9= new TFile(path6+"single_top4.root");TFile *file10= new TFile(path6+"single_top5.root");

TString pathqcd= "/home/muhammad/root_files/";
TFile *fileqcd= new TFile(pathqcd+"qcdPt_20toInf.root");
TString pathww= "/home/muhammad/root_files/";
TFile *fileww= new TFile(pathww+"ww_bkg.root");
TString pathzz= "/home/muhammad/root_files/";
TFile *filezz= new TFile(pathzz+"zz_bkg.root");

TString path10= "/home/muhammad/root_files/";
TFile *datafile= new TFile(path10+"data_singleMu.root");
//void plot_shape(TString name,TString dir, TString xtitle, TString ytitle, int binning = 1 )
void format_h(TH1F* h, int fillcolor){
     h->SetLineWidth(1.5);
     h->SetLineColor(fillcolor);
     h->SetMarkerStyle(20);
     h->SetFillColor(fillcolor);
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
      {"No Cut","HLT","elec=0","mu=1","DR(mu,j)","Jets >= 4","bjets >=2","ljets >=2"};
      TString nam=h2->GetName();
      if (nam=="events_eachCut"){
      for (int i=0;i<sizeof(naming)/sizeof(TString);i++){
      h2->GetXaxis()->SetBinLabel(i+1,naming[i]);
      }}}

void plot_shape(TString name, TString xtitle, TString ytitle, int binning = 1 )
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
   //c1->SetLogy(kTRUE);
   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();
   //gStyle->SetOptStat(1111111);
  /*      TH1F *hs1 = (TH1F*)file1->Get(name);
  format_h(hs1,kRed+2);
  TH1F *hs2 = (TH1F*)file2->Get(name);
  format_h(hs2, kBlue);
  TH1F *hs3 = (TH1F*)file3->Get(name);
  format_h(hs3, kRed);*/
  TH1F *htt = (TH1F*)ttfile->Get(name);
  format_h(htt,kGreen);
  TH1F *hwj = (TH1F*)wjfile->Get(name);
  format_h(hwj,kOrange);
  TH1F *hst_tch_at = (TH1F*)file6->Get(name);
  format_h(hst_tch_at,kViolet);
  TH1F *hst_tch_t = (TH1F*)file7->Get(name);
  format_h(hst_tch_t,kViolet);
  TH1F *hst_twch_at = (TH1F*)file8->Get(name);
  format_h(hst_twch_at,kViolet);
  TH1F *hst_twch_t = (TH1F*)file9->Get(name);
  format_h(hst_twch_t,kViolet);
  TH1F *hst_sch = (TH1F*)file10->Get(name);
    format_h(hst_sch,kViolet);

  TH1F *hqcd = (TH1F*)fileqcd->Get(name);
  format_h(hqcd,kBlue);
  TH1F *hww = (TH1F*)fileww->Get(name);
  format_h(hww,kYellow);
  TH1F *hzz = (TH1F*)filezz->Get(name);
  format_h(hzz,kCyan);

  TH1F *hdata = (TH1F*)datafile->Get(name);

  hdata->SetMarkerSize(0.78);
  hdata->SetLineColor(kBlack);
  hdata->SetLineWidth(1.5);
  hdata->SetMarkerColor(kBlack);
  hdata->SetMarkerStyle(20);
	hdata->Sumw2();
  //float lumi=41.0; int s_evnt=200000, tt_evnt=200000,wj_evnt=200000, st_tch_top_evnt=200000,st_tch_atop_evnt=200000,st_twch_top_evnt= 200000,st_twch_atop_evnt=200000;
  float lumi=552.67; 
  int s_evnt=200000; 
  int tt_evnt=19899492;
  int wj_evnt=24151209;
  int st_tch_top_evnt=3299799;
  int st_tch_atop_evnt=1695400;
  int st_twch_top_evnt= 995598;
  int st_twch_atop_evnt=1000000;
  int st_sch_evnt=984399;
  int qcd_evnt=13201689;
  int zz_evnt=996168;
  int ww_evnt=994416;

  TString nam=hdata->GetName();
  Int_t n_bins=hdata->GetSize();
  for (Int_t i=1; i<=n_bins-2; i++){
    if(nam=="events_eachCut"){
      cout<<"this is hqcd  "<<htt->GetBinContent(i)<<endl;
      }
      }

      htt->Scale((831.76*lumi)/tt_evnt); 
      hwj->Scale((61526.7*lumi)/wj_evnt);    
      hst_tch_at->Scale((80.95*lumi)/st_tch_atop_evnt);  
      hst_tch_t->Scale((136.02*lumi)/st_tch_top_evnt);
      hst_twch_at->Scale((35.6*lumi)/st_twch_atop_evnt);
      hst_twch_t->Scale((35.6*lumi)/st_twch_top_evnt);
      hst_sch->Scale((11.36*lumi)/st_sch_evnt);
      hqcd->Scale((381304*lumi)/qcd_evnt);
      hww->Scale((118.7*lumi)/ww_evnt);
      hzz->Scale((15.4*lumi)/zz_evnt);

      ofstream outputFile;
      outputFile.open("cutFlowTable.tex");

      outputFile<<"\\documentclass[english,a4paper,12pt]{report}
      \\usepackage[T1]{fontenc}
      \\usepackage{graphicx}
      \\usepackage[latin9]{inputenc}
      \\usepackage{graphicx}
      \\usepackage[showframe=false]{geometry}
      \\usepackage{changepage}
      \\usepackage{hhline}
      \\usepackage{booktabs}
      \\usepackage[table]{xcolor}
      \\usepackage[font=normalsize,labelfont=bf,justification=centering]{caption}
      \\begin{document}
      \\begin{table}[htb]
      \\rowcolors{1}{green}{pink}
      \\begin{adjustwidth}{-1.5cm}{}
      \\resizebox{0.6\\textwidth}{!}{
      \\begin{minipage}{\\linewidth}
      \\def\\arraystretch{2.2}%

      \\begin{tabular}{ *{13}{c}}
      \\toprule
      \\multicolumn{13}{c}{\\textbf{Cut Flow Table}}"<<"\\"<<"\\"<<endl;

      outputFile<<"\\textbf{Data Sets}&\\textbf{No Cut}&\\textbf{HLT}&      \\textbf{elec=0}&\\textbf{mu=1} &\\textbf{ Pt,Eta mu}&      \\textbf{ID muon}&\\textbf{muon Iso} &\\textbf{ Jet ID}&      \\textbf{DR(mu,j)}& \\textbf{jets $\\geq$ 4}&      \\textbf{ bjets $\\geq$ 2}& \\textbf{ljets $\\geq$ 2}"<< "\\"<<"\\"<< "\\midrule"<<endl;
      TString nam=hdata->GetName();
      Int_t n_bins=hdata->GetSize();
      outputFile<<" \n \\textbf{$T\\overline{T}$}"<<endl;
      for (Int_t i=1; i<=n_bins-2; i++){if (nam=="events_eachCut"){outputFile<<"&"<<htt->GetBinContent(i);}}
      outputFile<<"\\"<<"\\ \n \\textbf{ W+J}"<<endl;
      for (Int_t i=1; i<=n_bins-2; i++){if (nam=="events_eachCut"){outputFile<<"&"<<hwj->GetBinContent(i);}}
      outputFile<<"\\"<<"\\  \n \\textbf{S-T}"<<endl;
      for (Int_t i=1; i<=n_bins-2; i++){if (nam=="events_eachCut"){
                          outputFile<<"&"<<hst_twch_t->GetBinContent(i)+
                            hst_twch_at->GetBinContent(i)+
                              hst_tch_t->GetBinContent(i)+
                                hst_tch_at->GetBinContent(i)+
                                 hst_sch->GetBinContent(i);}
                                  }
      outputFile<<"\\"<<"\\ \n \\textbf{ZZ} "<<endl;
      for (Int_t i=1; i<=n_bins-2; i++){if (nam=="events_eachCut"){outputFile<<"&"<<hzz->GetBinContent(i);}}
      outputFile<<"\\"<<"\\ \n \\textbf{WW} "<<endl;
      for (Int_t i=1; i<=n_bins-2; i++){if (nam=="events_eachCut"){outputFile<<"&"<<hww->GetBinContent(i);}}
      outputFile<<"\\"<<"\\ \n \\textbf{QCD}"<<endl;
      for (Int_t i=1; i<=n_bins-2; i++){if (nam=="events_eachCut"){outputFile<<"&"<<hqcd->GetBinContent(i);}}
      outputFile<<"\\"<<"\\ \n \\textbf{Total bkg}"<<endl;
      for (Int_t i=1; i<=n_bins-2; i++){if (nam=="events_eachCut"){
                  outputFile<<"&"<<htt->GetBinContent(i)+hwj->GetBinContent(i)
                  +hst_twch_t->GetBinContent(i)+hst_twch_at->GetBinContent(i)
                  +hst_tch_t->GetBinContent(i)+hst_tch_at->GetBinContent(i)
                  +hzz->GetBinContent(i)+hww->GetBinContent(i)+hqcd->GetBinContent(i);}}
                  outputFile<<"\\"<<"\\ \n \\textbf{Data}"<<endl;
      for (Int_t i=1; i<=n_bins-2; i++){if (nam=="events_eachCut"){outputFile<<"&"<<hdata->GetBinContent(i);}}outputFile<<"\\"<<"\\ \n "<<endl;
      outputFile<<"\\bottomrule
                   \\end{tabular}
                   \\end{minipage}}
                   \\end{adjustwidth}
                   \\end{table}
                   \\end{document}"<<endl;
                    outputFile.close();
      

      TH1F *hs2_o=htt;
//      hdata->Scale();
      TH1* stackmc= (TH1D*)hs2_o->Clone("stackmc");
 /*     stackmc->Add(hs1);
        stackmc->Add(hs2);
        stackmc->Add(hs3);*/
     stackmc->Add(htt);
     stackmc->Add(hwj);
     stackmc->Add(hst_tch_at);
     stackmc->Add(hst_tch_t);
     stackmc->Add(hst_twch_at);
     stackmc->Add(hst_twch_t);
     stackmc->Add(hst_sch);
     stackmc->Add(hqcd);
     stackmc->Add(hww);
     stackmc->Add(hzz);

     double d_all=hdata->Integral()/stackmc->Integral();
     std::cout<<d_all<<endl;
//stackmc->Scale(hdata->Integral()/stackmc->Integral());
/*    hs1->Scale(d_all);
      hs2_o->Scale(d_all);
      hs2->Scale(d_all);
      hs3->Scale(d_all);*/
/*      htt->Scale(d_all);
      hwj->Scale(d_all);
      hst_tch_at->Scale(d_all);
      hst_tch_t->Scale(d_all);
      hst_twch_at->Scale(d_all);
      hst_twch_t->Scale(d_all);
      hqcd->Scale(d_all);
      hww->Scale(d_all);
      hzz->Scale(d_all);
*/
//hs->Add(stackmc);
      hs->Add(hqcd);
      hs->Add(hww);
      hs->Add(hzz);
      hs->Add(hst_twch_t);
      hs->Add(hst_twch_at);
      hs->Add(hst_tch_t);
      hs->Add(hst_tch_at);
      hs->Add(hst_sch);
      hs->Add(hwj);
      hs->Add(htt);

     c1->cd();
     TPad*    upperPad = new TPad("upperPad", "upperPad",0.043, 0.32, 0.97, 0.95);
     upperPad->SetLogy();
     upperPad->SetTopMargin(0);
     upperPad->SetBottomMargin(0);
 //    upperPad->SetGridx();                                                                                                                             
     upperPad->SetTickx(1);
     upperPad->SetTicky(1);
     upperPad->Draw();
     upperPad->cd();
     //  hs->Draw("HIST");
//     hdata->Draw("P,E");
//     hs->Draw("SAMES HIST");
     hs->Draw("HIST");
     hdata->Draw("P,SAMES,E");
//     hdata->SetTitle(" ");
     upperPad->Modified();
     upperPad->cd();
//      hs->SetTitle(xtitle);
hs->SetMaximum(hdata->GetMaximum()*1.6);
      hs->GetXaxis()->SetTitle("");
      hs->GetXaxis()->SetTitleSize(0);
      hs->GetXaxis()->SetLabelSize(0);
      hs->GetXaxis()->SetTitleOffset(2.035);
	hs->GetXaxis()->SetTitleOffset(2.80);
      hs->GetYaxis()->SetTitle(ytitle);
      hs->GetYaxis()->SetTitleSize(0.07);
      hs->GetYaxis()->SetLabelSize(0.05);
      hs->GetYaxis()->SetTitleOffset(0.80);

     // increase arg 1 move to right, arg2 start the lower edge, increase arg3 move to left, arg 4 is the upper limit,
     c1->cd();
     TPad*    lowerPad = new TPad("lowerPad", "lowerPad",0.043, 0.085, 0.97,0.32);
     // dividend->GetYaxis()->SetRangeUser(0.42,3.40);
     lowerPad->SetTopMargin(0);
     lowerPad->SetBottomMargin(0);
     lowerPad->SetGridy();
     lowerPad->SetTickx(1);
     lowerPad->SetTicky(1);
     lowerPad->Draw();
     lowerPad->cd();
     TH1F* dividend=new TH1F(*hdata);
     TH1F* totalMC=(TH1F*)htt->Clone();
     totalMC->Add(htt);totalMC->Add(hwj);totalMC->Add(hst_tch_at);
     totalMC->Add(hst_tch_t);totalMC->Add(hst_twch_at);totalMC->Add(hst_twch_t);
     totalMC->Add(hst_sch);
     totalMC->Add(hzz);totalMC->Add(hww);totalMC->Add(hqcd);
     dividend->Divide(totalMC);
//     hdata->GetXaxis()->SetLabelSize(0.0);

     bin_name(dividend);
     dividend->GetYaxis()->SetTitle("Data/MC");
     dividend->GetYaxis()->CenterTitle(true);
   //  dividend->GetXaxis()->SetTitle(xtitle);
     
     dividend->GetXaxis()->SetTitle(xtitle);
//dividend->GetXaxis()->SetTitle("t");
     dividend->GetXaxis()->SetLabelSize(0.17);
     dividend->GetXaxis()->SetLabelOffset(0.02);
     dividend->GetXaxis()->SetTitleSize(0.20);
     dividend->GetXaxis()->SetTitleOffset(10.730);
     dividend->GetYaxis()->SetLabelSize(0.17);
     dividend->GetYaxis()->SetTitleSize(0.20);
     dividend->GetYaxis()->SetTitleOffset(0.27);
     dividend->GetYaxis()->SetNdivisions(5);
     float max=dividend->SetMaximum(dividend->GetMaximum()*1.1);   
    // float min=dividend->SetMinimum(dividend->GetMinimum()*1.3);

     dividend->SetStats(0);
     format_h(dividend,kBlack);
     gStyle->SetOptStat(0);
     dividend->Draw();
     dividend->SetTitle("");
     //           c1->cd();
        upperPad->cd();
  TLatex *tex2 = new TLatex(0.61,1.01,"#bf{552.67 pb^{-1} at 13 TeV}");
  tex2->SetNDC();
  tex2->SetTextAngle(0);
  tex2->SetTextFont(42);
  tex2->SetTextAlign(11);
  tex2->SetTextSize(0.06);
  tex2->Draw();
        TLegend *leg = new TLegend(0.71,0.98,1.1,0.467,NULL,"brNDC");//1)move to left and expand 2)move down and vertically seperate 3)move left and shrink 4) shrink vertically
        leg->SetFillStyle ( 0);
        leg->SetFillColor ( 0);
        leg->SetBorderSize( 0);
 //       leg->AddEntry(hs1,"M400","L");
        leg->SetLineWidth(1);
        leg->SetTextSize(0.056);
        leg->SetFillColor(0);
        leg->Draw();
//        leg->AddEntry(hs2,"M600","F");
 //       leg->AddEntry(hs3,"M800","L");
        leg->AddEntry(hdata,"Data","LP");
        leg->AddEntry(htt,"T#bar{T}","F");
        leg->AddEntry(hwj,"W+J","F");
        leg->AddEntry(hst_tch_at,"S-T","F");
        leg->AddEntry(hzz,"ZZ","F");
        leg->AddEntry(hww,"WW","F");
        leg->AddEntry(hqcd,"QCD","F");
	
}



void table()
{
 plot_shape("events_eachCut","events_eachCut","No. of Events");
        c1->SaveAs("results/events_eachCut.png");
        c1->SaveAs("results/events_eachCut.root");
}



