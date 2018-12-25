#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"

         void signifi_massRange() {
         const int nBins = 9;

         Double_t ratio[nBins] = {40.8656, 45.5842, 48.4062, 49.6622, 50.3697, 50.5891, 50.4788, 50.3164, 50.0271};
         TCanvas *c1 = new TCanvas("c1","Fitting Demo",10,10,700,500);
         c1->SetFillColor(33);
         c1->SetFrameFillColor(41);
         c1->SetGrid();
         float hm_var=430.4,hm_min=50.,hm_max=210.4;
         for (int i=hm_min;i<hm_max;i+=20){
         cout<<hm_var-i<<"-"<<hm_var+i<<endl;
         }
         TString naming[]=
         {"380.4-480.4","360.4-500.4","340.4-520.4","320.4-540.4","300.4-560.4","280.4-580.4","260.4-600.4","240.4-620.4","220.4-640.4"};
         TH2F *histo = new TH2F("histo","  ",9,0.,9.,14,40.,53.);
         histo->SetMarkerStyle(21);
         histo->SetMarkerSize(0.8);
         histo->SetStats(0);
         histo->GetXaxis()->SetTitle("Higgs Mass Range");
         
         for(int i=0; i < nBins;  i++){
         cout<<"i is: "<<i<<"this is bins cont: "<<ratio[i]<<endl;
         histo->Fill(i,ratio[i]);}

         for (int i=0;i<sizeof(naming)/sizeof(TString);i++){
         histo->GetXaxis()->SetBinLabel(i+1,naming[i]);}
         histo->Draw();   
         
            TLatex *tex1 = new TLatex(0.03,.5,"#frac{S}{#sqrt{B}}");
            tex1->SetNDC();
            tex1->SetTextAngle(0);
            tex1->SetTextFont(42);
            tex1->SetTextSize(0.035);
            tex1->SetTextAlign(11);
            tex1->Draw();   
                           
            TLatex *tex2 = new TLatex(0.4,.95,"Higgs Mass Range v/s #frac{S}{#sqrt{B}}");
            tex2->SetNDC();
            tex2->SetTextAngle(0);
            tex2->SetTextFont(42);
            tex2->SetTextSize(0.035);
            tex2->SetTextAlign(11);
            tex2->Draw();                                                                                                                                                                          // draw the legend
         TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
         legend->SetTextFont(72);
         legend->SetTextSize(0.04);
         legend->AddEntry(histo,"Data","lpe");
                                                                                                                                                                                             //legend->Draw();

}
