{
   auto c1 = new TCanvas("c1","mass of #{chi}1 vs #{sigma}_{LO}}",200,10,700,500);
   ifstream inf1;
   inf1.open("scan_run_02_10.txt");
   Float_t a,b,c;
//   c1->SetFillColor(42);
//   c1->SetGrid();
   c1->SetLogy();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   const Int_t n = 9;
   Double_t x[n]  = {750, 800, 850, 900, 950, 1000,1050,1100,1150};
   Double_t y[n]  = {1.925175e-01,1.285268e-01,8.747191e-02,6.071729e-02,
                     4.257310e-02,3.038538e-02,2.204414e-02,1.584721e-02,1.150209e-02};
   Double_t ex[n] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
   Double_t ey[n] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
   auto gr = new TGraphErrors(n,x,y,ex,ey);
   gr->GetYaxis()->SetRangeUser(-1000, 1.0);
   gr->SetTitle("mass of #chi1 vs #sigma_{LO}      M_{S1} = 1250 GeV");
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("ALP");
   
   gr->GetXaxis()->SetTitle("M_{#chi1}");
   gr->GetXaxis()->SetTitleOffset(0.9);
   gr->GetXaxis()->SetTitleSize(0.05);
   gr->GetYaxis()->SetTitle("#sigma_{LO}");
   gr->GetYaxis()->SetTitleOffset(1.1);

TLegend *leg = new TLegend(0.50,0.70,0.70,0.85,NULL,"brNDC");//1)move to left and expand 2)move down and vertically seperate 3)move left and shrink 4) shrink vertically
//        leg-> SetNColumns(3);
//        leg->SetFillStyle ( 0);
//        leg->SetFillColor ( 0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(3);
        leg->SetTextSize(0.05);
        leg->SetFillColor(0);
        leg->Draw();
        leg->AddEntry(gr,"#bar{#chi1} #chi1","L");

  }
