    {
    //gROOT->Reset();
    TCanvas *c1 = new TCanvas("c1","PDFs",200,10,700,500);
    gStyle->SetOptStat(0);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(0);
    c1->SetFrameBorderMode(0);
    c1->SetFrameLineWidth(2);
    c1->SetTickx(1);
    gStyle->SetEndErrorSize(3);
    
    TH2F *hr =     new TH2F("","",9,0,9,14,40,54);
    hr->GetYaxis()->SetLabelOffset(1);
    hr->GetYaxis()->SetNdivisions(0);

    hr->SetXTitle("tqz cross section #sigma (pb)");
    hr->Draw();

    gPad->RedrawAxis();
       
    Int_t n = 5;
    //    Double_t x_sec[]  = {.02870,.01976, 0.00475, 0.04911, 0.02478,}; //vertical position of values
    //    Double_t x_sec[]  = {0.01456,.00383, -0.017565, 0.02274079, 0.0067816,}; 
    Int_t nBins=9;
    Double_t ratio[nBins]  = {40.8656, 45.5842, 48.4062, 49.6622, 50.3697, 50.5891, 50.4788, 50.3164 , 50.0271};
    Double_t ratio1[nBins]  = {1,2,3,4,5,6,7,8,9} ; 

    std::string names[nBins]=  {"380.4-480.4","360.4-500.4","340.4-520.4","320.4-540.4","300.4-560.4","280.4-580.4","260.4-600.4","240.4-620.4","220.4-640.4"};

    // TString naming[]=    {"380.4-480.4","360.4-500.4","340.4-520.4","320.4-540.4","300.4-560.4","280.4-580.4","260.4-600.4","240.4-620.4","220.4-640.4"};

 //      for (int i=0;i<sizeof(naming)/sizeof(TString);i++){
 //	histo->GetXaxis()->SetBinLabel(i+1,naming[i]);}

    //Double_t x_sec[]  = {0.024, 0.012,  0.035,-0.006, 0.006}; 

      //      for(int i=0; i < nBins;  i++){
      //	cout<<"i is: "<<i<<"this is bins cont: "<<ratio[i]<<endl;
      //	histo->Fill(i,ratio[i]);}
    
    gr3 = new TGraphErrors(nBins,ratio1,ratio,0,0);
    TAxis *ax = gr3->GetHistogram()->GetXaxis();
    Double_t x1 = ax->GetBinLowEdge(1);
    Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
    gr3->GetHistogram()->GetXaxis()->Set(9,x1,x2);



    gr3->SetMarkerColor(kRed);
    gr3->SetMarkerSize(0.9);
    gr3->SetMarkerStyle(20);
    gr3->SetLineWidth(2);
    gr3->SetLineColor(kRed);

    for(Int_t k=0;k<nBins;k++){
      gr3->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    }

    gr3->Draw("");
    
    c1->SaveAs("gulaaa.png");        
}
                                                                            
