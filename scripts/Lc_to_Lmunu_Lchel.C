{
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.07);
    gStyle->SetErrorX(0);
    double axisFontSize = 0.065;
    
    
    
    
    TString datapath = "../analysis/hmerge/";    
    TChain* ch1dat = new TChain("h1"); //without pi0
    ch1dat -> Add(datapath+"*.root");
    
    double lend=-1, rend=1, MLambdac=2.28646; //lend=2.21, rend=2.36
    int Nbins=10;
    TCanvas *c1 = new TCanvas("c1","Lambda_c invariant mass",1600,900);
    TH1D* hdat = new TH1D("hdat","#Lambda_{c} #rightarrow #Lambdae#nu",Nbins,lend,rend);
    TH1D* hsb = new TH1D("hsb","#Lambda_{c} #rightarrow #Lambda#mu#nui",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;
    
    
    double Ntot, Nsig, dNsig, Nbkg, dNbkg;
    Ntot = ch1dat -> Draw("hlc>>hdat","lcch == 3 && abs(ml-1.11568)<0.003 && ((tag!=11 && tag!=12) || abs(ml1-1.11568)<0.003) && mlc < 2.13 && plept>0.7 && abs(lepcost)<0.7 && fox>0.2 && plc>1.4 && abs(rmx-2.2901)<0.0438*2 && abs(rmvis*rmvis)<0.05","goff"); //
    ch1dat -> Draw("hlc>>hsb","lcch == 3 && abs(ml-1.11568)<0.003 && ((tag!=11 && tag!=12) || abs(ml1-1.11568)<0.003) && mlc < 2.13 && plept>0.7 && abs(lepcost)<0.7 && fox>0.2 && plc>1.4 && abs(rmx-2.2901)<0.0438*2 && abs(rmvis*rmvis)>0.15 && abs(rmvis*rmvis)<0.20","goff"); //abs(rmx-2.2969)>0.0468*3 && abs(rmx-2.2969)<0.0468*5
    
    
  
    
  

   
    
    hdat -> Sumw2();
    hsb  -> Sumw2();
    hdat -> Add(hsb,-1);
 
    
    TF1* fdat = new TF1("fdat","[0]*(1-x*x)+[1]",lend,rend);
    TFitResultPtr fitResult;
    fitResult = hdat -> Fit("fdat","L S M N","goff"); //L S M N Q
    double * par;
    par = fdat -> GetParameters();
    
     double chisq = fdat->GetChisquare(); 
    int ndf = Nbins - fdat ->GetNumberFreeParameters(); // nbins - npar
    cout << "ChiSquare / ndf: " << chisq << " / " << ndf << " = " << chisq/ndf << endl;
    
    hdat -> GetXaxis()-> SetTitle("-cos(p_{#Lambda},p_{e^{+}e^{-}}) in #Lambda_{c} system");
    hdat -> GetXaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetXaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitle(Form("Events /  %.1f ",binw));
    hdat -> GetYaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetYaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitleOffset(0.8);
    //hdat -> GetYaxis()->CenterTitle(true);
    hdat -> GetXaxis()->SetTickSize(0.04);
    hdat -> SetMarkerStyle(20);
    hdat -> SetMarkerSize(1.2);
    hdat -> SetMarkerColor(1);
    hdat -> SetLineColor(1);
    hdat -> SetLineWidth(2);
    hdat -> SetMinimum(0);
    hdat -> Draw("e p");
    
  
   /*   hsb -> SetLineColor(8);
    hsb -> SetLineWidth(4);
    hsb -> Draw("same");
     
   fbkg -> SetLineStyle(2);
    fbkg -> SetLineColor(12);
    fbkg -> SetLineWidth(4);
    fbkg -> DrawCopy("same");
    
    
    fsig -> SetLineColor(4);
    fsig -> SetLineWidth(4);
    //fsig -> Draw("same");
    
   
    fdat -> SetLineColor(2);
    fdat -> SetLineWidth(4);
    fdat-> DrawCopy("same");
   */  
    
    TLegend* leg = new TLegend(0.77,0.75,0.89,0.95);
    leg -> AddEntry("hdat","data with sideband subtraction","E P");
    //leg -> AddEntry("hsb","sideband in M(#Lambda#pi)");
    // leg->AddEntry("fsig","Signal","l");
    //leg -> AddEntry("fdat","fit","l");
    //leg -> AddEntry("fbkg","background","l");
    leg -> SetBorderSize(0);
    leg -> SetTextSize(axisFontSize);
    leg -> Draw("same");
  /*  
    TLatex *tstatfit = new TLatex();
    tstatfit -> SetNDC();
    tstatfit -> SetTextColor(1);
    tstatfit -> SetTextSize(axisFontSize);
    tstatfit -> SetTextAngle(0);
    tstatfit -> DrawLatex(0.67, 0.65, Form("N_{sig} = %0.lf #pm %0.lf",Nsig, dNsig)); //
   // tstatfit -> DrawLatex(0.67, 0.59, Form("N_{bkg} = %0.lf #pm %0.lf",Nbkg, dNbkg)); //
    tstatfit -> DrawLatex(0.67, 0.59, Form("Mean_{sig} = %0.4lf #pm %0.4lf",par[1], fdat -> GetParError(1)));
    tstatfit -> DrawLatex(0.67, 0.53, Form("#sigma_{sig} = %0.4lf #pm %0.4lf",par[2], fdat -> GetParError(2)));

   */
    
} 
