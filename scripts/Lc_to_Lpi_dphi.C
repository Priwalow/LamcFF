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
    
    double lend=0, rend=2*TMath::Pi(), MLambdac=2.28646; //lend=2.21, rend=2.36
    int Nbins=10;
    TCanvas *c1 = new TCanvas("c1","Lambda_c invariant mass",1600,900);
    TH1D* hdat = new TH1D("hdat","#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}",Nbins,lend,rend);
    TH1D* hsb = new TH1D("hsb","#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;
    
    
    double Ntot, Nsig, dNsig, Nbkg, dNbkg;
    
    Ntot = ch1dat -> Draw("philclam>>hdat","lcch==1 && abs(pvis)<0.05 && abs(ecms-sqrt(mvis*mvis+pvis*pvis))<0.05 && abs(mlc-2.28646)<0.015 && abs(ml-1.11568)<0.003 && abs(rmx-2.29)<0.1 && ((tag==3 && ((dstch==1 && abs(mdst-2.01026)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6)))) || (dstch==2 && abs(mdst-2.01026)<0.002 && abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3)))) || (tag==4 && dstch==1 && abs(mdst-2.00685)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))))","goff"); //"((tag!=11 && tag!=12) || abs(ml1-1.11568)<0.003) && pvis<0.05 && ecms-sqrt(mvis*mvis+pvis*pvis)<0.05 && abs(rmvis)<0.1
    ch1dat -> Draw("philclam>>hsb","lcch==1 && abs(pvis)<0.05 && abs(ecms-sqrt(mvis*mvis+pvis*pvis))<0.05 && abs(mlc-2.28646)>0.02 && abs(mlc-2.28646)<0.035 && abs(ml-1.11568)<0.003 && abs(rmx-2.29)<0.1 && ((tag==3 && ((dstch==1 && abs(mdst-2.01026)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6)))) || (dstch==2 && abs(mdst-2.01026)<0.002 && abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3)))) || (tag==4 && dstch==1 && abs(mdst-2.00685)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))))","goff"); //
    
    
  
    
  

   
    
   // hdat -> Sumw2();
   // hsb  -> Sumw2();
   // hdat -> Add(hsb,-1);
   //hsb  -> Scale(0.5);
    
    TF1* fdat = new TF1("fdat","[0]*(1-0.245*[1]*cos([2]+x))",lend,rend);
    fdat -> SetParameters(4,0.75,5.68);
    fdat -> SetParLimits(2,0,2*TMath::Pi());
    
    TFitResultPtr fitResult;
    fitResult = hdat -> Fit("fdat","L S M N","goff"); //L S M N Q
    double * par;
    par = fdat -> GetParameters();

     double chisq = fdat->GetChisquare(); 
    int ndf = Nbins - fdat ->GetNumberFreeParameters(); // nbins - npar
    cout << "ChiSquare / ndf: " << chisq << " / " << ndf << " = " << chisq/ndf << endl;
    
    hdat -> GetXaxis()-> SetTitle("#Delta in #Lambda_{c} system");
    hdat -> GetXaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetXaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitle(Form("Events /  %.1f #pi ",binw/TMath::Pi()));
    hdat -> GetYaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetYaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitleOffset(0.7);
    //hdat -> GetYaxis()->CenterTitle(true);
    hdat -> GetXaxis()->SetTickSize(0.04);
    hdat -> SetMarkerStyle(20);
    hdat -> SetMarkerSize(1.2);
    hdat -> SetMarkerColor(1);
    hdat -> SetLineColor(1);
    hdat -> SetLineWidth(2);
    hdat -> SetMinimum(0);
    hdat -> Draw("e p");
    
  
    hsb -> SetLineColor(8);
    hsb -> SetLineWidth(4);
    hsb -> Draw("hist same");
     /* 
     fbkg -> SetLineStyle(2);
    fbkg -> SetLineColor(12);
    fbkg -> SetLineWidth(4);
    fbkg -> DrawCopy("same");
    
    
    fsig -> SetLineColor(4);
    fsig -> SetLineWidth(4);
    //fsig -> Draw("same");
    
     */
    fdat -> SetLineColor(2);
    fdat -> SetLineWidth(4);
    fdat-> DrawCopy("same");
   
    
    TLegend* leg = new TLegend(0.77,0.75,0.89,0.95);
    leg -> AddEntry("hdat","data","E P"); //with sideband subtraction
    leg -> AddEntry("hsb","sideband in M(#Lambda#pi)");
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
