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
    
    double lend=0, rend=5, MLambdac=2.28646; //lend=2.21, rend=2.36
    int Nbins=25;
    TCanvas *c1 = new TCanvas("c1","Lambda_c invariant mass",1600,900);
    TH1D* hdat = new TH1D("hdat","#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}",Nbins,lend,rend);
    TH1D* hsb = new TH1D("hsb","#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;
    
    
    double Ntot, Nsig, dNsig, Nbkg, dNbkg;
    Ntot = ch1dat -> Draw("lcp2dcm>>hdat","lcch==2 && abs(pvis)<0.05 && abs(ecms-sqrt(mvis*mvis+pvis*pvis))<0.05 && abs(mlc-2.28646)<0.035 && abs(ml-1.11568)<0.003 && abs(rmx-2.29)<0.1 && ((tag==3 && ((dstch==1 && abs(mdst-2.01026)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6)))) || (dstch==2 && abs(mdst-2.01026)<0.002 && abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3)))) || (tag==4 && dstch==1 && abs(mdst-2.00685)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))))","goff"); //
    ch1dat -> Draw("lcp2dcm>>hsb","lcch==2 && abs(pvis)<0.05 && abs(ecms-sqrt(mvis*mvis+pvis*pvis))<0.05 && abs(mlc-2.28646)>0.03 && abs(mlc-2.28646)<0.055 && abs(ml-1.11568)<0.003 && abs(rmx-2.29)<0.1 && ((tag==3 && ((dstch==1 && abs(mdst-2.01026)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6)))) || (dstch==2 && abs(mdst-2.01026)<0.002 && abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3)))) || (tag==4 && dstch==1 && abs(mdst-2.00685)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))))","goff"); //

   // hdat -> Sumw2();
   // hsb  -> Sumw2();
   // hdat -> Add(hsb,-1);
   //hsb  -> Scale(0.5);
    
    hdat -> GetXaxis()-> SetTitle("p_{#pi^{+}} in e^{+}e^{-} system [GeV]");
    hdat -> GetXaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetXaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitle(Form("Events /  %.1f GeV",binw));
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
    //hsb -> Draw("hist same");

    
    TLegend* leg = new TLegend(0.77,0.75,0.89,0.95);
    leg -> AddEntry("hdat","data","E P"); //with sideband subtraction
    //leg -> AddEntry("hsb","sideband in M(#Lambda#pi)");
    leg -> SetBorderSize(0);
    leg -> SetTextSize(axisFontSize);
    leg -> Draw("same");

    
} 
