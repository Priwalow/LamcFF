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
    
    double lend=-0.5, rend=4, MLambdac=2.28646; //lend=2.21, rend=2.36
    int Nbins=90;
    TCanvas *c1 = new TCanvas("c1","Lambda_c invariant mass",1600,900);
    TH1D* hdat = new TH1D("hdat","#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}",Nbins,lend,rend);
    TH1D* hsb = new TH1D("hsb","sideband",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;

    
    double Ntot=0, Nsig, dNsig, Nbkg3s, dNbkg3s;
    TCut Mwindow = Form("ecms-sqrt(mvis*mvis+pvis*pvis) > %lf && ecms-sqrt(mvis*mvis+pvis*pvis) < %lf",lend,rend);
    TCut commonLcLpiCut = "lcch==1 && abs(ml-1.11568)<0.003 && abs(rmx-2.29)<0.1 && ((tag==3 && ((dstch==1 && abs(mdst-2.01026)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6)))) || (dstch==2 && abs(mdst-2.01026)<0.002 && abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3)))) || (tag==4 && dstch==1 && abs(mdst-2.00685)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))))";
    TCut rmxlsbCut = "mlc-2.28646 < -0.02 && mlc-2.28646 > -0.035";
    TCut rmxrsbCut = "mlc-2.28646 > 0.02 && mlc-2.28646 < 0.035";
    Ntot +=  ch1dat -> Draw("ecms-sqrt(mvis*mvis+pvis*pvis)>>+hdat","abs(mlc-2.28646)<0.015"+ commonLcLpiCut+Mwindow,"goff");
    ch1dat -> Draw("ecms-sqrt(mvis*mvis+pvis*pvis)>>+hsb",rmxrsbCut+commonLcLpiCut+Mwindow,"goff");
    ch1dat -> Draw("ecms-sqrt(mvis*mvis+pvis*pvis)>>+hsb",rmxlsbCut+commonLcLpiCut+Mwindow,"goff");
    
    
    
    cout << "Total number of events: " << Ntot << endl;
    
  
    
  
    hdat -> GetXaxis()-> SetTitle("E_{miss}(#Lambda#pi^{+}#bar{p}#bar{D}^{*0})||E_{miss}(#Lambda#pi^{+}#bar{p}D^{*-}#pi^{+})[GeV]");
    hdat -> GetXaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetXaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetXaxis()-> SetTitleOffset(1.1);
    hdat -> GetYaxis()-> SetTitle(Form("Events /  %.2f GeV ",binw));
    hdat -> GetYaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetYaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitleOffset(0.7);
    //hdat -> GetYaxis()->CenterTitle(true);
    hdat -> GetXaxis()->SetTickSize(0.04);
    hdat -> SetMarkerStyle(20);
    hdat -> SetMarkerSize(1.5);
    hdat -> SetMarkerColor(1);
    hdat -> SetLineColor(1);
    hdat -> SetLineWidth(2);
    hdat -> SetMinimum(0);
    hdat -> Draw("e p");
    
    hsb -> SetLineColor(8);
    hsb -> SetLineWidth(4);
    hsb -> Draw("same");


    TLegend* leg = new TLegend(0.77,0.75,0.89,0.95);
    leg -> AddEntry("hdat","data","e p");
    leg -> AddEntry("hsb","M(#Lambda#pi^{+}) sideband","l");
    leg -> SetBorderSize(0);
    leg -> SetTextSize(axisFontSize);
    leg -> Draw("same");
    

   
    
} 
 
