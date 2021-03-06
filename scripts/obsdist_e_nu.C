{
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.07);
    gStyle->SetErrorX(0);
    
    double axisFontSize = 0.065;
    
    //TString datapath = "../analysis/hmerge/";
    TString datapath = "../mc_analysis/hmerge/";    
    TChain* ch1dat = new TChain("h1"); //without pi0
    ch1dat -> Add(datapath+"*.root");
    
    
    double qsqmax = 1.370726;
    TCut gencut = "rmvis*rmvis<0.1 && abs(rmx-2.29)<0.1 && q*fabs(q)>0 &&  q*fabs(q)<1.37 && lcch==3 && abs(ml-1.11568)<0.003 && pvis>0.05 && ecms-sqrt(mvis*mvis+pvis*pvis)>0.05 && ((tag==3 && ((dstch==1 && abs(mdst-2.01026)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6)))) || (dstch==2 && abs(mdst-2.01026)<0.002 && abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3)))) || (tag==4 && dstch==1 && abs(mdst-2.00685)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))))";
    
    TCanvas *c0 = new TCanvas("c1","Lc -> Lenu observables",740,600);
    c0 -> Divide(2,2);
    
    TH1D* hqsq = new TH1D("hqsq","q^{2}/q^{2}_{max}",10,0,1);
    TH1D* hcosw = new TH1D("hcosw","cos#theta_{W}",10,-1,1);
    TH1D* hcosl = new TH1D("hcosl","cos#theta_{#Lambda}",10,-1,1);
    TH1D* hcoschi = new TH1D("hcoschi","#chi",10,-180,180);
    
    ch1dat -> Draw("q*fabs(q)/1.370726>>hqsq",gencut,"goff");
    ch1dat -> Draw("hw>>hcosw",gencut,"goff");
    ch1dat -> Draw("hl>>hcosl",gencut,"goff");
    ch1dat -> Draw(Form("chi*180/%lf>>hcoschi",TMath::Pi()),gencut,"goff");
    
    hqsq -> GetXaxis()->SetTitle("q^{2}/q^{2}_{max}");
    hqsq -> GetXaxis()-> SetTitleSize(axisFontSize);
    hqsq -> GetXaxis()-> SetLabelSize(axisFontSize);
    hqsq -> GetYaxis()-> SetTitle(Form("Events / 0.1"));
    hqsq -> GetYaxis()-> SetTitleSize(axisFontSize);
    hqsq -> GetYaxis()-> SetLabelSize(axisFontSize);
    hqsq -> GetYaxis()-> SetTitleOffset(0.8);
    //hdat -> GetYaxis()->CenterTitle(true);
    hqsq -> GetXaxis()->SetTickSize(0.04);
    hqsq -> SetMarkerStyle(20);
    hqsq -> SetMarkerSize(1.5);
    hqsq -> SetMarkerColor(1);
    hqsq -> SetLineColor(1);
    hqsq -> SetLineWidth(2);
    hqsq -> SetMinimum(0);
    c0 -> cd(1);
    hqsq -> Draw("e p");
    
    hcosw -> GetXaxis()->SetTitle("cos#theta_{W}");
    hcosw -> GetXaxis()-> SetTitleSize(axisFontSize);
    hcosw -> GetXaxis()-> SetLabelSize(axisFontSize);
    hcosw -> GetYaxis()-> SetTitle(Form("Events / 0.2"));
    hcosw -> GetYaxis()-> SetTitleSize(axisFontSize);
    hcosw -> GetYaxis()-> SetLabelSize(axisFontSize);
    hcosw -> GetYaxis()-> SetTitleOffset(0.8);
    //hdat -> GetYaxis()->CenterTitle(true);
    hcosw -> GetXaxis()->SetTickSize(0.04);
    hcosw -> SetMarkerStyle(20);
    hcosw -> SetMarkerSize(1.5);
    hcosw -> SetMarkerColor(1);
    hcosw -> SetLineColor(1);
    hcosw -> SetLineWidth(2);
    hcosw -> SetMinimum(0);
    c0 -> cd(2);
    hcosw -> Draw("e p");
    
    hcosl -> GetXaxis()->SetTitle("cos#theta_{#Lambda}");
    hcosl -> GetXaxis()-> SetTitleSize(axisFontSize);
    hcosl -> GetXaxis()-> SetLabelSize(axisFontSize);
    hcosl -> GetYaxis()-> SetTitle(Form("Events / 0.2"));
    hcosl -> GetYaxis()-> SetTitleSize(axisFontSize);
    hcosl -> GetYaxis()-> SetLabelSize(axisFontSize);
    hcosl -> GetYaxis()-> SetTitleOffset(0.8);
    //hdat -> GetYaxis()->CenterTitle(true);
    hcosl -> GetXaxis()->SetTickSize(0.04);
    hcosl -> SetMarkerStyle(20);
    hcosl -> SetMarkerSize(1.5);
    hcosl -> SetMarkerColor(1);
    hcosl -> SetLineColor(1);
    hcosl -> SetLineWidth(2);
    hcosl -> SetMinimum(0);
    c0 -> cd(3);
    hcosl -> Draw("e p");
    
    hcoschi -> GetXaxis()->SetTitle("#chi[#circ]");
    hcoschi -> GetXaxis()-> SetTitleSize(axisFontSize);
    hcoschi -> GetXaxis()-> SetLabelSize(axisFontSize);
    hcoschi -> GetYaxis()-> SetTitle(Form("Events / 0.2"));
    hcoschi -> GetYaxis()-> SetTitleSize(axisFontSize);
    hcoschi -> GetYaxis()-> SetLabelSize(axisFontSize);
    hcoschi -> GetYaxis()-> SetTitleOffset(0.8);
    //hdat -> GetYaxis()->CenterTitle(true);
    hcoschi -> GetXaxis()->SetTickSize(0.04);
    hcoschi -> SetMarkerStyle(20);
    hcoschi -> SetMarkerSize(1.5);
    hcoschi -> SetMarkerColor(1);
    hcoschi -> SetLineColor(1);
    hcoschi -> SetLineWidth(2);
    hcoschi -> SetMinimum(0);
    c0 -> cd(4);
    hcoschi -> Draw("e p");
    
    
    
}
 
