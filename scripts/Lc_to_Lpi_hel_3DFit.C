
TCut totCUT = "lcch==1 && abs(mlc-2.28646)<0.02 && abs(pvis)<0.05 && abs(ecms-sqrt(mvis*mvis+pvis*pvis))<0.05 && abs(ml-1.11568)<0.003 && abs(rmx-2.29)<0.1 && ((tag==3 && ((dstch==1 && abs(mdst-2.01026)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6)))) || (dstch==2 && abs(mdst-2.01026)<0.002 && abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3)))) || (tag==4 && dstch==1 && abs(mdst-2.00685)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))))";

double lend=-1, rend=1, MLambdac=2.28646, alphaLam=0.732, alphaLam_c = -0.84; 
int Nbins=3;

TString datapath = "../analysis/hmerge/";    
TString mcpath = "../mc_analysis/hmerge/";


TH3D* hmceff = new TH3D("hmceff","#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}",Nbins,-TMath::Pi(),TMath::Pi(),Nbins,lend,rend,Nbins,lend,rend);
TH3D* hmcgen = new TH3D("hmcgen","#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}",Nbins,-TMath::Pi(),TMath::Pi(),Nbins,lend,rend,Nbins,lend,rend);
TH3D* hmcsel = new TH3D("hmcsel","#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}",Nbins,-TMath::Pi(),TMath::Pi(),Nbins,lend,rend,Nbins,lend,rend);

double Ntot, Nmcgen, Nmcsel;

void histeffgen() 
{
    TChain* ch1mc = new TChain("h1");
    ch1mc -> Add(mcpath+"*.root");
    
    TChain* ch2mc = new TChain("h2");
    ch2mc -> Add(mcpath+"*.root");
    
    Nmcsel = ch1mc -> Draw("hlc:hl:philclam>>hmcsel",totCUT,"goff");
    Nmcgen = ch2mc -> Draw("hlc:hl:philclam>>hmcgen","lcch==1","goff");
    
    hmcsel->Sumw2();
    hmcgen->Sumw2();
    *hmceff = *hmcsel/(*hmcgen);
}

Double_t myfdat(Double_t* x, Double_t *par)
{
   double xx = x[0], yy = x[1], zz=x[2];
   int binx = hmceff->GetXaxis()->FindBin(xx), biny = hmceff->GetYaxis()->FindBin(yy), binz = hmceff->GetZaxis()->FindBin(zz);
   double eff = hmceff->GetBinContent(binx,biny,binz);
   return 1*par[0]*(1+alphaLam*alphaLam_c*yy+par[1]*(zz*(alphaLam_c+alphaLam*yy)-sqrt((1-yy*yy)*(1-zz*zz))*alphaLam*sqrt(1-alphaLam_c*alphaLam_c)*cos(par[2]+xx)));
}

Double_t myfdat1(Double_t* x, Double_t *par)
{
   double xx = x[0], yy = x[1];
   return 0;
}

Double_t myfdat2(Double_t* x, Double_t *par)
{
   double xx = x[0], yy = x[1];
   return 0;
}

Double_t myfdat3(Double_t* x, Double_t *par)
{
   double xx = x[0], yy = x[1];
   return 0;
}


void Lc2Lpi3Dfit()
{
    histeffgen();
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.07);
    gStyle->SetErrorX(0);
    double axisFontSize = 0.05;
    

    TChain* ch1dat = new TChain("h1");
    ch1dat -> Add(datapath+"*.root");
    

    TH3D* hmain = new TH3D("hmain","#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}",Nbins,-TMath::Pi(),TMath::Pi(),Nbins,lend,rend,Nbins,lend,rend);

    double hwidth = rend-lend, binw = hwidth/Nbins;
    
    
    

    Ntot = ch1dat -> Draw("hlc:hl:philclam>>hmain",totCUT,"goff"); //"((tag!=11 && tag!=12) || abs(ml1-1.11568)<0.003) && pvis<0.05 && ecms-sqrt(mvis*mvis+pvis*pvis)<0.05 && abs(rmvis)<0.1
    hmain -> Sumw2();
    

    
    cout << hmcsel->GetBinContent(1,1,1)<< " +- " << hmcsel->GetBinError(1,1,1) << " / " << hmcgen->GetBinContent(1,1,1)<< " +- " << hmcgen->GetBinError(1,1,1) << " = " << hmceff->GetBinContent(1,1,1)<< " +- " << hmceff->GetBinError(1,1,1) << endl;
    
    cout << hmain->GetBinContent(1,1,1)<< " +- " << hmain->GetBinError(1,1,1) << " / " << hmceff->GetBinContent(1,1,1)<< " +- " << hmceff->GetBinError(1,1,1) << " = ";
    //*hmain = *hmain/(*hmceff);
    cout << hmain->GetBinContent(1,1,1)<< " +- " << hmain->GetBinError(1,1,1) << endl;
    
    /*TEfficiency * pEff = new TEfficiency(*hmcsel,*hmcgen);
    int effbin;
    for(int i =0; i < Nbins; i++)
    {
        for(int j =0; j < Nbins; j++)
        {
            for(int k=0; k< Nbins; k++)
            {
                effbin = pEff -> GetGlobalBin(i+1,j+1,k+1);
                eff[i*Nbins*Nbins+j*Nbins+k] = pEff -> GetEfficiency(effbin);
                effUerr[i*Nbins*Nbins+j*Nbins+k] = pEff -> GetEfficiencyErrorUp(effbin);
                effLerr[i*Nbins*Nbins+j*Nbins+k] = pEff -> GetEfficiencyErrorLow(effbin);
                effconst+=pEff -> GetEfficiency(effbin);
                effBins[i*Nbins*Nbins+j*Nbins+k] = i*Nbins*Nbins+j*Nbins+k+1;
            }
        }
    }
    effconst=effconst/(Nbins*Nbins*Nbins);
    for(int i=0;i<Nbins*Nbins*Nbins;i++)
    {
        //eff[i]=eff[i]/effconst;
        cout << eff[i] << " " << effUerr[i] << " " <<effLerr[i] << endl;
    }
    
    TCanvas *c4 = new TCanvas("c4","efficiency",1600,900);
    TGraphAsymmErrors * gr = new TGraphAsymmErrors(Nbins*Nbins*Nbins,effBins,eff,0,0,effLerr,effUerr);
    gr->SetName("gr");
    gr->SetTitle("Efficiency");
    gr->GetXaxis()->SetTitle("bin");
    gr->GetYaxis()->SetTitle("#epsilon");
    gr->SetMarkerStyle(8);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerSize(1.2);
    gr->Draw("ap");
    
    TH3D* hmceff = new TH3D("hmceff","#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}",Nbins,0,2*TMath::Pi(),Nbins,lend,rend,Nbins,lend,rend);
    for(int i =0; i < Nbins; i++)
    {
        for(int j =0; j < Nbins; j++)
        {
            for(int k=0; k< Nbins; k++)
            {
                effbin = hmceff -> GetBin(i+1,j+1,k+1)
            }
        }
    }
    */
    
    
    TH1* hlcl = hmain -> Project3D("zy");//new TH2D("hlcl","#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}",Nbins,lend,rend,Nbins,lend,rend);
    TH1* hlcphi = hmain -> Project3D("zx"); //new TH2D("hlcphi","#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}",Nbins,0,2*TMath::Pi(),Nbins,lend,rend);
    TH1* hlphi = hmain -> Project3D("yx");//new TH2D("hlphi","#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}",Nbins,0,2*TMath::Pi(),Nbins,lend,rend);
  
    //ch1dat -> Draw("hlc:hl>>hlcl",totCUT,"goff");
    //ch1dat -> Draw("hlc:philclam>>hlcphi",totCUT,"goff");
    //ch1dat -> Draw("hl:philclam>>hlphi",totCUT,"goff");

    
    TF3* fmain = new TF3("fmain",myfdat,-TMath::Pi(),TMath::Pi(),lend,rend,lend,rend,3);
    fmain -> SetParameters(4,0.5,0);
    fmain -> SetParLimits(2,-TMath::Pi()-0.005,TMath::Pi()+0.005);
    
    
    TFitResultPtr fitResult;
    fitResult = hmain -> Fit("fmain","L S M N","goff"); //L S M N Q
    double * par;
    par = fmain -> GetParameters();
    double chisq = fmain->GetChisquare(), ndf = fmain->GetNDF(); 
    cout << "ChiSquare / ndf: " << chisq << " / " << ndf << " = " << chisq/ndf << "; Prob = " << fmain -> GetProb()<<  endl;
    
    
    TF2* flcl = new TF2("flcl",Form("[0]*(1+%lf*x+[1]*y*(%lf+%lf*x))",alphaLam*alphaLam_c,alphaLam_c,alphaLam),lend,rend,lend,rend);
    TF2* flcphi = new TF2("flcphi",Form("[0]*(1+[1]*(%lf*y-%lf*cos([2]+x)*sqrt(1-y*y)))",
                                        alphaLam_c,alphaLam*sqrt(1-alphaLam_c*alphaLam_c)*TMath::Pi()/4),-TMath::Pi(),TMath::Pi(),lend,rend);
    TF2* flphi = new TF2("flphi",Form("[0]*(1+%lf*y-[1]*%lf*cos([2]+x)*sqrt(1-y*y))",
                                      alphaLam_c*alphaLam,alphaLam*sqrt(1-alphaLam_c*alphaLam_c)*TMath::Pi()/4),-TMath::Pi(),TMath::Pi(),lend,rend);
    
    flcl -> SetParameters(2*TMath::Pi()*par[0],par[1]);
    flcphi -> SetParameters(2*par[0],par[1],par[2]);
    flphi -> SetParameters(2*par[0],par[1],par[2]);
    
    TCanvas *c1 = new TCanvas("c1","3Dfit_1",1600,900);
    hlcl -> GetXaxis()-> SetTitle("-cos(p_{p},p_{#Lambda_{c}}) in #Lambda system");
    hlcl -> GetXaxis()-> SetTitleSize(axisFontSize);
    hlcl -> GetXaxis()-> SetLabelSize(axisFontSize);
    hlcl -> GetXaxis()-> SetTitleOffset(2);
    hlcl -> GetYaxis()-> SetTitle("-cos(p_{#Lambda},p_{e^{+}e^{-}}) in #Lambda_{c} system");
    hlcl -> GetYaxis()-> SetTitleSize(axisFontSize);
    hlcl -> GetYaxis()-> SetLabelSize(axisFontSize);
    hlcl -> GetYaxis()-> SetTitleOffset(1.5);
    hlcl -> GetYaxis()->CenterTitle(true);
    hlcl -> GetZaxis()-> SetTitle(Form("Events /  %.1f #times %.1f",binw,binw));
    hlcl -> GetZaxis()-> SetTitleSize(axisFontSize);
    hlcl -> GetZaxis()-> SetLabelSize(axisFontSize);
    hlcl -> GetZaxis()-> SetTitleOffset(0.8);
    hlcl -> GetXaxis()->SetTickSize(0.04);
    hlcl -> SetMarkerStyle(20);
    hlcl -> SetMarkerSize(1.2);
    hlcl -> SetMarkerColor(1);
    hlcl -> SetLineColor(1);
    hlcl -> SetLineWidth(2);
    hlcl -> SetMinimum(0);
    hlcl -> Draw("lego");    
    
    flcl -> SetLineColor(2);
    flcl -> SetLineWidth(4);
    flcl-> DrawCopy("same");
    
    TCanvas *c2 = new TCanvas("c2","3Dfit_2",1600,900);
    hlcphi -> GetXaxis()-> SetTitle("#Delta in #Lambda_{c} system");
    hlcphi -> GetXaxis()-> SetTitleSize(axisFontSize);
    hlcphi -> GetXaxis()-> SetLabelSize(axisFontSize);
    hlcphi -> GetXaxis()-> SetTitleOffset(1.3);
    hlcphi -> GetXaxis()->CenterTitle(true);
    hlcphi -> GetYaxis()-> SetTitle("-cos(p_{#Lambda},p_{e^{+}e^{-}}) in #Lambda_{c} system");
    hlcphi -> GetYaxis()-> SetTitleSize(axisFontSize);
    hlcphi -> GetYaxis()-> SetLabelSize(axisFontSize);
    hlcphi -> GetYaxis()-> SetTitleOffset(1.5);
    hlcphi -> GetYaxis()->CenterTitle(true);
    hlcphi -> GetZaxis()-> SetTitle(Form("Events /  %.1f#pi #times %.1f",binw,binw));
    hlcphi -> GetZaxis()-> SetTitleSize(axisFontSize);
    hlcphi -> GetZaxis()-> SetLabelSize(axisFontSize);
    hlcphi -> GetZaxis()-> SetTitleOffset(0.8);
    hlcphi -> GetXaxis()->SetTickSize(0.04);
    hlcphi -> SetMarkerStyle(20);
    hlcphi -> SetMarkerSize(1.2);
    hlcphi -> SetMarkerColor(1);
    hlcphi -> SetLineColor(1);
    hlcphi -> SetLineWidth(2);
    hlcphi -> SetMinimum(0);
    hlcphi -> Draw("lego");
    
    flcphi -> SetLineColor(2);
    flcphi -> SetLineWidth(4);
    flcphi-> DrawCopy("same");
    
    TCanvas *c3 = new TCanvas("c3","3Dfit_3",1600,900);
    hlphi -> GetXaxis()-> SetTitle("#Delta in #Lambda_{c} system");
    hlphi -> GetXaxis()-> SetTitleSize(axisFontSize);
    hlphi -> GetXaxis()-> SetLabelSize(axisFontSize);
    hlphi -> GetXaxis()-> SetTitleOffset(1.3);
    hlphi -> GetXaxis()->CenterTitle(true);
    hlphi -> GetYaxis()-> SetTitle("-cos(p_{p},p_{#Lambda_{c}}) in #Lambda system");
    hlphi -> GetYaxis()-> SetTitleSize(axisFontSize);
    hlphi -> GetYaxis()-> SetLabelSize(axisFontSize);
    hlphi -> GetYaxis()-> SetTitleOffset(1.5);
    hlphi -> GetYaxis()->CenterTitle(true);
    hlphi -> GetZaxis()-> SetTitle(Form("Events /  %.1f#pi #times %.1f",binw,binw));
    hlphi -> GetZaxis()-> SetTitleSize(axisFontSize);
    hlphi -> GetZaxis()-> SetLabelSize(axisFontSize);
    hlphi -> GetZaxis()-> SetTitleOffset(0.8);
    
    hlphi -> GetXaxis()->SetTickSize(0.04);
    hlphi -> SetMarkerStyle(20);
    hlphi -> SetMarkerSize(1.2);
    hlphi -> SetMarkerColor(1);
    hlphi -> SetLineColor(1);
    hlphi -> SetLineWidth(2);
    hlphi -> SetMinimum(0);
    hlphi -> Draw("lego");
    
    flphi -> SetLineColor(2);
    flphi -> SetLineWidth(4);
    flphi-> DrawCopy("same");
} 
