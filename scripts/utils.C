#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <TMath.h>
#include <TRandom.h>



void BrFracMeas(double B1, double dB1, double B2, double dB2)
{
    double A, dA;
    A = B1/B2;
    dA = A*sqrt(dB1*dB1/(B1*B1)+dB2*dB2/(B2*B2));
    TRandom r; 
    r.SetSeed(0);
    double buff;
    TH1D* H = new TH1D("H", "BrFrac", 100,A-6*dA,A+6*dA);
    
    for(int i =0; i<100000; i++)
    {
        buff = r.Gaus(B1,dB1)/r.Gaus(B2,dB2);
        H -> Fill(buff);
    }
    H->DrawCopy();
    delete H;
}



void LctoLpi() 
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
    
    double lend=2.21, rend=2.36, MLambdac=2.28646; //lend=2.21, rend=2.36
    int Nbins=50;
    TCanvas *c1 = new TCanvas("c1","Lambda_c invariant mass",1600,900);
    TH1D* hdat = new TH1D("hdat","#Lambda_{c} #rightarrow #Lambda#pi",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;
    
    
    double Ntot, Nsig, dNsig, Nbkg, dNbkg;
    TCut Mwindow = Form("mlc > %lf && mlc < %lf",lend,rend);
    Ntot = ch1dat -> Draw("mlc>>+hdat","lcch==1 && tag==1 && abs(ml-1.1156)<0.002"+Mwindow,"goff"); //"lcch == 1 && ml>1.1 && ml<1.12"
    
    
    TF1* fdat = new TF1("fdat",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)+[3]+[4]*(x-2.3)",binw),lend,rend);
    TF1* fsig = new TF1("fsig",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)",binw),lend,rend);
    TF1* fbkg = new TF1("fbkg","[0]+[1]*(x-2.3)",lend,rend);    
    

    fdat -> SetParameters(1000,MLambdac,0.1,1400,-300);
    fdat -> SetParLimits(1,MLambdac-0.1,MLambdac+0.1);
    fdat -> SetParLimits(0,0,1e5);

    TFitResultPtr fitResult;
    fitResult = hdat -> Fit("fdat","L S M N","goff"); //L S M N Q
    TMatrixDSym covFit = fitResult->GetCovarianceMatrix();
    TMatrixDSym covSignal, covBackground;
    covFit.GetSub(0,2,covSignal);
    covFit.GetSub(3,4,covBackground);
    double * par;
    par = fdat -> GetParameters();
  
   
    fsig -> SetParameters(par);
    fbkg -> SetParameters(par+3);
    
    Nsig = fsig -> Integral(lend,rend)/binw; 
    dNsig = fsig -> IntegralError(lend,rend,par,covSignal.GetMatrixArray())/binw;
    Nbkg = fbkg -> Integral(lend,rend)/binw; 
    dNbkg = fbkg -> IntegralError(lend,rend,par+3,covBackground.GetMatrixArray())/binw;
    
    cout << "Total number of events: " << Ntot << endl;
    cout << "N Lambda_c: " << (int) Nsig <<" +- " << (int) dNsig << endl;
    cout << "N background: " << (int) Nbkg << " +- " << (int) dNbkg << endl;
 
    double chisq = fdat->GetChisquare(); 
    int ndf = Nbins - fdat ->GetNumberFreeParameters(); // nbins - npar
    cout << "ChiSquare / ndf: " << chisq << " / " << ndf << " = " << chisq/ndf << endl;
    
  
    
  
    hdat -> GetXaxis()-> SetTitle("M(#Lambda#pi) [GeV]");
    hdat -> GetXaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetXaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitle(Form("Events / ( %.2f )",binw));
    hdat -> GetYaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetYaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitleOffset(1);
    //hdat -> GetYaxis()->CenterTitle(true);
    hdat -> GetXaxis()->SetTickSize(0.04);
    hdat -> SetMarkerStyle(20);
    hdat -> SetMarkerSize(1.5);
    hdat -> SetMarkerColor(1);
    hdat -> SetLineColor(1);
    hdat -> SetLineWidth(2);
    hdat -> SetMinimum(0);
    hdat -> DrawCopy("e p");
    
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
    
    
    TLegend* leg = new TLegend(0.77,0.75,0.89,0.95);
    leg -> AddEntry("hdat","Data","E P");
   // leg->AddEntry("fsig","Signal","l");
    leg -> AddEntry("fdat","Fit","l");
    leg -> AddEntry("fbkg","Background","l");
    leg -> SetBorderSize(0);
    leg -> SetTextSize(axisFontSize);
    leg -> Draw("same");
    
    TLatex *tstatfit = new TLatex();
    tstatfit -> SetNDC();
    tstatfit -> SetTextColor(1);
    tstatfit -> SetTextSize(axisFontSize);
    tstatfit -> SetTextAngle(0);
    tstatfit -> DrawLatex(0.67, 0.65, Form("N_{sig} = %0.lf #pm %0.lf",Nsig, dNsig)); //
   // tstatfit -> DrawLatex(0.67, 0.59, Form("N_{bkg} = %0.lf #pm %0.lf",Nbkg, dNbkg)); //
    tstatfit -> DrawLatex(0.67, 0.59, Form("Mean_{sig} = %0.4lf #pm %0.4lf",par[1], fdat -> GetParError(1)));
    tstatfit -> DrawLatex(0.67, 0.53, Form("#sigma_{sig} = %0.4lf #pm %0.4lf",par[2], fdat -> GetParError(2)));
    
    delete hdat;

} 

void MergeExp(TString DataINPUT, TString DataOUTPUT)
{
    for(int i = 7; i<=73; i+=2)
    {   
        TChain chdat("h1"); 
        TString fn = Form("%d",i);
        if(i<10) fn = "0"+fn;
        chdat.Add(DataINPUT+"/"+fn+"*.root");
        chdat.Merge(DataOUTPUT+"/"+fn+".root");
    }
}

void SBpol2(double p3, double p4, double p5, double m, double s) //something is wrong!
{
    double a, b, c, A, B, C, D, p, q, Q, alpha, beta;
    
    a = p5;
    b = p3-2*p5*2.286;
    c = p4-2.286*p3+2.286*2.286*p5;
    
    A = 2./3*a;
    B = 6*a*s;
    C = 2*a*m*m+18*a*s*s+2*b*m+2*c;
    D = -a/3*(16*s*s*s+12*m*m*s)-4*b*m*s-4*s*c;
    
    p = (3*A*C-B*B)/(3*A*A);
    q = (2*B*B*B-9*A*B*C+27*A*A*D)/(27*A*A*A);
    
    Q = p*p*p/27+q*q/4;
    
    cout << "Q = " << Q << endl;
    
    alpha = -q/2+sqrt(Q);
    if (alpha < 0) alpha = -pow(-alpha,1./3);
        else alpha = pow(alpha,1./3);
    beta = -q/2-sqrt(Q);
    if (beta < 0) beta = -pow(-beta,1./3);
        else beta = pow(beta,1./3);
    
    cout << "d = " <<  alpha + beta << endl;
}


void Lctonumu() 
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
    
    double lend=2.21, rend=2.36, MLambdac=2.28646; //lend=2.21, rend=2.36
    int Nbins=50;
    TCanvas *c1 = new TCanvas("c1","Lambda_c invariant mass",1600,900);
    TH1D* hdat = new TH1D("hdat","#Lambda_{c} #rightarrow #Lambda#pi",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;
    
    
    double Ntot, Nsig, dNsig, Nbkg, dNbkg;
    TCut Mwindow = Form("mlc > %lf && mlc < %lf",lend,rend);
    Ntot = ch1dat -> Draw("mlc>>+hdat","lcch==1 && tag==1 && abs(ml-1.1156)<0.002"+Mwindow,"goff"); //"lcch == 1 && ml>1.1 && ml<1.12"
    
    
    TF1* fdat = new TF1("fdat",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)+[3]+[4]*(x-2.3)",binw),lend,rend);
    TF1* fsig = new TF1("fsig",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)",binw),lend,rend);
    TF1* fbkg = new TF1("fbkg","[0]+[1]*(x-2.3)",lend,rend);    
    

    fdat -> SetParameters(1000,MLambdac,0.1,1400,-300);
    fdat -> SetParLimits(1,MLambdac-0.1,MLambdac+0.1);
    fdat -> SetParLimits(0,0,1e5);

    TFitResultPtr fitResult;
    fitResult = hdat -> Fit("fdat","L S M N","goff"); //L S M N Q
    TMatrixDSym covFit = fitResult->GetCovarianceMatrix();
    TMatrixDSym covSignal, covBackground;
    covFit.GetSub(0,2,covSignal);
    covFit.GetSub(3,4,covBackground);
    double * par;
    par = fdat -> GetParameters();
  
   
    fsig -> SetParameters(par);
    fbkg -> SetParameters(par+3);
    
    Nsig = fsig -> Integral(lend,rend)/binw; 
    dNsig = fsig -> IntegralError(lend,rend,par,covSignal.GetMatrixArray())/binw;
    Nbkg = fbkg -> Integral(lend,rend)/binw; 
    dNbkg = fbkg -> IntegralError(lend,rend,par+3,covBackground.GetMatrixArray())/binw;
    
    cout << "Total number of events: " << Ntot << endl;
    cout << "N Lambda_c: " << (int) Nsig <<" +- " << (int) dNsig << endl;
    cout << "N background: " << (int) Nbkg << " +- " << (int) dNbkg << endl;
 
    double chisq = fdat->GetChisquare(); 
    int ndf = Nbins - fdat ->GetNumberFreeParameters(); // nbins - npar
    cout << "ChiSquare / ndf: " << chisq << " / " << ndf << " = " << chisq/ndf << endl;
    
  
    
  
    hdat -> GetXaxis()-> SetTitle("M(#Lambda#pi) [GeV]");
    hdat -> GetXaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetXaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitle(Form("Events / ( %.2f )",binw));
    hdat -> GetYaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetYaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitleOffset(1);
    //hdat -> GetYaxis()->CenterTitle(true);
    hdat -> GetXaxis()->SetTickSize(0.04);
    hdat -> SetMarkerStyle(20);
    hdat -> SetMarkerSize(1.5);
    hdat -> SetMarkerColor(1);
    hdat -> SetLineColor(1);
    hdat -> SetLineWidth(2);
    hdat -> SetMinimum(0);
    hdat -> DrawCopy("e p");
    
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
    
    
    TLegend* leg = new TLegend(0.77,0.75,0.89,0.95);
    leg -> AddEntry("hdat","Data","E P");
   // leg->AddEntry("fsig","Signal","l");
    leg -> AddEntry("fdat","Fit","l");
    leg -> AddEntry("fbkg","Background","l");
    leg -> SetBorderSize(0);
    leg -> SetTextSize(axisFontSize);
    leg -> Draw("same");
    
    TLatex *tstatfit = new TLatex();
    tstatfit -> SetNDC();
    tstatfit -> SetTextColor(1);
    tstatfit -> SetTextSize(axisFontSize);
    tstatfit -> SetTextAngle(0);
    tstatfit -> DrawLatex(0.67, 0.65, Form("N_{sig} = %0.lf #pm %0.lf",Nsig, dNsig)); //
   // tstatfit -> DrawLatex(0.67, 0.59, Form("N_{bkg} = %0.lf #pm %0.lf",Nbkg, dNbkg)); //
    tstatfit -> DrawLatex(0.67, 0.59, Form("Mean_{sig} = %0.4lf #pm %0.4lf",par[1], fdat -> GetParError(1)));
    tstatfit -> DrawLatex(0.67, 0.53, Form("#sigma_{sig} = %0.4lf #pm %0.4lf",par[2], fdat -> GetParError(2)));
    
    delete hdat;

} 
