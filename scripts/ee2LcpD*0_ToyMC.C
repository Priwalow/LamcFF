#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <TMath.h>
#include <TRandom.h>

TLorentzVector boostT( TLorentzVector p, TLorentzVector p_boost);
double heli (  TLorentzVector   P4A, TLorentzVector  P4B, TLorentzVector P4system);
//three body decay

double M_L = 1.11568, M_Lc = 2.28646, //GeV 
M_p = 0.93827, M_K = 0.49768, 
M_pi = 0.13957, M_mu = 0.10566,
M_nu = 0,       M_pi0 = 0.13498,
M_Xi0 = 1.31486, M_Xic = 2.46793,
M_Dp = 1.86965,  M_Dst0 = 2.00685;     

double c = 29979245800, ctau_L=7.89,                         //ctau in cm,
ctau_Xi0 = 8.69,
ctau_Dp = 0.0312;

void PhaseSpace() 
{
    double dphi, scalprod, ppi, hlc, hl;
    TFile *f = new TFile("ToyLc3bDecay.root", "UPDATE");
    TTree *tree = new TTree("h2", "Phase space");
    tree->Branch("dphi", &dphi, "dphi/D");
    tree->Branch("scalprod", &scalprod, "scalprod/D");
    tree->Branch("ppi", &ppi, "ppi/D");
    tree->Branch("hlc", &hlc, "hlc/D");
    tree->Branch("hl", &hl, "hl/D");
    
    TRandom r; 
    r.SetSeed(0);
    
    double E_e = 8.0, E_p = 3.5;
    TLorentzVector pUPS= TLorentzVector (E_e*sin(0.022), 0.,E_e*cos(0.022)-E_p, E_e+E_p);
    TVector3 CMsystBoost = pUPS.BoostVector();
    //pUPS.Boost(-CMsystBoost);
    
    double InitMasses[3] = {M_Lc,M_Dst0,M_p},
    Lc_dec_prod_masses[2] = {M_L,M_pi}, 
    L_dec_prod_masses[2] = {M_p,M_pi};  
    
    TGenPhaseSpace ee2XLcEvent, LcDecEvent, LDecEvent;
    TLorentzVector P_Lc, P_L, P_p, P_pi;
    ee2XLcEvent.SetDecay(pUPS, 3, InitMasses);
    
    for (Int_t n=0;n<1e8;n++) 
    {
        
        ee2XLcEvent.Generate();
        
        P_Lc = *ee2XLcEvent.GetDecay(0);
        LcDecEvent.SetDecay(P_Lc,2,Lc_dec_prod_masses);
        LcDecEvent.Generate();
        
        P_L = *LcDecEvent.GetDecay(0);
        LDecEvent.SetDecay(P_L,2,L_dec_prod_masses);     
        LDecEvent.Generate();
        
        P_p = *LDecEvent.GetDecay(0);
        P_pi = *LDecEvent.GetDecay(1);
        
        TVector3 norm_lam_c, norm_lam;
        norm_lam_c = (boostT(P_L, P_Lc).Vect()).Cross(boostT(pUPS, P_Lc).Vect());
        norm_lam_c = norm_lam_c*(1./norm_lam_c.Mag());
        norm_lam = (boostT(P_p, P_Lc).Vect()).Cross(boostT(P_pi, P_Lc).Vect());
        norm_lam = norm_lam*(1./norm_lam.Mag());
        dphi=norm_lam_c.Angle(norm_lam);
        scalprod = boostT(P_p, P_Lc).Vect().Dot(norm_lam_c);
        if(scalprod < 0) dphi = 2*TMath::Pi()-dphi;
        
        ppi = P_pi.Vect().Mag();
        
        
        hlc = -cos(heli(P_L, pUPS,  P_Lc));
        hl = -cos(heli(P_p, P_Lc,  P_L));
        
        tree -> Fill();
    }
    
    tree -> Write();
    f->Close();
    
}

void AngDist()
{
    double dphi, hlc, hl, plccms, clccms;
    TFile *f = new TFile("ToyLc3bDecay.root", "recreate");
    TTree *tree = new TTree("h1", "Toy MC ROOT tree");
    tree->Branch("philclam", &dphi, "philclam/D");
    tree->Branch("hlc", &hlc, "hlc/D");
    tree->Branch("hl", &hl, "hl/D");
    tree->Branch("plamccms", &plccms, "plamccms/D");
    tree->Branch("coslccms", &clccms, "coslccms/D");
    
    TTree *t2 = new TTree("h2", "Toy MC ROOT tree");
    t2->Branch("philclam", &dphi, "philclam/D");
    t2->Branch("hlc", &hlc, "hlc/D");
    t2->Branch("hl", &hl, "hl/D");
    t2->Branch("plamccms", &plccms, "plamccms/D");
    t2->Branch("coslccms", &clccms, "coslccms/D");
    
    TRandom r; 
    r.SetSeed(0);
    double P = -0.7, delta = -0.5*TMath::Pi(), alphaLam=0.732, alphaLam_c=-0.84;
    
    //alphaLam = r.Gaus(0.732,0.014);
    //alphaLam_c = r.Gaus(-0.84,0.09);
    
    TF3 fmain("fmain",Form("1+[0]*y+%lf*(z*([1]+[2]*y)-sqrt((1-y*y)*(1-z*z))*[3]*cos(%lf+x))",P,delta),-TMath::Pi(),TMath::Pi(),-1,1,-1,1);
    fmain.SetParameters(alphaLam*alphaLam_c,alphaLam_c,alphaLam,alphaLam*sqrt(1-alphaLam_c*alphaLam_c));
    
    
    for (Int_t n=0;n<1e2;n++) 
    {

        fmain.GetRandom3(dphi,hl,hlc);
        plccms = r.Uniform(0,5);
        clccms = r.Uniform(-1,1);
        tree -> Fill();
        t2->Fill();
        
    }
    
    
    tree -> Write();
    t2 -> Write();
    f->Close();
    
    
    TChain* ch1dat = new TChain("h1");
    ch1dat -> Add("ToyLc3bDecay.root");
    
    double lend=-1, rend=1; //lend=2.21, rend=2.36
    int Nbins=10;
    TH3D* hmain = new TH3D("hmain","#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}",3,-TMath::Pi(),TMath::Pi(),2,lend,rend,2,lend,rend);
    ch1dat -> Draw("hlc:hl:philclam>>hmain","","goff"); 
    
    TF3* fdat = new TF3("fdat",Form("[0]*(1+%lf*y+[1]*(z*(%lf+%lf*y)-sqrt((1-y*y)*(1-z*z))*%lf*cos([2]+x)))",
                                    alphaLam*alphaLam_c,alphaLam_c,alphaLam,alphaLam*sqrt(1-alphaLam_c*alphaLam_c)),-TMath::Pi(),TMath::Pi(),lend,rend,lend,rend);
    fdat -> SetParameters(4,0.4,3);
    fdat -> SetParLimits(2,-TMath::Pi()-0.1,TMath::Pi()+0.1);
    TFitResultPtr fitResult;
    fitResult = hmain -> Fit("fdat","L S M N","goff");
    double chisq = fdat->GetChisquare(), ndf = fdat->GetNDF(); 
    cout << "ChiSquare / ndf: " << chisq << " / " << ndf << " = " << chisq/ndf << "; Prob = " << fdat -> GetProb()<<  endl;
    
}


TLorentzVector boostT( TLorentzVector p, TLorentzVector p_boost)
{                                                                
    TLorentzVector tmp = p;
    tmp.Boost(-p_boost.BoostVector());
    return tmp;
}


double heli (  TLorentzVector   P4A, TLorentzVector  P4B, TLorentzVector P4system) // return angle btw P4A and P4B in P4system rest frame
{
        return (boostT(P4A, P4system).Vect()).Angle(boostT(P4B, P4system).Vect());
}
