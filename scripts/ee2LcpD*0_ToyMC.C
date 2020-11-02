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
    double dphi, scalprod, ppi;
    TFile *f = new TFile("ToyLc3bDecay.root", "recreate");
    TTree *tree = new TTree("h1", "Toy MC ROOT tree");
    tree->Branch("dphi", &dphi, "dphi/D");
    tree->Branch("scalprod", &scalprod, "scalprod/D");
    tree->Branch("ppi", &ppi, "ppi/D");
    
    
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
    
    for (Int_t n=0;n<1e6;n++) 
    {
        ee2XLcEvent.SetDecay(pUPS, 3, InitMasses);
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
        
        tree -> Fill();
    }
    
    tree -> Write();
    f->Close();
    
}

TLorentzVector boostT( TLorentzVector p, TLorentzVector p_boost)
{                                                                
    TLorentzVector tmp = p;
    tmp.Boost(-p_boost.BoostVector());
    return tmp;
}
