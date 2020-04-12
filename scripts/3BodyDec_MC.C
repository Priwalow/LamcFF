#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <TMath.h>
#include <TRandom.h>


//three body decay
    
 double M_L = 1.11568, M_Lc = 2.28646, //GeV 
        M_p = 0.93827, M_K = 0.49768, 
        M_pi = 0.13957, M_mu = 0.10566,
        M_nu = 0,       M_pi0 = 0.13498,
        M_Xi0 = 1.31486, M_Xic = 2.46793;        
 
 double c = 29979245800, ctau_L=7.89,                         //ctau in cm,
                         ctau_Xi0 = 8.69;


 void PhaseSpace() {
 
    double pnu, pvis, mlmu, qsq, mx, px, ldec, rhodec, 
           rmx, rmvis, plc, mvis, ptmu, ptl; 
    int npi, npi0;
    TFile *f = new TFile("ToyLc3bDecay.root", "recreate");
    TTree *tree = new TTree("h1", "Toy MC ROOT tree");
    tree->Branch("pnu", &pnu, "pnu/D");
    tree->Branch("ptmu", &ptmu, "ptmu/D");
    tree->Branch("ptl", &ptl, "ptl/D");
    tree->Branch("mlmu", &mlmu, "mlmu/D");
    tree->Branch("plc",&plc,"plc/D"); //plc - p(Lambda+lepton) in CM
    tree->Branch("qsq", &qsq, "qsq/D");
    tree->Branch("npi", &npi, "npi/I");
    tree->Branch("npi0", &npi0, "npi0/I");
    tree->Branch("mx", &mx, "mx/D");
    tree->Branch("px", &px, "px/D");
    tree->Branch("pvis", &pvis, "pvis/D");
    tree->Branch("rmx", &rmx, "rmx/D");
    tree->Branch("rmvis", &rmvis, "rmvis/D");
    tree->Branch("ldec",&ldec,"ldec/D");
    tree->Branch("rhodec",&rhodec,"rhodec/D");
    
    TRandom r; 
    r.SetSeed(0);
    
    double E_e = 8.0, E_p = 3.5;
    TLorentzVector pUPS= TLorentzVector (E_e*sin(0.022), 0.,E_e*cos(0.022)-E_p, E_e+E_p);
    TVector3 CMsystBoost = pUPS.BoostVector();
    pUPS.Boost(-CMsystBoost);
    //cout << pUPS.Px() << "  " << pUPS.Py() << "  " << pUPS.Pz() << "    " << pUPS.E() << endl;
    
    //(Momentum, Energy units are Gev/C, GeV)
   

    double Lc_dec_prod_masses[3] = {M_L,M_mu,M_nu};  
    
    TGenPhaseSpace XLcEvent, LcDecEvent;
    double* masses = NULL;
   
    for (Int_t n=0;n<100000;n++) {
       
        npi = TMath::Nint(2*(r.Poisson(2.42))+1);
        npi0  = TMath::Nint(r.Exp(2.44));
        
        masses = new double [npi+npi0+3]; 
        masses[0] = M_Lc;
        masses[1] = M_p;
        masses[2] = M_K;
        for ( int i = 3; i< npi+3; i++)
            masses[i]=M_pi;
        for ( int i = npi+3; i< npi+npi0+3; i++)
            masses[i]=M_pi0;
        
        if (!XLcEvent.SetDecay(pUPS, npi+npi0+3, masses))
        {
            n--;
            continue;
        };
        
        Double_t weight = XLcEvent.Generate();
 
       TLorentzVector pLc = *XLcEvent.GetDecay(0);
       TLorentzVector pX(0,0,0,0);
       for(int i = 1; i<npi+3; i++)
           pX+=*XLcEvent.GetDecay(i);
       
       LcDecEvent.SetDecay(pLc,3,Lc_dec_prod_masses);
       LcDecEvent.Generate();
       TLorentzVector pL =  *LcDecEvent.GetDecay(0);
       TLorentzVector pMu = *LcDecEvent.GetDecay(1);
       TLorentzVector pNu = *LcDecEvent.GetDecay(2);
       
       
       pnu = pNu.P();
       pvis = (pL+pMu+pX).P();
       mlmu = (pL+pMu).M();
       qsq = (pNu+pMu).M2();
       mx = pX.M();
       px = pX.P();
       rmx = (pUPS-pX).M();
       rmvis = (pUPS-pX-pL-pMu).M();
       plc = (pL+pMu).P();
       ptmu = pMu.P()*sin(pMu.Theta());
       ptl = pL.P()*sin(pL.Theta());
       
       pL.Boost(CMsystBoost);
       ldec = r.Exp(ctau_L*pL.Gamma()*pL.Beta());
       rhodec = sin(pL.Theta())*ldec;
       tree -> Fill();
       
       delete [] masses;
    }
   tree -> Write();
   f->Close();
   
 }
