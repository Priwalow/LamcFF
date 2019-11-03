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

double dE=0.01, dtheta=0.01, dphi=0.01;  // in %
double M=2.286, m_1=1.116, m_2=0.105, m_3=1e-15; // m_1 + m_2 + m_3 < M !!! m_1 > m_2 > m_3
double Pi = TMath::Pi();



inline double sq(double x)
{
    return x*x;
}


void p2calc(double p1, double & p2, double & ct)
{
    double A, B, C, D, e1 = sqrt(p1*p1+m_1*m_1), p2p, p2m, minabscos12;
    
    
    C = 0.5*(M*M+m_1*m_1+m_2*m_2-2*e1*M-m_3*m_3);
    A = C/(M-e1);
    B = p1/(M-e1);
    minabscos12 = (m_2*m_2-A*A)/sq(m_2*B);
    
    
    D=-1;
    TRandom r; 
    r.SetSeed(0);
    do {
        if (minabscos12<0) 
            ct = r.Uniform(-1,1);
        else
            ct = (2*r.Binomial(1,0.5)-1)*r.Uniform(minabscos12,1);
      D = (sq(m_2*B*ct)-m_2*m_2+A*A);
    } while (D<0);
    
    p2p = (sqrt(D)-A*B*ct)/(1-B*ct*B*ct);
    if (p2p<0)
    {
       // cout << "p2p < 0, correcting cos12 -> -cos12" << endl;
        ct = -ct;
        p2p = (sqrt(D)-A*B*ct)/(1-B*ct*B*ct);
    }
    p2=p2p;
}


void test(unsigned int Nevtot)
{
    TFile *f = new TFile("Nu_MC.root", "recreate");
    TTree *tmc = new TTree("mc", "True values");
    TTree *td = new TTree("data", "Values with detector response");
    
    double rm, m12, m23;
    
    tmc->Branch("m12", &m12, "m12/D");
    tmc->Branch("m23", &m23, "m23/D");
    tmc->Branch("rm", &rm, "rm/D");
    
    td->Branch("m12", &m12, "m12/D");
    td->Branch("m23", &m23, "m23/D");
    td->Branch("rm", &rm, "sqrm/D");
    
    
    TLorentzVector P_0, P_1, P_2, P_3; //P_0 = P_1 + P_2 + P_3
    
    P_0.SetXYZM(0,0,0,M);
    
    TRandom r; 
    r.SetSeed(0);
    
    double p_1_max, E_1_max = (sq(M)+sq(m_1)-sq(m_3+m_2))/(2*M);
    p_1_max = sqrt(sq(E_1_max)-sq(m_1));
    
    TVector3 p_1, p_2, rndmAxis;
    double E_1, E_2, pmag_2, pmag_1, rndmAng, minabscos12, cos12; // angle between momenta of particles 1 and 2
    
    
    
    for(int i =0; i<Nevtot; i++)
    {
        pmag_1 = r.Uniform(0,p_1_max); //Set particle 1
        P_1.SetXYZM(0,0,pmag_1,m_1); 
        E_1 = P_1.E(); 
       
       p2calc(pmag_1, pmag_2, cos12); //Set particle 2
       P_2.SetXYZM(0,0,pmag_2,m_2);
       P_2.RotateY(acos(cos12)); 
       
       r.Sphere(rndmAxis(0),rndmAxis(1),rndmAxis(2),1); // Rotate plane of particles 1 and 2
       rndmAng =  r.Uniform(2*Pi);
       P_1.Rotate(rndmAng,rndmAxis);
       P_2.Rotate(rndmAng,rndmAxis);
    
       
       P_3 = P_0 - P_1 - P_2; //Obtain particle 3
    
       rm = P_3.M(); 
       m12 = (P_1+P_2).M2();
       m23 = (P_0-P_1).M2();
       tmc -> Fill(); //Fill tree with real data
       
       
       
       P_1.SetPhi(P_1.Phi()*(1+r.Gaus(0,dphi)));
       P_1.SetTheta(P_1.Theta()*(1+r.Gaus(0,dtheta)));
       P_1.SetE(P_1.E()*(1+r.Gaus(0,dE)));
       
       P_2.SetPhi(P_2.Phi()*(1+r.Gaus(0,dphi)));
       P_2.SetTheta(P_2.Theta()*(1+r.Gaus(0,dtheta)));
       P_2.SetE(P_2.E()*(1+r.Gaus(0,dE)));
       
       P_3 = P_0 - P_1 - P_2;
       
       rm = P_3.M(); 
       m12 = (P_1+P_2).M2();
       m23 = (P_0-P_1).M2();
       td -> Fill(); //Fill tree with real data
    }
    tmc -> Write();
    td -> Write();
    f->Write();
    f->Close();
    
}
