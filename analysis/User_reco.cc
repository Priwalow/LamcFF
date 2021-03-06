#include "User_reco.h"
#include "userinfo.h"
#include EVTCLS_H //R2 distribution

//e+ e- -> x+Lc

using namespace Belle;
namespace Belle {
    void User_reco::hist_def( void )
    {
        
        extern BelleTupleManager* BASF_Histogram;
        t1 = BASF_Histogram->ntuple ("data","tag dch dstch md mdst rmx rmvis rmvis_nopi0 mvis px plamclab plamccms coslclab coslccms pvis fox ecms mks ch_tag lcch ml mlc hlc hl q hw chi lcp2dcm lcp2dlab philclam plslc addpi addpi0 mcwadpi0 totcharg" ); // not ALL momenta in CMS! 	lepton cosTheta in CMS, rholam, rholamcms	
        //t2 = BASF_Histogram->ntuple ("withGamma","lcch tag ml mlc mx mvis npi npi0 ngamma ecms egammatot rmx rmvis plc px pvis ml1 hl hlc phi q fox hw chi" ); // ALL momenta in CMS! 
        //"tag dch dstch mlc ml md rmx rmvis px npi npi0 nk ngam fox pvis ecms hlc hl hw chi q lcch"
    };
    //***********************************************************************
    void doMassVertexFit(class vector<Particle> &p_list, double mass=-1);
    void doMassVertexFit(class Particle &P, double mass=-1);
    void doVertexFit(class vector<Particle> &p_list);
    void doVertexFit(class Particle &P);
    void setPi0Error(Particle &p);
    void setPi0Error(std::vector<Particle> &p_list);
    void setGammaError(Particle &p, const HepPoint3D &gVertex, const HepSymMatrix &errGVertex);
    void setGammaError(std::vector<Particle> &p_list);
    double heli (  HepLorentzVector  * P4A, HepLorentzVector * P4B, HepLorentzVector* P4system);
    double heli (  HepLorentzVector   P4A, HepLorentzVector P4B, HepLorentzVector P4system);
    HepLorentzVector boostT( HepLorentzVector p, HepLorentzVector p_boost);
    HepLorentzVector boostT( HepLorentzVector *p, HepLorentzVector *p_boost);
    void eraseLeptons(std::vector<Particle> &list);
    
    void withdRdZcut(class std::vector<Particle> &p_list,double ip_position=0., double drcut = 2., double dzcut = 4.);
    void makeLam(std::vector<Particle> &lam0, std::vector<Particle> &lam0b);
    void makeProton (std::vector<Particle> &p, std::vector<Particle> &pbar);
    
    int calcuCharge (Particle *p );
    
    void User_reco::event ( BelleEvent* evptr, int* status )
    {
        static int nevent=0,skimmed=0,skimmedPi0=0, skimmed2=0, skimmed3=0, skimmed4=0;
        nevent++;
        
        
        *status=0;
        
        
        Belle_runhead_Manager& rhdmgr = Belle_runhead_Manager::get_manager();
        Belle_runhead_Manager::const_iterator rhd = rhdmgr.begin();
        
        
        double elec=8.0, posi=3.5;
        if(rhd!=rhdmgr.end() && rhd->EHER() > 5.){
            elec = rhd->EHER();
            posi = rhd->ELER();
        }
        HepLorentzVector pUPS= HepLorentzVector (elec*sin(0.022), 0.,elec*cos(0.022)-posi, elec+posi);
        
        Evtcls_hadron_info_Manager&  ehimgr = Evtcls_hadron_info_Manager::get_manager();
        
        
        std::vector<Evtcls_hadron_info>::iterator iti = ehimgr.begin();
        const HepPoint3D&   runIp     = IpProfile::position();
        
        double fox=0;
        int ntrk = 0;
        if( iti !=  ehimgr.end() && *iti )
        {
            fox = (*iti).R2();
            //       ntrk = (*iti).Ntrk();
        }
        std::cout << "INITIALIZING!" << std::endl;
        
        
        //------------------------MAKE PARTICLE LISTINGS----------------------------------
        //protons
        std::vector<Particle> p_p, p_m; 
        makeProton(p_p,p_m);
        
        if(p_p.size()+p_m.size()==0) return; 
        
        
        //kaons and pions
        std::vector<Particle>  k_p, k_m, pi_p, pi_m, pions;
        makeKPi(k_p, k_m, pi_p, pi_m,1);

        
        ntrk=k_p.size()+k_m.size();
        
        withdRdZcut(k_p,runIp.z());
        withdRdZcut(pi_p,runIp.z());
        withdRdZcut(k_m,runIp.z());
        withdRdZcut(pi_m,runIp.z());
        
        for(std::vector<Particle>::iterator l = pi_m.begin(); l!=pi_m.end(); ++l)
            pions.push_back(*l);
        for(std::vector<Particle>::iterator l = pi_p.begin(); l!=pi_p.end(); ++l)
            pions.push_back(*l);
        
        withKaonId(k_p,0.6,3,1,5,3,4);
        withKaonId(k_p,0.6,3,1,5,3,2);
        withKaonId(k_m,0.6,3,1,5,3,4);
        withKaonId(k_m,0.6,3,1,5,3,2);
        
        
        
        
        /*for(std::vector<Particle>::iterator l = pi_m.begin(); l!=pi_m.end(); ++l)
            pions.push_back(*l);
        for(std::vector<Particle>::iterator l = pi_p.begin(); l!=pi_p.end(); ++l)
            pions.push_back(*l);*/
        
        //Ks mesons
        std::vector<Particle> k_s;
        makeKs(k_s);
        for(std::vector<Particle>::iterator l = k_s.begin(); l!=k_s.end(); ++l)
        {
            HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
            Vector3 P(l->px(),l->py(),0);
            V=V-runIp;
            V.setZ(0.);
            if (abs(l->mass()-0.4976)>0.015 || V.perp()<0.1 ||
                V.angle(P)>0.01 || l->mdstVee2().z_dist()>10. ) {
                k_s.erase(l); --l;}
        }
        doMassVertexFit(k_s);
        
        /*double k_s_chisq, bufchisq=1000000;
         *        while(k_s.size()>1)
         *        {
         *            for(std::vector<Particle>::iterator l = k_s.begin(); l!=k_s.end(); ++l)
         *            {
         *                k_s_chisq = dynamic_cast<UserInfo&>(l->userInfo()).vchisq();
         *                if((k_s_chisq<0) || (k_s_chisq > bufchisq))
         *                {
         *                    k_s.erase(l); 
         *                    --l;
         *                    continue;
    }
    bufchisq = k_s_chisq;
    }
    }
    
    if (!(nevent%1000))std::cout<<nevent<<"Best k_s candidate's chisq/ndf = " << k_s_chisq << '\n';
    */
        //Pi0 mesons
        std::vector<Particle>  pi0;
        makePi0(pi0);
        
        for(std::vector<Particle>::iterator i=pi0.begin(); i!=pi0.end();++i)
            if(i->mdstPi0().gamma(0).ecl().energy()<0.05||
                i->mdstPi0().gamma(1).ecl().energy()<0.05||
                abs(i->mdstPi0().mass()-.135)>0.015)
            {pi0.erase(i); --i;}
        setPi0Error(pi0);        
        
        
        //photons
        /*std::vector<Particle> photons;
        Mdst_gamma_Manager& GamMgr = Mdst_gamma_Manager::get_manager();
        for (std::vector<Mdst_gamma>::iterator it=GamMgr.begin();it!=GamMgr.end(); it++) 
        {
            Particle prtcl(*it);
            bool gam_from_pi0 = false;
            for (std::vector<Particle>::iterator pi=pi0.begin(); pi!=pi0.end();++pi)
            {
                if (checkSame(*it,*pi))
                {
                    gam_from_pi0 = true;
                    break;
                }
            }
            if  ( (prtcl.e()>0.05) && (!gam_from_pi0) )
            {
                photons.push_back(prtcl);
            }
        }*/
        
        //Lambda hyperons
        std::vector<Particle> lam, lamb;
        makeLam(lam,lamb);
        
        setUserInfo(lam,  11 ); 
        setUserInfo(lamb, 12 ); 
        doMassVertexFit(lam);
        doMassVertexFit(lamb);
        
        //leptons (e, mu) 
        std::vector<Particle>  e_p,e_m,mu_p,mu_m;
        makeLepton(e_p,e_m,mu_p,mu_m);
        
        withMuId(mu_p);
        withMuId(mu_m);
        withEId(e_p);
        withEId(e_m);
        
        
        std::vector<Particle> D0, D0_b, D_p, D_m, Dst_p, Dst_m, Dst0, Dst0_b;
        
        //######################################    TAG SIDE
    
        combination (D0_b,ptype_D0B, k_p, pi_m, 0.04);
        combination (D0,ptype_D0, k_m, pi_p, 0.04);
        setUserInfo(D0_b, 1); 
        setUserInfo(D0, 1); 
        combination (D0_b,ptype_D0B, pi_m, pi_m, k_p, pi_p, 0.04);
        combination (D0,ptype_D0, pi_p, pi_p, k_m, pi_m, 0.04);
        setUserInfo(D0_b, 2); 
        setUserInfo(D0, 2); 
        combination (D0_b,ptype_D0B, k_s, pi_p, pi_m, 0.04);
        combination (D0,ptype_D0, k_s, pi_p, pi_m, 0.04);
        setUserInfo(D0_b, 3);
        setUserInfo(D0, 3); 
        combination (D0_b,ptype_D0B, k_p, pi0, pi_m,  0.04);
        combination (D0,ptype_D0, k_m, pi0, pi_p,  0.04);
        setUserInfo(D0_b, 4);
        setUserInfo(D0, 4); 
        
        //combination (D0_b,ptype_D0B, k_p, pi0, pi_m, pi_p, pi_m, 0.06);
        //combination (D0,ptype_D0, k_m, pi0, pi_p, pi_p, pi_m, 0.06);
        //setUserInfo(D0_b, 5);
        //setUserInfo(D0, 5);
        combination (D0_b,ptype_D0B, k_s, pi0, pi_p, pi_m, 0.04);
        combination (D0,ptype_D0, k_s, pi0, pi_p, pi_m, 0.04);
        setUserInfo(D0_b, 6);
        setUserInfo(D0, 6);
        
        combination (D_m,ptype_Dm, k_p, pi_m, pi_m, 0.04);
        combination (D_p,ptype_Dp, k_m, pi_p, pi_p, 0.04);
        setUserInfo(D_m, 1);
        setUserInfo(D_p, 1);
        combination (D_m,ptype_Dm, k_s, pi_m, 0.04);
        combination (D_p,ptype_Dp, k_s, pi_p, 0.04);
        setUserInfo(D_m, 2);
        setUserInfo(D_p, 2);
        combination (D_m,ptype_Dm, pi_m, pi_m, k_s, pi_p, 0.04);
        combination (D_p,ptype_Dp, pi_p, pi_p, k_s, pi_m,0.04);
        setUserInfo(D_m, 3);
        setUserInfo(D_p, 3);
        combination (D_m,ptype_Dm, k_p, pi_m, k_m, 0.04);
        combination (D_p,ptype_Dp, k_m, pi_p, k_p, 0.04);
        setUserInfo(D_m, 4);
        setUserInfo(D_p, 4);
        
        if(D0_b.size()+D_m.size()+D0.size()+D_p.size()==0) return;
        
        doMassVertexFit(D0_b);
        doMassVertexFit(D_m);
        doMassVertexFit(D0);
        doMassVertexFit(D_p);        
        
        
        combination (Dst0_b,ptype_DstB, D0_b, pi0, 0.04);
        setUserInfo(Dst0_b, 1);
        //combination (Dst0_b,ptype_DstB, D0_b, photons, 0.025);
        //setUserInfo(Dst0_b, 2);
        
        for (std::vector<Particle>::iterator i=Dst0_b.begin(); i!=Dst0_b.end();++i)
        {
            Particle dst0b(*i);
            if(abs(dst0b.mass()-dst0b.child(0).mass()-0.142)>0.01)
            {
                Dst0_b.erase(i); 
                --i;
            }
        }        
        
        combination (Dst_m,ptype_Dstm, D0_b, pi_m, 0.04);
        setUserInfo(Dst_m, 1);
        combination (Dst_m,ptype_Dstm, D_m, pi0, 0.04);
        setUserInfo(Dst_m, 2);
            
        for (std::vector<Particle>::iterator i=Dst_m.begin(); i!=Dst_m.end();++i)
        {
            Particle dstm(*i);
            if(abs(dstm.mass()-dstm.child(0).mass()-0.1425)>0.01)
            {
                Dst_m.erase(i); 
                --i;
            }
        }
        
        combination (Dst0,ptype_Dst0, D0, pi0, 0.04);
        setUserInfo(Dst0, 1);
        //combination (Dst0,ptype_Dst0, D0, photons, 0.025);
        //setUserInfo(Dst0, 2);
        
        for (std::vector<Particle>::iterator i=Dst0.begin(); i!=Dst0.end();++i)
        {
            Particle dst0(*i);
            if(abs(dst0.mass()-dst0.child(0).mass()-0.142)>0.01)
            {
                Dst0.erase(i); 
                --i;
            }
        }
        
        combination (Dst_p,ptype_Dstp, D0, pi_p, 0.04);
        setUserInfo(Dst_p, 1);
        combination (Dst_p,ptype_Dstp, D_p, pi0, 0.04);
        setUserInfo(Dst_p, 2);
 
        for (std::vector<Particle>::iterator i=Dst_p.begin(); i!=Dst_p.end();++i)
        {
            Particle dstp(*i);
            if(abs(dstp.mass()-dstp.child(0).mass()-0.1425)>0.01)
            {
                Dst_p.erase(i); 
                --i;
            }
        }
        
        
        doMassVertexFit(Dst0_b);
        doMassVertexFit(Dst_m);
        doMassVertexFit(Dst0);
        doMassVertexFit(Dst_p);
        
        
        
        
        std::vector <Particle> L_, L_b;
        combination (L_b,ptype_Lamc, p_m, D0_b);
        setUserInfo(L_b, 1);
        //combination (L_b,ptype_Lamc, p_m, D0_b, pi_p, pi_m);
        //setUserInfo(L_b, 11);
        //combination (L_b,ptype_Lamc, p_m, D0_b, pi0);
        //setUserInfo(L_b, 12);
        combination (L_b,ptype_Lamc, p_m, D_m, pi_p);
        setUserInfo(L_b, 2);
        //combination (L_b,ptype_Lamc, p_m, D_m, pi_p, pi0);
        //setUserInfo(L_b, 22);
        combination (L_b,ptype_Lamc, p_m, Dst_m, pi_p);
        setUserInfo(L_b, 3);
        //combination (L_b,ptype_Lamc, p_m, Dst_m, pi_p, pi0);
        //setUserInfo(L_b, 32);
        combination (L_b,ptype_Lamc, p_m, Dst0_b);
        setUserInfo(L_b, 4);
        //combination (L_b,ptype_Lamc, p_m, Dst0_b, pi_p, pi_m);
        //setUserInfo(L_b, 41);
        //combination (L_b,ptype_Lamc, p_m, Dst0_b, pi0);
        //setUserInfo(L_b, 42);
        
        combination (L_,ptype_ALamc, p_p, D0);
        setUserInfo(L_, 1);
        //combination (L_,ptype_ALamc, p_p, D0, pi_m, pi_p);
        //setUserInfo(L_, 11);
        //combination (L_,ptype_ALamc, p_p, D0, pi0);
        //setUserInfo(L_, 12);
        combination (L_,ptype_ALamc, p_p, D_p, pi_m);
        setUserInfo(L_, 2);
        //combination (L_,ptype_ALamc, p_p, D_p, pi_m, pi0);
        //setUserInfo(L_, 22);
        combination (L_,ptype_ALamc, p_p, Dst_p, pi_m);
        setUserInfo(L_, 3);
        //combination (L_,ptype_ALamc, p_p, Dst_p, pi_m, pi0);
        //setUserInfo(L_, 32);
        combination (L_,ptype_ALamc, p_p, Dst0);
        setUserInfo(L_, 4);
        //combination (L_,ptype_ALamc, p_p, Dst0, pi_m, pi_p);
        //setUserInfo(L_, 41);
        //combination (L_,ptype_ALamc, p_p, Dst0, pi0);
        //setUserInfo(L_, 42);
        
        
        //######################################    SIGNAL SIDE
        
        std::vector <Particle> Lc, Lcb; 
        
        combination (Lc,ptype_Lamc, lam, e_p);
        combination (Lcb,ptype_Lamc, lamb, e_m);
        setUserInfo(Lc,3);
        setUserInfo(Lcb,3);  
        combination (Lc,ptype_Lamc, lam, e_m);     //fake
        combination (Lcb,ptype_Lamc, lamb, e_p);   //fake
        setUserInfo(Lc,300);
        setUserInfo(Lcb,300);
        
        combination (Lc,ptype_Lamc, lam, mu_p);
        combination (Lcb,ptype_Lamc, lamb, mu_m);
        setUserInfo(Lc,4);
        setUserInfo(Lcb,4);
        combination (Lc,ptype_Lamc, lam, mu_m);     //fake
        combination (Lcb,ptype_Lamc, lamb, mu_p);   //fake
        setUserInfo(Lc,400);
        setUserInfo(Lcb,400);
        
        combination (Lc,ptype_Lamc, lam, e_p, pi0);
        combination (Lcb,ptype_Lamc, lamb, e_m, pi0);
        setUserInfo(Lc,6);
        setUserInfo(Lcb,6);  
        
        combination (Lc,ptype_Lamc, lam, mu_p, pi0);
        combination (Lcb,ptype_Lamc, lamb, mu_m, pi0);
        setUserInfo(Lc,7);
        setUserInfo(Lcb,7);
       
     /*
        combination (Lc,ptype_Lamc, lam, e_p, pi_p, pi_m);
        combination (Lcb,ptype_Lamc, lamb, e_m, pi_p, pi_m);
        setUserInfo(Lc,8);
        setUserInfo(Lcb,8);  
        
        combination (Lc,ptype_Lamc, lam, mu_p, pi_p, pi_m);
        combination (Lcb,ptype_Lamc, lamb, mu_m, pi_p, pi_m);
        setUserInfo(Lc,9);
        setUserInfo(Lcb,9);
        */
        
        for(std::vector<Particle>::iterator l = Lc.begin(); l!=Lc.end(); ++l)
            if (l->mass()>ptype_Lamc.mass()+0.1)
                {Lc.erase(l); --l;}
        for(std::vector<Particle>::iterator l = Lcb.begin(); l!=Lcb.end(); ++l)
            if (l->mass()>ptype_Lamc.mass()+0.1)
                {Lcb.erase(l); --l;}
                
                
        combination (Lc,ptype_Lamc, lam, pi_p,0.05);
        combination (Lcb,ptype_Lamc, lamb, pi_m,0.05);
        setUserInfo(Lc,1);
        setUserInfo(Lcb,1);
        
        combination (Lc,ptype_Lamc, lam, pi_p, pi0, 0.08);
        combination (Lcb,ptype_Lamc, lamb, pi_m, pi0, 0.08);
        setUserInfo(Lc,2);
        setUserInfo(Lcb,2);
        
        combination (Lc,ptype_Lamc, p_p, k_m, pi_p, 0.05);
        combination (Lcb,ptype_Lamc, p_m, k_p, pi_m, 0.05);
        setUserInfo(Lc,5);
        setUserInfo(Lcb,5);
        
        if(!(nevent%100))std::cout<<nevent<<" event. Number of candidates p = " << p_p.size() << "; pbar = " << p_m.size() << "; pi+ = "<< pi_p.size() << "; pi- = "<< pi_m.size() << "; K+ = "<< k_p.size() << "; K- = "<< k_m.size() << "; K_S = "<< k_s.size() << "; pi0 = "<< pi0.size() << "; D0 = " << D0.size() << "; D0bar = " << D0_b.size() << "; D+ = " << D_p.size() << "; D- = "<< D_m.size() << "; Dst0 = " << Dst0.size() << "; Dst0_b = " << Dst0_b.size() << "; D*+ = " << Dst_p.size() << "; D*- = " << Dst_m.size() << "; Lam = " << lam.size() << "; Lam_bar = "<< lamb.size() << "; Lam_c = " << Lc.size() << "; Lam_c_bar = " << Lcb.size() << "; e+ = " << e_p.size() <<"; e- = " << e_m.size() <<"; mu+ = " << mu_p.size() <<"; mu- = " << mu_m.size() << "; Number of recoil candidates L_ = " << L_.size() << "; L_b = " << L_b.size() << '\n';
        
     
        //######################################   FINAL SELECTION
            
        for (std::vector<Particle>::iterator a=L_.begin(); a!=L_.end();++a)
        {
            Particle &ALamC=*a;
            HepLorentzVector momentum=ALamC.p();
            int charge_tag= calcuCharge (&ALamC);
            
            double rmx = (pUPS-momentum).mag(), rm=-1000, Mvis=-1000, rmNopi0 = -1000;//, rm =(pUPS-(momentum+LamC.p())).mag();
            
            if (rmx > 1.5 && rmx < 2.6) 
            {
                //std::cout<<nevent<<" Selected!" << endl;
                int tag = dynamic_cast<UserInfo&>(ALamC.userInfo()).channel(),
                dstch=-1, dch=-1, lcch=0, mconsaddpi0 = 0;
                double mD=-1, mDst=-1, mKs=-1, mL=-1, mLc=-1, hl = -10, hlc = -10, cosW = -10, angchi = -10, 
                pvis=-10, qW=1000, p_2d_from_lamc_cms=-1, p_2d_from_lamc_labs=-1, dphi_lc_lam = -1000, plamlcs = -10, plamclab = -10, plamccms=-10, coslamclab=-10, coslamccms=-10;
        
                if (tag==1 || tag==11 || tag==12)
                {
                    dch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                    mD =  dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                    if( dch==3 || dch==6 ) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                    //mKs = ALamC.child(1).child(0).mass();
                }
                else if (tag==2 || tag==22)
                {
                    dch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                    mD =  dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                    if( dch==2) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                        else if (dch == 3) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(2).userInfo()).mass();
                }
                else if (tag==3 || tag==32)
                {
                    dstch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                    mDst = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                    dch = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).channel();
                    mD = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                    if( dstch==1 && (dch==3 || dch==6)) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(0).userInfo()).mass();
                    else if (dstch==2 && dch==2) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(0).userInfo()).mass();
                    else if (dstch==2 && dch == 3) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(2).userInfo()).mass();
                }
                else if (tag==4 || tag==41 || tag==42)
                {
                    dstch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                    mDst = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                    dch = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).channel();
                    mD = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                    if( dch==3 || dch==6 ) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(0).userInfo()).mass();
                }

                if (Lcb.size()>0) for(std::vector<Particle>::iterator l=Lcb.begin(); l!=Lcb.end();++l)
                {
                    if ( checkSame(*a,*l) ) continue;
                    Particle &LamC=*l;
                    
                    int addpi=0, addpi0=0, totcharge=calcuCharge(&(*l))+calcuCharge(&(*a));
                    for(std::vector<Particle>::iterator i=pi0.begin(); i!=pi0.end();++i)
                    {
                        if( !(checkSame(*l,*i)||checkSame(*a,*i)) ) 
                        {
                            addpi0++;
                            if ( (LamC.p()+i->p()).mag()<2.265) mconsaddpi0++;
                        }
                    }
                    for(std::vector<Particle>::iterator pi = pions.begin(); pi!=pions.end(); ++pi)
                    {
                        if( !(checkSame(*l,*pi)||checkSame(*a,*pi)) ) 
                        {   addpi++;
                            totcharge += pi->charge();
                        }
                    }
                    
                    rm =(pUPS-(momentum+LamC.p())).mag();
                    Mvis = (momentum+LamC.p()).mag();
                    lcch = dynamic_cast<UserInfo&>(LamC.userInfo()).channel(); 
                    mLc = LamC.mass();          
                   
                    pvis = pStar(momentum+LamC.p(),elec,posi).vect().mag();
                    if(lcch==2 || lcch==6 || lcch==7)
                    {
                        rmNopi0 = (pUPS-(momentum+LamC.child(0).p()+LamC.child(1).p())).mag();
                    }
                    
                    // lamc heli
                    HepLorentzVector p_lamc;
                    if (lcch==1 || lcch==2 || lcch==5) p_lamc=LamC.p();
                    else p_lamc = pUPS-momentum;
                    hlc = -cos(heli(LamC.child(0).p(),pUPS,p_lamc)); 
                    
                    
                    //lam heli
                    HepLorentzVector p_proton_from_lam, p_pi_from_lam; 
                    if(lcch!=5)
                    {
                        mL  = dynamic_cast<UserInfo&>(LamC.child(0).userInfo()).mass();
                        if (abs(LamC.child(0).child(0).lund())>1000)
                        {
                            p_proton_from_lam=LamC.child(0).child(0).p(); 
                            p_pi_from_lam=LamC.child(0).child(1).p();
                        }
                        else
                        {
                            p_proton_from_lam=LamC.child(0).child(1).p();
                            p_pi_from_lam=LamC.child(0).child(0).p();
                        }
                        hl  = -cos(heli(p_proton_from_lam, p_lamc,  LamC.child(0).p())); 
                        
                        Hep3Vector norm_lam_c, norm_lam;
                        norm_lam_c = (boostT(LamC.child(0).p(), p_lamc).vect()).cross(boostT(pUPS, p_lamc).vect());
                        norm_lam_c = norm_lam_c*(1./norm_lam_c.mag());
                        norm_lam = (boostT(p_proton_from_lam, p_lamc).vect()).cross(boostT(p_pi_from_lam, p_lamc).vect());
                        norm_lam = norm_lam*(1./norm_lam.mag());
                        dphi_lc_lam=norm_lam_c.angle(norm_lam);
                        
                        if(boostT(p_proton_from_lam, p_lamc).vect().dot(norm_lam_c) < 0) dphi_lc_lam = -dphi_lc_lam;
                    }
            
                    //q = sqrt((P_Lc - P_L)^2) OR sqrt((P_UPS-P_X-P_L)^2)
                    HepLorentzVector p_W_from_lamc;
                    p_W_from_lamc = pUPS-LamC.child(0).p()-momentum;
                    
                    
                    if ((lcch==1) || (lcch==2) || (lcch==5)) qW = (LamC.p()-LamC.child(0).p()).mag(); 
                    else qW =  p_W_from_lamc.mag();	
                    
                    //W heli, angle chi and additional pi0 and pi+pi-
                    if (lcch==3 || lcch==4 || lcch==6 || lcch==7)
                    {
                        HepLorentzVector p_l, p_nu;
                        p_l = LamC.child(1).p();
                        
                        cosW = -cos(heli(p_l,p_lamc,p_W_from_lamc));
            
                        p_nu = p_W_from_lamc-p_l;
            
                        Hep3Vector norm_lambda, norm_W;
                    
                        norm_lambda = (boostT(p_proton_from_lam, p_lamc).vect()).cross(boostT(p_pi_from_lam, p_lamc).vect());
                        norm_lambda = norm_lambda*(1./norm_lambda.mag());
                        norm_W = (boostT(p_nu, p_lamc).vect()).cross(boostT(p_l, p_lamc).vect()); 
                        norm_W = norm_W*(1./norm_W.mag());
            
                        angchi = norm_lambda.angle(norm_W);
                        if(boostT(p_nu, p_lamc).vect().dot(norm_lambda) < 0) angchi = -angchi;
                    }
                    
                    //Lam_c 2nd daughter's momentum
                    HepLorentzVector P_2d_from_LamC=LamC.child(1).p();
                    p_2d_from_lamc_cms=pStar(P_2d_from_LamC,elec,posi).vect().mag(); 
                    p_2d_from_lamc_labs=P_2d_from_LamC.vect().mag();
                    
                    
                    //Lam momentum in Lam_c system
                    plamlcs = boostT(LamC.child(0).p(), p_lamc).vect().mag();
                    
                    //lamc visible momentum
                    plamclab = p_lamc.vect().mag();
                    plamccms = pStar(p_lamc,elec,posi).vect().mag();
                    coslamclab = p_lamc.vect().cosTheta();
                    coslamccms = pStar(p_lamc,elec,posi).vect().cosTheta();
                    
                    t1 -> column("tag", tag);
                    t1 -> column("rmx", rmx);
                    t1 -> column("dstch", dstch);
                    t1 -> column("mdst", mDst);
                    t1 -> column("dch", dch);
                    t1 -> column("md", mD);
                    t1 -> column("px", pStar(momentum,elec,posi).vect().mag());
                    t1 -> column("fox", fox);  
                    t1 -> column("ecms",pUPS.mag());
                    t1 -> column("mks",mKs);
                    t1 -> column("ch_tag", charge_tag);
                
                    t1 -> column("lcch",lcch);
                    t1 -> column("ml", mL);
                    t1 -> column("mlc", mLc);
                    t1->column("hl", hl);
                    t1 -> column("hlc", hlc);
                
                    t1 -> column("q", qW);
                    t1->column("hw",cosW);
                    t1->column("chi",angchi);
                    t1->column("plslc",plamlcs);
                
                    t1 -> column("pvis",pvis);
                    t1 -> column("rmvis", rm);
                    t1 -> column("mvis", Mvis);
                    t1 -> column("rmvis_nopi0", rmNopi0);
                    
                    t1 -> column("philclam",dphi_lc_lam);
                    t1 -> column("lcp2dcm",p_2d_from_lamc_cms);
                    t1 -> column("lcp2dlab",p_2d_from_lamc_labs);
                    
                    t1 -> column("plamclab",plamclab);
                    t1 -> column("plamccms",plamccms);
                    t1 -> column("coslclab",coslamclab);
                    t1 -> column("coslccms",coslamccms);
                    
                    t1 -> column("addpi",addpi);
                    t1 -> column("addpi0",addpi0);
                    t1 -> column("mcwadpi0",mconsaddpi0);
                    t1 -> column("totcharg",totcharge);
                    t1->dumpData();
                }
                else
                {
                    t1 -> column("tag", tag);
                    t1 -> column("rmx", rmx);
                    t1 -> column("dstch", dstch);
                    t1 -> column("mdst", mDst);
                    t1 -> column("dch", dch);
                    t1 -> column("md", mD);
                    t1 -> column("px", pStar(momentum,elec,posi).vect().mag());
                    t1 -> column("fox", fox);  
                    t1 -> column("ecms",pUPS.mag());
                    t1 -> column("mks",mKs);
                    t1 -> column("ch_tag", charge_tag);
                
                    t1 -> column("lcch",lcch);
                    t1 -> column("ml", mL);
                    t1 -> column("mlc", mLc);
                    t1->column("hl", hl);
                    t1 -> column("hlc", hlc);
                
                    t1 -> column("q", qW);
                    t1->column("hw",cosW);
                    t1->column("chi",angchi);
                    t1->column("plslc",plamlcs);
                
                    t1 -> column("pvis",pvis);
                    t1 -> column("rmvis", rm);
                    t1 -> column("mvis", Mvis);
                    t1 -> column("rmvis_nopi0", rmNopi0);
                    
                    t1 -> column("philclam",dphi_lc_lam);
                    t1 -> column("lcp2dcm",p_2d_from_lamc_cms);
                    t1 -> column("lcp2dlab",p_2d_from_lamc_labs);
                    
                    t1 -> column("plamclab",plamclab);
                    t1 -> column("plamccms",plamccms);
                    t1 -> column("coslclab",coslamclab);
                    t1 -> column("coslccms",coslamccms);
                    
                    t1 -> column("addpi",-1);
                    t1 -> column("addpi0",-1);
                    t1 -> column("mcwadpi0",mconsaddpi0);
                    t1 -> column("totcharg",-100);
                    t1->dumpData();
                }
            }
        }
        
        
        
        
        for (std::vector<Particle>::iterator a=L_b.begin(); a!=L_b.end();++a)
        {
            Particle &ALamC=*a;
            HepLorentzVector momentum=ALamC.p();
            int charge_tag= calcuCharge (&ALamC);
            
            double rmx = (pUPS-momentum).mag(), rm=-1000, Mvis=-1000, rmNopi0 = -1000;//, rm =(pUPS-(momentum+LamC.p())).mag();
            
            if (rmx > 1.5 && rmx < 2.6) 
            {
                int tag = dynamic_cast<UserInfo&>(ALamC.userInfo()).channel(),
                dstch=-1, dch=-1, lcch=0, mconsaddpi0=0;
                double mD=-1, mDst=-1, mKs=-1, mL=-1, mLc=-1, hl = -10, hlc = -10, cosW = -10, angchi = -10, 
                       pvis=-10, qW=1000, p_2d_from_lamc_cms=-1, p_2d_from_lamc_labs=-1, dphi_lc_lam = -1000, plamlcs = -10,plamclab = -10, plamccms=-10, coslamclab=-10, coslamccms=-10;
        
                if (tag==1 || tag==11 || tag==12)
                {
                    dch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                    mD =  dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                    if( dch==3 || dch==6 ) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                    //mKs = ALamC.child(1).child(0).mass();
                }
                else if (tag==2 || tag==22)
                {
                    dch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                    mD =  dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                     if( dch==2) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                        else if (dch == 3) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(2).userInfo()).mass();
                }
                else if (tag==3 || tag==32)
                {
                    dstch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                    mDst = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                    dch = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).channel();
                    mD = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                    if( dstch==1 && (dch==3 || dch==6)) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(0).userInfo()).mass();
                    else if (dstch==2 && dch==2) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(0).userInfo()).mass();
                    else if (dstch==2 && dch == 3) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(2).userInfo()).mass();
                }
                else if (tag==4 || tag==41 || tag==42)
                {
                    dstch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                    mDst = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                    dch = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).channel();
                    mD = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                    if( dch==3 || dch==6 ) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(0).userInfo()).mass();
                }
                
                if (Lc.size()>0) for(std::vector<Particle>::iterator l=Lc.begin(); l!=Lc.end();++l)
                {
                    if ( checkSame(*a,*l) ) continue;
                    Particle &LamC=*l;
                    
                    int addpi=0, addpi0=0, totcharge=calcuCharge(&(*l))+calcuCharge(&(*a));
                    for(std::vector<Particle>::iterator i=pi0.begin(); i!=pi0.end();++i)
                    {
                        if( !(checkSame(*l,*i)||checkSame(*a,*i)) ) 
                        {
                            addpi0++;
                            if ( (LamC.p()+i->p()).mag()<2.265) mconsaddpi0++;
                        }
                    }
                    for(std::vector<Particle>::iterator pi = pions.begin(); pi!=pions.end(); ++pi)
                    {
                        if( !(checkSame(*l,*pi)||checkSame(*a,*pi)) ) 
                        {   
                            addpi++;
                            totcharge += pi->charge();
                        }
                    }
                    
                    rm =(pUPS-(momentum+LamC.p())).mag();
                    Mvis = (momentum+LamC.p()).mag();
                    

                    lcch = dynamic_cast<UserInfo&>(LamC.userInfo()).channel(); 
                    mLc = LamC.mass();          
                    pvis = pStar(momentum+LamC.p(),elec,posi).vect().mag();
                    if(lcch==2 || lcch==6 || lcch==7)
                    {
                        rmNopi0 = (pUPS-(momentum+LamC.child(0).p()+LamC.child(1).p())).mag();
                    }
                    
                    // lamc heli
                    HepLorentzVector p_lamc;
                    if (lcch==1 || lcch==2 || lcch==5) p_lamc=LamC.p();
                    else p_lamc = pUPS-momentum;
                    hlc = -cos(heli(LamC.child(0).p(),pUPS,p_lamc)); 
                    
                    
                    //lam heli and dphi_lamc_lam
                    HepLorentzVector p_proton_from_lam, p_pi_from_lam; 
                    if(lcch!=5)
                    {
                        mL  = dynamic_cast<UserInfo&>(LamC.child(0).userInfo()).mass();
                        if (abs(LamC.child(0).child(0).lund())>1000)
                        {
                            p_proton_from_lam=LamC.child(0).child(0).p(); 
                            p_pi_from_lam=LamC.child(0).child(1).p();
                        }
                        else
                        {
                            p_proton_from_lam=LamC.child(0).child(1).p();
                            p_pi_from_lam=LamC.child(0).child(0).p();
                        }
                        hl  = -cos(heli(p_proton_from_lam, p_lamc,  LamC.child(0).p()));
                        
                        Hep3Vector norm_lam_c, norm_lam;
                        norm_lam_c = (boostT(LamC.child(0).p(), p_lamc).vect()).cross(boostT(pUPS, p_lamc).vect());
                        norm_lam_c = norm_lam_c*(1./norm_lam_c.mag());
                        norm_lam = (boostT(p_proton_from_lam, p_lamc).vect()).cross(boostT(p_pi_from_lam, p_lamc).vect());
                        norm_lam = norm_lam*(1./norm_lam.mag());
                        
                        dphi_lc_lam=norm_lam_c.angle(norm_lam);
                        if(boostT(p_proton_from_lam, p_lamc).vect().dot(norm_lam_c) < 0) dphi_lc_lam = -dphi_lc_lam;
                    }
            
                    //q = sqrt((P_Lc - P_L)^2) OR sqrt((P_UPS-P_X-P_L)^2)
                    HepLorentzVector p_W_from_lamc;
                    p_W_from_lamc = pUPS-LamC.child(0).p()-momentum;
                    
                    
                    if ((lcch==1) || (lcch==2) || (lcch==5)) qW = (LamC.p()-LamC.child(0).p()).mag(); 
                    else qW =  p_W_from_lamc.mag();	
                    
                    //W heli and angle chi
                    if (lcch==3 || lcch==4 || lcch==6 || lcch==7)
                    {
                        HepLorentzVector p_l, p_nu;
                        p_l = LamC.child(1).p();
                        
                        cosW = -cos(heli(p_l,p_lamc,p_W_from_lamc));
            
                        p_nu = p_W_from_lamc-p_l;
            
                        Hep3Vector norm_lambda, norm_W;
                    
                        norm_lambda = (boostT(p_proton_from_lam, p_lamc).vect()).cross(boostT(p_pi_from_lam, p_lamc).vect());
                        norm_lambda = norm_lambda*(1./norm_lambda.mag());
                        norm_W = (boostT(p_nu, p_lamc).vect()).cross(boostT(p_l, p_lamc).vect()); 
                        norm_W = norm_W*(1./norm_W.mag());
            
                        angchi = norm_lambda.angle(norm_W);
                        if(boostT(p_nu, p_lamc).vect().dot(norm_lambda) < 0) angchi = -angchi;
                    }
                    
                    //Lam_c 2nd daughter's momentum
                    HepLorentzVector P_2d_from_LamC=LamC.child(1).p();
                    p_2d_from_lamc_cms=pStar(P_2d_from_LamC,elec,posi).vect().mag(); 
                    p_2d_from_lamc_labs=P_2d_from_LamC.vect().mag();
                    
                    //Lam momentum in Lam_c system
                    plamlcs = boostT(LamC.child(0).p(), p_lamc).vect().mag();
                    
                    //lamc visible momentum
                    plamclab = p_lamc.vect().mag();
                    plamccms = pStar(p_lamc,elec,posi).vect().mag();
                    coslamclab = p_lamc.vect().cosTheta();
                    coslamccms = pStar(p_lamc,elec,posi).vect().cosTheta();
                    
                    t1 -> column("tag", tag);
                    t1 -> column("rmx", rmx);
                    t1 -> column("dstch", dstch);
                    t1 -> column("mdst", mDst);
                    t1 -> column("dch", dch);
                    t1 -> column("md", mD);
                    t1 -> column("px", pStar(momentum,elec,posi).vect().mag());
                    t1 -> column("fox", fox);  
                    t1 -> column("ecms",pUPS.mag());
                    t1 -> column("mks",mKs);
                    t1 -> column("ch_tag", charge_tag);
                
                    t1 -> column("lcch",lcch);
                    t1 -> column("ml", mL);
                    t1 -> column("mlc", mLc);
                    t1->column("hl", hl);
                    t1 -> column("hlc", hlc);
                
                    t1 -> column("q", qW);
                    t1->column("hw",cosW);
                    t1->column("chi",angchi);
                    t1->column("plslc",plamlcs);
                
                    t1 -> column("pvis",pvis);
                    t1 -> column("rmvis", rm);
                    t1 -> column("mvis", Mvis);
                    t1 -> column("rmvis_nopi0", rmNopi0);
                    
                    t1 -> column("philclam",dphi_lc_lam);
                    t1 -> column("lcp2dcm",p_2d_from_lamc_cms);
                    t1 -> column("lcp2dlab",p_2d_from_lamc_labs);
                    
                    t1 -> column("plamclab",plamclab);
                    t1 -> column("plamccms",plamccms);
                    t1 -> column("coslclab",coslamclab);
                    t1 -> column("coslccms",coslamccms);
                    
                    t1 -> column("addpi",addpi);
                    t1 -> column("addpi0",addpi0);
                    t1 -> column("mcwadpi0",mconsaddpi0);
                    t1 -> column("totcharg",totcharge);
                    t1->dumpData();
                }
                else
                {
                    t1 -> column("tag", tag);
                    t1 -> column("rmx", rmx);
                    t1 -> column("dstch", dstch);
                    t1 -> column("mdst", mDst);
                    t1 -> column("dch", dch);
                    t1 -> column("md", mD);
                    t1 -> column("px", pStar(momentum,elec,posi).vect().mag());
                    t1 -> column("fox", fox);  
                    t1 -> column("ecms",pUPS.mag());
                    t1 -> column("mks",mKs);
                    t1 -> column("ch_tag", charge_tag);
                
                    t1 -> column("lcch",lcch);
                    t1 -> column("ml", mL);
                    t1 -> column("mlc", mLc);
                    t1->column("hl", hl);
                    t1 -> column("hlc", hlc);
                
                    t1 -> column("q", qW);
                    t1->column("hw",cosW);
                    t1->column("chi",angchi);
                    t1->column("plslc",plamlcs);
                
                    t1 -> column("pvis",pvis);
                    t1 -> column("rmvis", rm);
                    t1 -> column("mvis", Mvis);
                    t1 -> column("rmvis_nopi0", rmNopi0);
                    
                    t1 -> column("philclam",dphi_lc_lam);
                    t1 -> column("lcp2dcm",p_2d_from_lamc_cms);
                    t1 -> column("lcp2dlab",p_2d_from_lamc_labs);
                    
                    t1 -> column("plamclab",plamclab);
                    t1 -> column("plamccms",plamccms);
                    t1 -> column("coslclab",coslamclab);
                    t1 -> column("coslccms",coslamccms);
                    
                    t1 -> column("addpi",-1);
                    t1 -> column("addpi0",-1);
                    t1 -> column("mcwadpi0",mconsaddpi0);
                    t1 -> column("totcharg",-100);
                    t1->dumpData();
                }
                
            }
        }
        
    }
    
    void withdRdZcut(std::vector<Particle> &p,double ip_position, double drcut, double dzcut)
    {
        for(std::vector<Particle>::iterator i=p.begin(); i!=p.end();++i)
            if (abs(i->mdstCharged().trk().mhyp(3).helix(0))>drcut||
                abs(i->mdstCharged().trk().mhyp(3).helix(3)-ip_position)>dzcut)
            {p.erase(i);--i;}
    }
    
    
    void makeLam(std::vector<Particle> &lam, std::vector<Particle> &lamb)
    {
        Mdst_vee2_Manager &veeMgr = Mdst_vee2_Manager::get_manager();
        for(std::vector<Mdst_vee2>::iterator i = veeMgr.begin(); i != veeMgr.end(); ++i)
        {
            // lambda or anti-lambda
            if (!(i->kind()==2 || i->kind()==3))
                continue;
            
            //create particle
            Particle tmp(*i);
            
            //  check proton id
            double prob_kpr=1, prob_pipr=1; 
            if(abs(tmp.child(0).lund())>1000) 
            {
                prob_kpr = atc_pid(3, 1, 5, 3, 4).prob(tmp.child(0).mdstCharged()); 
                prob_pipr = atc_pid(3, 1, 5, 2, 4).prob(tmp.child(0).mdstCharged());
            }
            if(abs(tmp.child(1).lund())>1000) 
            {
                prob_kpr=atc_pid(3, 1, 5, 3, 4).prob(tmp.child(1).mdstCharged());
                prob_pipr = atc_pid(3, 1, 5, 2, 4).prob(tmp.child(1).mdstCharged());
            }
            //if( (prob_kpr > 0.4) || (prob_pipr > 0.4) ) continue;
            if(prob_kpr > 0.9 ) continue;
            
            
            // check mass and flight dist
            HepPoint3D V(tmp.mdstVee2().vx(),tmp.mdstVee2().vy(),0);
            Vector3 P(tmp.px(),tmp.py(),0);
            V=V-IpProfile::position();
            V.setZ(0.);
            if (abs(tmp.mass()-1.115683)>0.015 || V.perp()<0.1 ||
                V.angle(P)>0.01 || tmp.mdstVee2().z_dist()>10. ) 
                continue;
            
            //add to list
            if(i->kind() == 2)
                lam.push_back(tmp);
            if(i->kind() == 3)
                lamb.push_back(tmp);
        }
    }
    
    
    void makeProton (std::vector<Particle> &p, std::vector<Particle> &pbar)
    {
        Mdst_charged_Manager& ChgMgr = Mdst_charged_Manager::get_manager();
        for(std::vector<Mdst_charged>::iterator it = ChgMgr.begin(); it != ChgMgr.end(); ++it)
        {
            if(good_charged(*it) == 0) continue;
            double prob_kpr = atc_pid(3, 1, 5, 3, 4).prob(*it);
            double prob_pipr = atc_pid(3, 1, 5, 2, 4).prob(*it);
            if( (prob_kpr > 0.4) || (prob_pipr > 0.4) ) continue;
            Ptype ptype_PP("P+");
            Ptype ptype_AP("AP+");
            
            if((*it).charge() > 0.)
            {
                Particle tmp1(*it, ptype_PP);
                p.push_back(tmp1);
            }
            else
            {
                Particle tmp1(*it, ptype_AP);
                pbar.push_back(tmp1);
            }
        }
    }
    
    
    
    
    
    void setPi0Error(Particle &p){
        if( !p.child(0) || !p.child(1) ) return;
        if( !p.child(0).mdstGamma() ||
            !p.child(1).mdstGamma()    ) return;
        for(int i=0; i<2; i++){
            double E     = p.child(i).mdstGamma().ecl().energy();
            double phi   = p.child(i).mdstGamma().ecl().phi();
            double theta = p.child(i).mdstGamma().ecl().theta();
            
            double E2  = E*E;
            double E4  = E2*E2;
            double ct2 = cos(theta)*cos(theta);
            double st2 = sin(theta)*sin(theta);
            double st4 = st2*st2;
            
            const HepPoint3D pivot(0.,0.,0.);
            HepMatrix  tmp_a(5,1);
            tmp_a[0][0] = 0.;
            tmp_a[1][0] = phi-M_PI/2;
            tmp_a[2][0] = 1/E/sin(theta);
            tmp_a[3][0] = 0.;
            tmp_a[4][0] = tan(M_PI/2-theta);
            HepVector  a(tmp_a);
            
            double errE     = p.child(i).mdstGamma().ecl().error(0);
            double errPHI   = p.child(i).mdstGamma().ecl().error(2);
            double errTHETA = p.child(i).mdstGamma().ecl().error(5);
            
            HepSymMatrix Ea(5,0);
            Ea[0][0] = 1.;
            Ea[1][1] = errPHI;
            Ea[2][2] = errE/E4/st2 + errTHETA*ct2/E2/st4;
            Ea[3][3] = 1.;
            Ea[4][2] = errTHETA*cos(theta)/E/st4;
            Ea[4][4] = errTHETA/st4;
            
            Helix helix(pivot, a, Ea);
            
            HepLorentzVector momentum;
            HepPoint3D position;
            HepSymMatrix error(7,0);
            
            momentum = helix.momentum(0.,0.,position,error);
            p.child(i).momentum().momentumPosition(momentum,position,error);
        }
        doMassVertexFit(p);
    }
    
    void setGammaError(Particle &g)
    {
        if( !g.mdstGamma() ) return;
        double E     = g.mdstGamma().ecl().energy();
        double phi   = g.mdstGamma().ecl().phi();
        double theta = g.mdstGamma().ecl().theta();
        
        double E2  = E*E;
        double E4  = E2*E2;
        double ct2 = cos(theta)*cos(theta);
        double st2 = sin(theta)*sin(theta);
        double st4 = st2*st2;
        
        const HepPoint3D pivot(0.,0.,0.);
        HepMatrix  tmp_a(5,1);
        tmp_a[0][0] = 0.;
        tmp_a[1][0] = phi-M_PI/2;
        tmp_a[2][0] = 1/E/sin(theta);
        tmp_a[3][0] = 0.;
        tmp_a[4][0] = tan(M_PI/2-theta);
        HepVector  a(tmp_a);
        
        double errE     = g.mdstGamma().ecl().error(0);
        double errPHI   = g.mdstGamma().ecl().error(2);
        double errTHETA = g.mdstGamma().ecl().error(5);
        
        HepSymMatrix Ea(5,0);
        Ea[0][0] = 1.;
        Ea[1][1] = errPHI;
        Ea[2][2] = errE/E4/st2 + errTHETA*ct2/E2/st4;
        Ea[3][3] = 1.;
        Ea[4][2] = errTHETA*cos(theta)/E/st4;
        Ea[4][4] = errTHETA/st4;
        
        Helix helix(pivot, a, Ea);
        
        HepLorentzVector momentum;
        HepPoint3D position;
        HepSymMatrix error(7,0);
        
        momentum = helix.momentum(0.,0.,position,error);
        g.momentum().momentumPosition(momentum,position,error);
    }
    
    
    void doMassVertexFit(class vector<Particle> &p_list, double mass)
    {
        for(vector<Particle>::iterator i=p_list.begin(); i!=p_list.end();++i)
            doMassVertexFit(*i, mass);
    }
    void doVertexFit(class vector<Particle> &p_list)
    {
        for(vector<Particle>::iterator i=p_list.begin(); i!=p_list.end();++i)
            doVertexFit(*i);
    }
    
    void setPi0Error(vector<Particle> &p_list)
    {
        for(vector<Particle>::iterator i=p_list.begin(); i!=p_list.end(); ++i)
            setPi0Error(*i);
    }
    
    void setGammaError(vector<Particle> &p_list)
    {
        for(vector<Particle>::iterator i=p_list.begin(); i!=p_list.end(); ++i)
            setGammaError(*i);
    }
    
    
    void doMassVertexFit(class Particle & P, double mass)
    {
        
        if(! & P.userInfo() ) setUserInfo(P);
        unsigned err;
        dynamic_cast<UserInfo&>(P.userInfo()).mass(P.mass());
        kmassvertexfitter kmvf;
        if (mass==-1)
            kmvf.invariantMass(P.pType().mass());
        else
            kmvf.invariantMass(mass);
        for (int j=0; j<P.nChildren();++j)
            addTrack2fit(kmvf,P.child(j));
        err = kmvf.fit();  //do "fitting".                                                                                                                                                                                                                                                                                                                                                       
        if(err == 0)
        {
            makeMother(kmvf,P);
            dynamic_cast<UserInfo&>(P.userInfo()).vchisq(kmvf.chisq()/kmvf.dgf());
        }
        else dynamic_cast<UserInfo&>(P.userInfo()).vchisq(-10);
    }
    
    
    
    void doVertexFit(class Particle & P)
    {
        if(! & P.userInfo() ) setUserInfo(P);
        unsigned err;
        dynamic_cast<UserInfo&>(P.userInfo()).mass(P.mass());
        kvertexfitter kvf;
        for (int j=0; j<P.nChildren();++j)
            addTrack2fit(kvf,P.child(j));
        err = kvf.fit();  //do "fitting".                                                                                                                                                                                                                                                                                                                                                        
        if(err == 0)
        {
            makeMother(kvf,P);
            dynamic_cast<UserInfo&>(P.userInfo()).vchisq(kvf.chisq()/kvf.dgf());
        }
        else dynamic_cast<UserInfo&>(P.userInfo()).vchisq(-10);
    }
    
    double heli (  HepLorentzVector  * P4A, HepLorentzVector * P4B, HepLorentzVector* P4system) // return angle btw P4A and P4B in P4system rest frame
    {
        return (boostT(P4A, P4system).vect()).angle(boostT(P4B, P4system).vect());
    }
    
    double heli (  HepLorentzVector   P4A, HepLorentzVector  P4B, HepLorentzVector P4system) // return angle btw P4A and P4B in P4system rest frame
    {
        return (boostT(P4A, P4system).vect()).angle(boostT(P4B, P4system).vect());
    }
    
    
    HepLorentzVector boostT( HepLorentzVector p, HepLorentzVector p_boost) //p --boost--> p_boost;  
    {                                                                            
        HepLorentzVector tmp(p);
        tmp.boost(-p_boost.boostVector());
        
        return tmp;
    }
    
    
    HepLorentzVector boostT( HepLorentzVector *p, HepLorentzVector *p_boost)
    {                                                                
        HepLorentzVector tmp(*p);
        tmp.boost(-p_boost->boostVector());
        return tmp;
    }
    int calcuCharge (Particle *p )
    {
        if (p->nChildren()==0)
            return p->charge();
        int ch=0;
        for (int k=0 ; k!=p->nChildren(); k++)
            ch+=calcuCharge(&(p->child(k)));
        return ch;
    }
    
    void eraseLeptons(std::vector<Particle> &list)
    {
        for(int i=0;i<(int)list.size();++i)
        {
            if(list[i].mdstCharged())
            {
                eid post_elid(list[i].mdstCharged());
                Muid_mdst muid(list[i].mdstCharged());    	 
                if((post_elid.prob() > 0.9)  || (muid.Muon_likelihood() > 0.95)) 
                {
                    list.erase(list.begin()+i);
                    --i;
                }
            }
            else
            {
                list.erase(list.begin()+i);
                --i;
            }
        }
    }
    
}
