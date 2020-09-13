#include "User_reco.h"
#include "userinfo.h"
#include EVTCLS_H //R2 distribution

//e+ e- -> x+Lc

using namespace Belle;
namespace Belle {
    void User_reco::hist_def( void )
    {
        
        extern BelleTupleManager* BASF_Histogram;
        t1 = BASF_Histogram->ntuple ("data","tag dch dstch md mdst rmx px fox ecms mks mpi0_d mpi0_dst ch_tag" ); // not ALL momenta in CMS! 	lepton cosTheta in CMS, rholam, rholamcms	
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

        eraseLeptons(p_p);
        eraseLeptons(p_m);
        
        //kaons and pions
        std::vector<Particle>  k_p, k_m, pi_p, pi_m, pions;
        makeKPi(k_p, k_m, pi_p, pi_m,1);
        
        ntrk=k_p.size()+k_m.size();
        
      /*
        withdRdZcut(k_p,runIp.z());
        withdRdZcut(pi_p,runIp.z());
        withdRdZcut(k_m,runIp.z());
        withdRdZcut(pi_m,runIp.z());
       */  
        
        withKaonId(k_p,0.6,3,1,5,3,4);
        withKaonId(k_p,0.6,3,1,5,3,2);
        withKaonId(k_m,0.6,3,1,5,3,4);
        withKaonId(k_m,0.6,3,1,5,3,2);
        
        withPionId(pi_p,0.4,3,1,5,3,2);
        withPionId(pi_p,0.4,3,1,5,4,2);
        withPionId(pi_m,0.4,3,1,5,3,2);
        withPionId(pi_m,0.4,3,1,5,4,2);

        eraseLeptons(k_p);
        eraseLeptons(k_m);
        eraseLeptons(pi_p);
        eraseLeptons(pi_m);
        
        
        for(std::vector<Particle>::iterator l = pi_m.begin(); l!=pi_m.end(); ++l)
            pions.push_back(*l);
        for(std::vector<Particle>::iterator l = pi_p.begin(); l!=pi_p.end(); ++l)
            pions.push_back(*l);
        
        
        
        
        
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
        while(k_s.size()>1)
        {
            for(std::vector<Particle>::iterator l = k_s.begin(); l!=k_s.end(); ++l)
            {
                k_s_chisq = dynamic_cast<UserInfo&>(l->userInfo()).vchisq();
                if((k_s_chisq<0) || (k_s_chisq > bufchisq))
                {
                    k_s.erase(l); 
                    --l;
                    continue;
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
            //leptons    
            /* std::vector<Particle>  e_p,e_m,mu_p,mu_m;
             *       makeLepton(e_p,e_m,mu_p,mu_m);
             *       
             *       withMuId(mu_p);
             *       withMuId(mu_m);
             *       withEId(e_p);
             *       withEId(e_m);
             */

            
            
            
        //photons
        std::vector<Particle> photons;
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
        }
            
            
            
            
            std::vector<Particle> D0, D0_b, D_p, D_m, Dst_p, Dst_m, Dst0, Dst0_b, pi_pm;
            
            //######################################    TAG SIDE
            combination (pi_pm, ptype_D0B, pi_p, pi_m);
            
            combination (D0_b,ptype_D0B, k_p, pi_m, 0.06);
            combination (D0,ptype_D0, k_m, pi_p, 0.06);
            setUserInfo(D0_b, 1); 
            setUserInfo(D0, 1); 
            combination (D0_b,ptype_D0B, k_p, pi_pm, pi_m, 0.05);
            combination (D0,ptype_D0, k_m, pi_pm, pi_p, 0.05);
            setUserInfo(D0_b, 2); 
            setUserInfo(D0, 2); 
            combination (D0_b,ptype_D0B, k_s, pi_pm, 0.05);
            combination (D0,ptype_D0, k_s, pi_pm, 0.05);
            setUserInfo(D0_b, 3);
            setUserInfo(D0, 3); 
            combination (D0_b,ptype_D0B, k_p, pi0, pi_m,  0.06);
            combination (D0,ptype_D0, k_m, pi0, pi_p,  0.06);
            setUserInfo(D0_b, 4);
            setUserInfo(D0, 4); 
            
            combination (D0_b,ptype_D0B, k_p, pi0, pi_m, pi_pm, 0.06);
            combination (D0,ptype_D0, k_m, pi0, pi_p, pi_pm, 0.06);
            setUserInfo(D0_b, 5);
            setUserInfo(D0, 5);
            combination (D0_b,ptype_D0B, k_s, pi0, pi_pm, 0.06);
            combination (D0,ptype_D0, k_s, pi0, pi_pm, 0.06);
            setUserInfo(D0_b, 6);
            setUserInfo(D0, 6);
            
            combination (D_m,ptype_Dm, k_p, pi_m, pi_m, 0.05);
            combination (D_p,ptype_Dp, k_m, pi_p, pi_p, 0.05);
            setUserInfo(D_m, 1);
            setUserInfo(D_p, 1);
            combination (D_m,ptype_Dm, k_s, pi_m, 0.05);
            combination (D_p,ptype_Dp, k_s, pi_p, 0.05);
            setUserInfo(D_m, 2);
            setUserInfo(D_p, 2);
            combination (D_m,ptype_Dm, k_s, pi_pm, pi_m, 0.05);
            combination (D_p,ptype_Dp, k_s, pi_pm, pi_p, 0.05);
            setUserInfo(D_m, 3);
            setUserInfo(D_p, 3);
            combination (D_m,ptype_Dm, k_p, pi_m, k_m, 0.05);
            combination (D_p,ptype_Dp, k_m, pi_p, k_p, 0.05);
            setUserInfo(D_m, 4);
            setUserInfo(D_p, 4);
            
            //doMassVertexFit(D0_b);
            //doMassVertexFit(D_m);
            
            if(D0_b.size()+D_m.size()+D0.size()+D_p.size()==0) return;
            
            combination (Dst0_b,ptype_DstB, D0_b, pi0, 0.2);
            setUserInfo(Dst0_b, 1);
            combination (Dst0_b,ptype_DstB, D0_b, photons, 0.2);
            setUserInfo(Dst0_b, 2);
            
            
            combination (Dst_m,ptype_Dstm, D0_b, pi_m, 0.2);
            setUserInfo(Dst_m, 1);
            combination (Dst_m,ptype_Dstm, D_m, pi0, 0.2);
            setUserInfo(Dst_m, 2);
            
            /*for (std::vector<Particle>::iterator i=Dst0_b.begin(); i!=Dst0_b.end();++i)
            {
                Particle dst0b(*i);
                if(abs(dst0b.mass()-dst0b.child(0).mass()-0.1425)>0.025)
                {
                    Dst0_b.erase(i); 
                    --i;
                }
            }
            */
            
            doMassVertexFit(Dst0_b);
            
            /*for (std::vector<Particle>::iterator i=Dst_m.begin(); i!=Dst_m.end();++i)
            {
                Particle dstm(*i);
                if(abs(dstm.mass()-dstm.child(0).mass()-0.143)>0.025)
                {
                    Dst_m.erase(i); 
                    --i;
                }
            }
            */
            
            doMassVertexFit(Dst_m);
            
            combination (Dst0,ptype_Dst0, D0, pi0, 0.2);
            setUserInfo(Dst0, 1);
            combination (Dst0,ptype_Dst0, D0, photons, 0.2);
            setUserInfo(Dst0, 2);
            
            
            combination (Dst_p,ptype_Dstp, D0, pi_p, 0.2);
            setUserInfo(Dst_p, 1);
            combination (Dst_p,ptype_Dstp, D_p, pi0, 0.2);
            setUserInfo(Dst_p, 2);
            
            /*for (std::vector<Particle>::iterator i=Dst0.begin(); i!=Dst0.end();++i)
            {
                Particle dst0(*i);
                if(abs(dst0.mass()-dst0.child(0).mass()-0.1425)>0.025)
                {
                    Dst0.erase(i); 
                    --i;
                }
            }
            */
            doMassVertexFit(Dst0);
            
            /*for (std::vector<Particle>::iterator i=Dst_p.begin(); i!=Dst_p.end();++i)
            {
                Particle dstp(*i);
                if(abs(dstp.mass()-dstp.child(0).mass()-0.143)>0.025)
                {
                    Dst_p.erase(i); 
                    --i;
                }
            }
            */
            doMassVertexFit(Dst_p);
            
            std::vector <Particle> L_, L_b;
            combination (L_b,ptype_Lamc, p_m, D0_b);
            setUserInfo(L_b, 1);
            combination (L_b,ptype_Lamc, p_m, D_m, pi_p);
            setUserInfo(L_b, 2);
            combination (L_b,ptype_Lamc, p_m, Dst_m, pi_p);
            setUserInfo(L_b, 3);
            combination (L_b,ptype_Lamc, p_m, Dst0_b);
            setUserInfo(L_b, 4);
            
            combination (L_,ptype_ALamc, p_p, D0);
            setUserInfo(L_, 1);
            combination (L_,ptype_ALamc, p_p, D_p, pi_m);
            setUserInfo(L_, 2);
            combination (L_,ptype_ALamc, p_p, Dst_p, pi_m);
            setUserInfo(L_, 3);
            combination (L_,ptype_ALamc, p_p, Dst0);
            setUserInfo(L_, 4);
            
            std::cout<<nevent<<" event. Number of candidates p = " << p_p.size() << "; pbar = " << p_m.size() << "; pi+ = "<< pi_p.size() << "; pi- = "<< pi_m.size() << "; K+ = "<< k_p.size() << "; K- = "<< k_m.size() << "; K_S = "<< k_s.size() << "; pi0 = "<< pi0.size() << "; D0 = " << D0.size() << "; D0bar = " << D0_b.size() << "; D+ = " << D_p.size() << "; D- = "<< D_m.size() << "; gamma = " << photons.size() << "; Dst0 = " << Dst0.size() << "; Dst0_b = " << Dst0_b.size() << "; D*+ = " << Dst_p.size() << "; D*- = " << Dst_m.size() << "; Number of recoil candidates L_ = " << L_.size() << "; L_b = " << L_b.size() << '\n';
            

            
            
            /*doVertexFit(L_b);
             *       double recoil_chisq;
             *       bufchisq=1000000;
             *       while(L_b.size()>1)
             *       {
             *           for(std::vector<Particle>::iterator l = L_b.begin(); l!=L_b.end(); ++l)
             *           {
             *               recoil_chisq = dynamic_cast<UserInfo&>(l->userInfo()).vchisq();
             *               if((recoil_chisq < 0) || (recoil_chisq > bufchisq))
             *               {
             *                   L_b.erase(l); 
             *                   --l;
             *                   continue;
    }
    bufchisq = recoil_chisq;
    }
    }
    
    if (L_b.size()+L_.size()==0) return; 
    if (!(nevent%1000))std::cout<<nevent<<"Best recoil candidate's chisq/ndf = " << recoil_chisq << '\n';*/
            
            
            /*   //######################################   SIGNAL SIDE
             *       
             *       std::vector <Particle> Lc, Lcb; 
             *       
             *       for(std::vector<Particle>::iterator l = lam.begin(); l!=lam.end(); ++l) Lc.push_back(*l);
             *       for(std::vector<Particle>::iterator l = lamb.begin(); l!=lamb.end(); ++l) Lcb.push_back(*l);
             *       setUserInfo(Lc,0);
             *       setUserInfo(Lcb,0);
             *       
             *       combination (Lc,ptype_Lamc, lam, e_p);
             *       combination (Lcb,ptype_Lamc, lamb, e_m);
             *       setUserInfo(Lc,3);
             *       setUserInfo(Lcb,3);  
             *       combination (Lc,ptype_Lamc, lam, e_m);     //fake
             *       combination (Lcb,ptype_Lamc, lamb, e_p);   //fake
             *       setUserInfo(Lc,300);
             *       setUserInfo(Lcb,300);
             *       
             *       combination (Lc,ptype_Lamc, lam, mu_p);
             *       combination (Lcb,ptype_Lamc, lamb, mu_m);
             *       setUserInfo(Lc,4);
             *       setUserInfo(Lcb,-4);
             *       combination (Lc,ptype_Lamc, lam, mu_m);     //fake
             *       combination (Lcb,ptype_Lamc, lamb, mu_p);   //fake
             *       setUserInfo(Lc,400);
             *       setUserInfo(Lcb,400);
             *       
             *       
             *       for(std::vector<Particle>::iterator l = Lc.begin(); l!=Lc.end(); ++l)
             *           if (l->mass()>ptype_Lamc.mass()+0.1)
             *           {Lc.erase(l); --l;}
             *       for(std::vector<Particle>::iterator l = Lcb.begin(); l!=Lcb.end(); ++l)
             *           if (l->mass()>ptype_Lamc.mass()+0.1)
             *           {Lcb.erase(l); --l;}
             *               
             *               
             *       combination (Lc,ptype_Lamc, lam, pi_p,0.1);
             *       combination (Lcb,ptype_Lamc, lamb, pi_m,0.1);
             *       setUserInfo(Lc,1);
             *       setUserInfo(Lcb,1);
             *       
             *       combination (Lc,ptype_Lamc, lam, pi_p, pi0, 0.1);
             *       combination (Lcb,ptype_Lamc, lamb, pi_m, pi0, 0.1);
             *       setUserInfo(Lc,2);
             *       setUserInfo(Lcb,2);
             *       
             *       
             *       
             *        //######################################   FINAL SELECTION
             *       
             *       std::vector<Particle> A; 
             *       combination(A,ptype_UPS4,Lc,L_b);
             *       combination(A,ptype_UPS4,Lcb,L_);
             */
            
            for (std::vector<Particle>::iterator a=L_.begin(); a!=L_.end();++a)
            {
                Particle &ALamC=*a;
                
                HepLorentzVector momentum=ALamC.p();
                int charge_tag= calcuCharge (&ALamC);
                /*   int charge= calcuCharge (&All);
                 *           
                 *           int n_pi = 0, n_pi0 = 0, n_k = 0, n_e = 0, n_mu = 0, n_p = 0, n_gam = 0 ;
                 *           for (std::vector<Particle>::iterator pi=pions.begin(); pi!=pions.end();++pi)
                 *           {
                 *               //pion cuts
                 *               if ( checkSame(*a,*pi) ) continue;
                 *               //
                 *               n_pi++;
                 *               charge+=pi->charge();
                 *               momentum+=pi->p();
            }
            
            for (std::vector<Particle>::iterator k=kaons.begin(); k!=kaons.end();++k)
            {
            //pion cuts
            if ( checkSame(*a,*k) ) continue;
            //
            n_k++;
            charge+=k->charge();
            momentum+=k->p();
            }
            
            HepLorentzVector momentum0 = momentum;
            
            for (std::vector<Particle>::iterator pi=pi0.begin(); pi!=pi0.end();++pi)
            {
            //pion cuts
            if ( checkSame(*a,*pi) ) continue;
            //
            n_pi0++;
            momentum0+=pi->p();
            }
            
            for (std::vector<Particle>::iterator g=photons.begin(); g!=photons.end();++g)
            {
            //pion cuts
            if ( checkSame(*a,*g) ) continue;
            //
            n_gam++;
            momentum0+=g->p();
            }
            
            // FINAL SELECTION
            if (charge!=0) continue;
            */
                double rmx = (pUPS-momentum).mag();//, rm =(pUPS-(momentum+LamC.p())).mag();
                
                if (abs(rmx-2.286)<1.29) 
                {
                    std::cout<<nevent<<" Selected!" << endl;
                    int //lcch = dynamic_cast<UserInfo&>(LamC.userInfo()).channel(), 
                    tag = dynamic_cast<UserInfo&>(ALamC.userInfo()).channel(),
                    dstch=-1, dch;
                    double mD=-1, mDst=-1, mPi0_D=-1, mPi0_Dst=-1, mKs=-1;
                    
                   
                    
                    if (tag<3)
                    {
                        std::cout<<nevent<<" 1" << endl;
                        dch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                        // mD = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                        std::cout<<nevent<<" 2" << endl;
                        mD = ALamC.child(1).mass();
                        //mKs = ALamC.child(1).child(0).mass();
                        std::cout<<nevent<<" 3" << endl;
                        if((dch>3) && (dch<7) && (tag==1))  mPi0_D = dynamic_cast<UserInfo&>(ALamC.child(1).child(1).userInfo()).mass();
                        std::cout<<nevent<<" 4" << endl;
                        if( (((dch==3) || (dch==6)) && (tag==1)) ||  (((dch==2) || (dch==3)) && (tag==2)) ) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                    }
                    else
                    {
                        std::cout<<nevent<<" 5" << endl;
                        dstch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                        std::cout<<nevent<<" 6" << endl;
                        mDst = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                        std::cout<<nevent<<" 7" << endl;
                        dch = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).channel();
                        //mD = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                        std::cout<<nevent<<" 8" << endl;
                        mD = ALamC.child(1).child(0).mass();
                        
                        std::cout<<nevent<<" 9" << endl;
                        if (dstch == 1) mPi0_Dst = ALamC.child(1).child(1).mass();
                        std::cout<<nevent<<" 10" << endl;
                        if((dch>3) && (dch<7) && (tag==4))  mPi0_D = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(1).userInfo()).mass();
                        std::cout<<nevent<<" 11" << endl;
                        if( (((dch==3) || (dch==6)) && (tag==4) && (dstch==1)) ||  (((dch==2) || (dch==3)) && (tag==3) && (dstch==2))) 
                            mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(0).userInfo()).mass();
                        
                        //mKs = ALamC.child(1).child(0).child(0).mass();   
                    }
                    
                    
                    //
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
                    t1 -> column("mpi0_d",mPi0_D);
                    t1 -> column("mpi0_dst",mPi0_Dst);
                    t1 -> column("ch_tag", charge_tag);
                    std::cout<<nevent<<" 12" << endl;
                    /*
                     *               t1 -> column("lcch", lcch);
                     *               t1 -> column("rmvis", rm);
                     *               t1 -> column("npi", n_pi);
                     *               t1 -> column("nk", n_k);
                     *               t1 -> column("npi0", n_pi0);  
                     *               t1 -> column("ngam", n_gam);  
                     *               t1 -> column("pvis",pStar(momentum+LamC.p(),elec,posi).vect().mag());
                     *               if (lcch==0) t1 -> column("ml", LamC.mass());
                     *               
                     *               if (lcch!=0)
                     *               {
                     *                   t1 -> column("ml", dynamic_cast<UserInfo&>(LamC.child(0).userInfo()).mass());
                     *                   t1 -> column("mlc", LamC.mass());
                     *                   
                     *                   // lamc heli
                     *                   HepLorentzVector p_lamc;
                     *                   if (lcch==1) 
                     *                       p_lamc=LamC.p();
                     *                   else 
                     *                       p_lamc = pUPS-momentum;
                     *               
                     *                   t1 -> column("hlc",-cos(heli(LamC.child(0).p(),pUPS,p_lamc)));
                     *               
                     *                   //lam heli
                     *                   HepLorentzVector p_proton_from_lam, p_pi_from_lam; 
                     *                   if (abs(LamC.child(0).child(0).lund())>1000)
                     *                   {
                     *                       p_proton_from_lam=LamC.child(0).child(0).p(); 
                     *                       p_pi_from_lam=LamC.child(0).child(1).p();
                }
                else
                {
                p_proton_from_lam=LamC.child(0).child(1).p();
                p_pi_from_lam=LamC.child(0).child(0).p();
                }
                t1->column("hl",-cos(heli (p_proton_from_lam, p_lamc,  LamC.child(0).p())));
                
                
                //q = sqrt((P_Lc - P_L)^2) OR sqrt((P_UPS-P_X-P_L)^2)
                HepLorentzVector p_W_from_lamc;
                p_W_from_lamc = pUPS-LamC.child(0).p()-momentum;
                
                if ((lcch==1) || (lcch==2))
                    t1 -> column("q",(LamC.p()-LamC.child(0).p()).mag());
                else
                    t1 -> column("q",(p_W_from_lamc).mag()); 	
                
                
                //W heli and angle chi
                HepLorentzVector p_l, p_nu;
                p_l = LamC.child(1).p();
                
                if ((lcch==3) || (lcch==4))
                {
                t1->column("hw",-cos( heli(p_l,p_lamc,p_W_from_lamc)));
                }
                else 
                {
                t1->column("hw",-999);
                }
                
                p_nu = p_W_from_lamc-p_l;
                
                Hep3Vector norm_lambda, norm_W;
                norm_lambda = (boostT(p_proton_from_lam, p_lamc).vect()).cross(boostT(p_pi_from_lam, p_lamc).vect());
                norm_lambda = norm_lambda*(1./norm_lambda.mag());
                norm_W = (boostT(p_nu, p_lamc).vect()).cross(boostT(p_l, p_lamc).vect()); 
                norm_W = norm_W*(1./norm_W.mag());
                
                if ((lcch==3) || (lcch==4))
                {
                t1->column("chi",norm_lambda.angle(norm_W));
                }
                else 
                {
                t1->column("chi",-999);
                }
                }
                */
                    t1->dumpData();
                }
            }
                std::cout<<nevent<<" 13" << endl;
                for (std::vector<Particle>::iterator a=L_b.begin(); a!=L_b.end();++a)
            {
                Particle &ALamC=*a;
                
                HepLorentzVector momentum=ALamC.p();
                int charge_tag= calcuCharge (&ALamC);
                
                double rmx = (pUPS-momentum).mag();//, rm =(pUPS-(momentum+LamC.p())).mag();
                
                if (abs(rmx-2.286)<1.29) 
                {
                    std::cout<<nevent<<" Selected!" << endl;
                    int //lcch = dynamic_cast<UserInfo&>(LamC.userInfo()).channel(), 
                    tag = dynamic_cast<UserInfo&>(ALamC.userInfo()).channel(),
                    dstch=-1, dch;
                    double mD=-1, mDst=-1, mPi0_D=-1, mPi0_Dst=-1, mKs=-1;
                    
                    
                    
                 
                    
                    if (tag<3)
                    {
                        std::cout<<nevent<<" 14" << endl;
                        dch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                        // mD = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                       std::cout<<nevent<<" 15" << endl;
                        mD = ALamC.child(1).mass();
                        //mKs = ALamC.child(1).child(0).mass();
                        std::cout<<nevent<<" 16" << endl;
                        if((dch>3) && (dch<7) && (tag==1))  mPi0_D = dynamic_cast<UserInfo&>(ALamC.child(1).child(1).userInfo()).mass();
                        std::cout<<nevent<<" 17" << endl;
                        if( (((dch==3) || (dch==6)) && (tag==1)) ||  (((dch==2) || (dch==3)) && (tag==2)) ) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                    }
                    else
                    {
                        std::cout<<nevent<<" 18" << endl;
                        dstch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                        std::cout<<nevent<<" 19" << endl;
                        mDst = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                        std::cout<<nevent<<" 20" << endl;
                        dch = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).channel();
                        std::cout<<nevent<<" 21" << endl;
                        //mD = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                        mD = ALamC.child(1).child(0).mass();
                        
                        std::cout<<nevent<<" 22" << endl;
                        if (dstch == 1) mPi0_Dst = ALamC.child(1).child(1).mass();
                        std::cout<<nevent<<" 23" << endl;
                        if((dch>3) && (dch<7) && (tag==4))  mPi0_D = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(1).userInfo()).mass();
                        std::cout<<nevent<<" 24" << endl;
                        if( (((dch==3) || (dch==6)) && (tag==4) && (dstch==1)) ||  (((dch==2) || (dch==3)) && (tag==3) && (dstch==2))) 
                            mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(0).userInfo()).mass();
                        
                        //mKs = ALamC.child(1).child(0).child(0).mass();   
                    }
                    
                    std::cout<<nevent<<" 25" << endl;
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
                    t1 -> column("mpi0_d",mPi0_D);
                    t1 -> column("mpi0_dst",mPi0_Dst);
                    t1 -> column("ch_tag", charge_tag);
                      
                    t1->dumpData();
                }
                std::cout<<nevent<<" 26" << endl;
                
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
            double prob_kpr=1; 
            if(abs(tmp.child(0).lund())>1000) 
                prob_kpr=atc_pid(3, 1, 5, 3, 4).prob(tmp.child(0).mdstCharged()); 
            if(abs(tmp.child(1).lund())>1000) 
                prob_kpr=atc_pid(3, 1, 5, 3, 4).prob(tmp.child(1).mdstCharged()); 
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
