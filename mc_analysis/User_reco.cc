#include "User_reco.h"
#include "userinfo.h"
#include "benergy/BeamEnergy.h"
#include EVTCLS_H //R2 distribution

//e+ e- -> x+Lc

using namespace Belle;
namespace Belle {
    void User_reco::hist_def( void )
    {
        
        extern BelleTupleManager* BASF_Histogram;
        t1 = BASF_Histogram->ntuple ("data_mc","tag dch dstch md mdst rmx rmvis rmvis_nopi0 mvis px plamclab plamccms pvis fox addpi addpi0 totcharg ecms mks ch_tag lcch ml mlc hlc hl q hw chi lcp2dcm lcp2dlab philclam plslc mc_lcch mc_pnu mc_plamc mc_angnv mcanglcx mc_rmx" ); // not ALL momenta in CMS! 	lepton cosTheta in CMS, rholam, rholamcms	
        t2 = BASF_Histogram->ntuple ("gen_mc","fox ecms mc_ecms mlc ch_lamc lcch hlc hl q hw chi lcp2dcm lcp2dlab philclam plslc plamclab plamccms" ); // not ALL momenta in CMS!
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
        static int nevent=0, nlamc_gen=0, nlamc_rec=0;
        nevent++;
        
        
        *status=0;
        
        double elec=8.0, posi=3.5;
        Belle_runhead_Manager& rhdmgr = Belle_runhead_Manager::get_manager();
        Belle_runhead_Manager::const_iterator rhd = rhdmgr.begin();

        if(rhd!=rhdmgr.end() && rhd->EHER() > 5.){
            elec = rhd->EHER();
            posi = rhd->ELER();
        }
        
        
        
        /*posi=BeamEnergy::E_LER();
        elec=BeamEnergy::E_HER();*/
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
        //std::cout << "INITIALIZING!" << std::endl;
        
        //------------------------MONTE CARLO DATA ANALYSIS----------------------------
        
        int lam_c_gen = 0, lam_c_rec=0, idlamc=-1, id=1, idHEP, nlamc_daughters=-1, lam_daid1=-1, lam_daid2=-1, mc_lcch = -1;
        
        Particle mc_LamC, mc_lam, mc_pi, mc_pi0, mc_mu, mc_numu, mc_e, mc_nue, mc_pfromlam, mc_pifromlam;
        
        HepLorentzVector mc_pUPS = HepLorentzVector(0., 0.,0.,0.);
        
        bool lamhere = 0, pihere = 0, pi0here = 0, numuhere = 0, muhere = 0, nuehere = 0, ehere = 0, lamBranch=0, pinlam=0, piinlam=0, lamchere=0;
        
        
        
        
        Gen_hepevt_Manager &evt_manager = Gen_hepevt::get_manager();   
        for (std::vector<Gen_hepevt>::iterator evt = evt_manager.begin(); evt != evt_manager.end(); ++evt) 
        {
            if (!(evt->mother() && evt->mother()==1 && (evt->mother().idhep()==10022 || evt->mother().idhep()==553 || evt->mother().idhep()==100553 || evt->mother().idhep()==200553 || evt->mother().idhep()==300553 || evt->mother().idhep() == 9000553))) continue; 
            if (abs(evt->idhep())>=22)  mc_pUPS+= HepLorentzVector(evt->PX(),evt->PY(),evt->PZ(),evt->E());
        }
  
        //std::cout<<mc_pUPS.px()<<" "<<mc_pUPS.e()<<" "<<mc_pUPS.mag()<<std::endl;
        
        
        
        for (std::vector<Gen_hepevt>::iterator evt = evt_manager.begin(); evt != evt_manager.end(); ++evt) 
        {              
            idHEP = evt->idhep();
            
            if (abs(idHEP)==4122 && lam_c_gen==0) // Lamc=4122   anti-Lamc=-4122  
            {
                /// lam_c is found! 
                lam_c_gen++;
                idlamc = id;
                mc_LamC = Particle(*evt);
                nlamc_daughters = evt -> da(1) - evt -> da(0)+1;
            }
            
            if(evt->mo(0) == idlamc)
            {
                switch(abs(idHEP))
                {
                    case 11:
                        mc_e = Particle(*evt);
                        ehere = 1;
                        break;
                    case 12:
                        mc_nue = Particle(*evt);
                        nuehere = 1;
                        break;
                    case 13:
                        mc_mu = Particle(*evt);
                        muhere = 1;
                        break;
                    case 14:
                        mc_numu = Particle(*evt);
                        numuhere = 1;
                        break;
                    case 111:
                        mc_pi0 = Particle(*evt);
                        pi0here = 1;
                        break;
                    case 211:
                        mc_pi = Particle(*evt);
                        pihere = 1;
                        break;
                    case 3122:
                        mc_lam = Particle(*evt);
                        lam_daid1 = evt->da(0);
                        lam_daid2 = evt->da(1);
                        break;
                    default:
                        return;
                }
            }
            
            if(id == lam_daid1 || id == lam_daid2)
            {
                if(abs(idHEP) == 2212) {pinlam = 1; mc_pfromlam = Particle(*evt);}
                    else if(abs(idHEP) == 211) {piinlam = 1;  mc_pifromlam = Particle(*evt);}
                            else return;
            }
    
            id++;
        }
        if (lam_c_gen==0) return;
        
        if(pinlam && piinlam && lam_daid2-lam_daid1+1==2) lamhere =1;
        
        if(lamhere && ((pihere && !pi0here && !muhere && !numuhere && !ehere && !nuehere && nlamc_daughters==2) 
            || (pihere && pi0here && !muhere && !numuhere && !ehere && !nuehere && nlamc_daughters==3)
            || (!pi0here && !pihere && muhere && numuhere && !ehere && !nuehere && nlamc_daughters==3) 
            || (!pi0here && !pihere && ehere && nuehere && !muhere && !numuhere && nlamc_daughters==3)))
        {
            if(pihere && !muhere && !numuhere && !ehere && !nuehere && !pi0here && nlamc_daughters==2) 
            {   
                mc_lcch=1; 
               // cout << "Lambda pi mode: m_Lambda = " << mc_lam.mass() << "; m_pi = " << mc_pi.mass() << "; Minv = " << (mc_lam.p()+mc_pi.p()).mag() << endl;
            }
            if(pihere && !muhere && !numuhere && !ehere && !nuehere && pi0here && nlamc_daughters==3) 
            {
                mc_lcch=2; 
               // cout << "Lambda pi pi0 mode: m_Lambda = " << mc_lam.mass() << "; m_pi = " << mc_pi.mass() << "; m_pi0 = "<< mc_pi0.mass() << "; Minv = " << (mc_lam.p()+mc_pi.p()+mc_pi0.p()).mag() << endl;
            }
            if(!pihere && !pi0here && ehere && nuehere && !muhere && !numuhere && nlamc_daughters==3) 
            {
                mc_lcch=3;
               // cout << "Lambda e nu_e mode: m_Lambda = " << mc_lam.mass() << "; m_e = " << mc_e.mass() << "; m_nue = "<< mc_nue.mass() << "; Minv = " << (mc_lam.p()+mc_e.p()+mc_nue.p()).mag() << endl;
            }
            if(!pihere && !pi0here && muhere && numuhere && !ehere && !nuehere && nlamc_daughters==3) 
            {
                mc_lcch=4;
               // cout << "Lambda mu nu_mu mode: m_Lambda = " << mc_lam.mass() << "; m_mu = " << mc_mu.mass() << "; m_numu = "<< mc_numu.mass() << "; Minv = " << (mc_lam.p()+mc_mu.p()+mc_numu.p()).mag() << endl;
            }
        }
        else return;
        
        nlamc_gen+=lam_c_gen;
        
        if(!(nevent%1000)) cout << "mc_pUPS: " << mc_pUPS <<  ";  pUPS: "<< pUPS << endl;
        pUPS = mc_pUPS; //temp
        
        //*****************FILLING GEN MC TREE**********************
         
                    
                    
        //lam heli

        Hep3Vector mc_norm_lam_c, mc_norm_lam;
        
        mc_norm_lam_c = (boostT(mc_lam.p(), mc_LamC.p()).vect()).cross(boostT(pUPS, mc_LamC.p()).vect());
        mc_norm_lam_c = mc_norm_lam_c*(1./mc_norm_lam_c.mag());
        
        mc_norm_lam = (boostT(mc_pfromlam.p(), mc_LamC.p()).vect()).cross(boostT(mc_pifromlam.p(), mc_LamC.p()).vect());
        mc_norm_lam = mc_norm_lam*(1./mc_norm_lam.mag());
        
        double mc_dphi_lc_lam=mc_norm_lam_c.angle(mc_norm_lam);
        if(boostT(mc_pfromlam.p(), mc_LamC.p()).vect().dot(mc_norm_lam_c) < 0) mc_dphi_lc_lam = -mc_dphi_lc_lam;
                    
        //W heli and angle chi
        double mc_cosW = -10, mc_angchi=-10;

        HepLorentzVector mc_p_l, mc_p_nu;
        if (mc_lcch==3 || mc_lcch==4)
        {
            
            
            if(mc_lcch==3)
            {
                mc_p_l = mc_e.p();
                mc_p_nu = mc_nue.p();
            }
            else
            {
                mc_p_l = mc_mu.p();
                mc_p_nu = mc_numu.p();
            }
            mc_cosW = -cos(heli(mc_p_l,mc_LamC.p(),mc_LamC.p()-mc_lam.p()));
            
            
            Hep3Vector mc_norm_W;
                    
            mc_norm_W = (boostT(mc_p_nu, mc_LamC.p()).vect()).cross(boostT(mc_p_l, mc_LamC.p()).vect()); 
            mc_norm_W = mc_norm_W*(1./mc_norm_W.mag());
            
            mc_angchi = mc_norm_lam.angle(mc_norm_W);
            if(boostT(mc_p_nu, mc_LamC.p()).vect().dot(mc_norm_lam) < 0) mc_angchi = -mc_angchi;
        }
                    
        //Lam_c 2nd daughter's momentum
        HepLorentzVector mc_P_2d_from_LamC;
        switch(mc_lcch)
        {
            case 1:
                mc_P_2d_from_LamC = mc_pi.p();
                break;
            case 2:
                mc_P_2d_from_LamC = mc_pi.p();
                break;
            case 3:
                mc_P_2d_from_LamC = mc_e.p();
                break;
            case 4:
                mc_P_2d_from_LamC = mc_mu.p();
                break;
            default:
                return;
        }
                    
                    
        //Lam momentum in Lam_c system
        //plamlcs = ;
                    
        t2 -> column("fox", fox);  
        t2 -> column("ecms",pUPS.mag());
        t2 -> column("mc_ecms",mc_pUPS.mag());
        t2 -> column("ch_lamc", mc_LamC.charge());
                
        t2 -> column("lcch",mc_lcch);
        t2 -> column("mlc", mc_LamC.mass());
        t2 -> column("hl", -cos(heli(mc_pfromlam.p(), mc_LamC.p(),  mc_lam.p())));
        t2 -> column("hlc", -cos(heli(mc_lam.p(),pUPS,mc_LamC.p())) );
                
        t2 -> column("q", (mc_LamC.p()-mc_lam.p()).mag());
        t2 -> column("hw",mc_cosW);
        t2 -> column("chi",mc_angchi);
        t2 -> column("plslc",boostT(mc_lam.p(), mc_LamC.p()).vect().mag());
                
        t2 -> column("philclam",mc_dphi_lc_lam);
        t2 -> column("lcp2dcm",pStar(mc_P_2d_from_LamC,elec,posi).vect().mag());
        t2 -> column("lcp2dlab",mc_P_2d_from_LamC.vect().mag());
        
        t2 -> column("plamclab",mc_LamC.p().vect().mag());
        t2 -> column("plamccms",pStar(mc_LamC.p(),elec,posi).vect().mag() );
        t2 -> dumpData();
        
        //***************************************************
        
        if(!(nevent%1000)) cout << "total Lambda_c generated in studied modes (1,2,3,4): " << nlamc_gen << "; total Lambda_c reconstructed: " << nlamc_rec << endl;
        
        
        
        //------------------------EXPERIMENTAL DATA ANALYSIS-----------------------------
        //--------------------------------------------------------------------------------
        //------------------------MAKE PARTICLE LISTINGS----------------------------------
        //protons
        std::vector<Particle> p_p, p_m; 
        makeProton(p_p,p_m);
        
        if(p_p.size()+p_m.size()==0) return; 
        
        
        //kaons and pions
        std::vector<Particle>  k_p, k_m, pi_p, pi_m, pions;
        makeKPi(k_p, k_m, pi_p, pi_m,1);
        
        
        for(std::vector<Particle>::iterator l = pi_m.begin(); l!=pi_m.end(); ++l)
            pions.push_back(*l);
        for(std::vector<Particle>::iterator l = pi_p.begin(); l!=pi_p.end(); ++l)
            pions.push_back(*l);
        
        withKaonId(k_p,0.6,3,1,5,3,4);
        withKaonId(k_p,0.6,3,1,5,3,2);
        withKaonId(k_m,0.6,3,1,5,3,4);
        withKaonId(k_m,0.6,3,1,5,3,2);
        
        ntrk=k_p.size()+k_m.size();

        withdRdZcut(k_p,runIp.z());
        withdRdZcut(pi_p,runIp.z());
        withdRdZcut(k_m,runIp.z());
        withdRdZcut(pi_m,runIp.z());
        
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
        //combination (L_b,ptype_Lamc, p_m, D0_b);
        //setUserInfo(L_b, 1);
        //combination (L_b,ptype_Lamc, p_m, D_m, pi_p);
        //setUserInfo(L_b, 2);
        combination (L_b,ptype_Lamc, p_m, Dst_m, pi_p);
        setUserInfo(L_b, 3);
        combination (L_b,ptype_Lamc, p_m, Dst0_b);
        setUserInfo(L_b, 4);

        
        //combination (L_,ptype_ALamc, p_p, D0);
        //setUserInfo(L_, 1);
        //combination (L_,ptype_ALamc, p_p, D_p, pi_m);
        //setUserInfo(L_, 2);
        combination (L_,ptype_ALamc, p_p, Dst_p, pi_m);
        setUserInfo(L_, 3);
        combination (L_,ptype_ALamc, p_p, Dst0);
        setUserInfo(L_, 4);
        
        
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
        
        if(!(nevent%1000))std::cout<<nevent<<" event. Number of candidates p = " << p_p.size() << "; pbar = " << p_m.size() << "; pi+ = "<< pi_p.size() << "; pi- = "<< pi_m.size() << "; K+ = "<< k_p.size() << "; K- = "<< k_m.size() << "; K_S = "<< k_s.size() << "; pi0 = "<< pi0.size() << "; D0 = " << D0.size() << "; D0bar = " << D0_b.size() << "; D+ = " << D_p.size() << "; D- = "<< D_m.size() << "; Dst0 = " << Dst0.size() << "; Dst0_b = " << Dst0_b.size() << "; D*+ = " << Dst_p.size() << "; D*- = " << Dst_m.size() << "; Lam = " << lam.size() << "; Lam_bar = "<< lamb.size() << "; Lam_c = " << Lc.size() << "; Lam_c_bar = " << Lcb.size() << "; e+ = " << e_p.size() <<"; e- = " << e_m.size() <<"; mu+ = " << mu_p.size() <<"; mu- = " << mu_m.size() << "; Number of recoil candidates L_ = " << L_.size() << "; L_b = " << L_b.size() << '\n';
        
     
        //######################################   FINAL SELECTION
            
        for (std::vector<Particle>::iterator a=L_.begin(); a!=L_.end();++a)
        {
            Particle &ALamC=*a;
            HepLorentzVector momentum=ALamC.p();
            int charge_tag= calcuCharge (&ALamC);
            
            double rmx = (pUPS-momentum).mag(), rm=-1000, Mvis=-1000, rmNopi0 = -1000;//, rm =(pUPS-(momentum+LamC.p())).mag();
            
            if (rmx > 0 && rmx < 2.6) 
            {
                //std::cout<<nevent<<" Selected!" << endl;
                int tag = dynamic_cast<UserInfo&>(ALamC.userInfo()).channel(),
                dstch=-1, dch=-1, lcch=0;
                double mD=-1, mDst=-1, mKs=-1, mL=-1, mLc=-1, hl = -10, hlc = -10, cosW = -10, angchi = -10, 
                pvis=-10, qW=1000, p_2d_from_lamc_cms=-1, p_2d_from_lamc_labs=-1, dphi_lc_lam = -1000, plamlcs = -10, plamclab = -10, plamccms=-10;
        
                if (tag==1)
                {
                    dch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                    mD =  dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                    if( dch==3 || dch==6 ) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                    //mKs = ALamC.child(1).child(0).mass();
                }
                else if (tag==2)
                {
                    dch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                    mD =  dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                    if( dch==2) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                        else if (dch == 3) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(2).userInfo()).mass();
                }
                else if (tag==3)
                {
                    dstch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                    mDst = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                    dch = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).channel();
                    mD = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                    if( dstch==1 && (dch==3 || dch==6)) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(0).userInfo()).mass();
                    else if (dstch==2 && dch==2) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(0).userInfo()).mass();
                    else if (dstch==2 && dch == 3) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(2).userInfo()).mass();
                }
                else if (tag==4)
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
                    lam_c_rec++;
                    
                    
                    int addpi=0, addpi0=0, totcharge=calcuCharge(&(*l))+calcuCharge(&(*a));
                    for(std::vector<Particle>::iterator i=pi0.begin(); i!=pi0.end();++i)
                    {
                        if( !(checkSame(*l,*i)||checkSame(*a,*i)) ) addpi0++;
                    }
                    for(std::vector<Particle>::iterator pi = pions.begin(); pi!=pions.end(); ++pi)
                    {
                        if( !(checkSame(*l,*pi)||checkSame(*a,*pi)) ) addpi++;
                        totcharge += pi->charge();
                    }
                    
                    
                    
                    rm =(pUPS-(momentum+LamC.p())).mag();
                    Mvis = (momentum+LamC.p()).mag();
                    lcch = dynamic_cast<UserInfo&>(LamC.userInfo()).channel(); 
                    mLc = LamC.mass();          
                   
                    pvis = pStar(momentum+LamC.p(),elec,posi).vect().mag();
                    if(lcch==2)
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
                    
                    //W heli and angle chi
                    if (lcch==3 || lcch==4)
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
                    
                    t1 -> column("mc_lcch", mc_lcch);
                    t1 -> column("mc_pnu", pStar(mc_p_nu,elec,posi).vect().mag());
                    t1 -> column("mc_plamc", pStar(mc_LamC.p(),elec,posi).vect().mag());
                    t1 -> column("mc_angnv", pStar(mc_p_nu,elec,posi).vect().angle(-pStar(momentum+LamC.p(),elec,posi).vect()));  
                    t1 -> column("mcanglcx", pStar(mc_LamC.p(),elec,posi).vect().angle(-pStar(momentum,elec,posi).vect())); 
                    t1 -> column("mc_rmx",(mc_pUPS-momentum).mag());
                    
                    t1 -> column("addpi",addpi);
                    t1 -> column("addpi0",addpi0);
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
                    t1 -> column("mc_lcch", mc_lcch);
                    t1 -> column("mc_pnu", -1);
                    t1 -> column("mc_plamc", -1);
                    t1 -> column("mc_angnv", -10);
                    t1 -> column("mcanglcx", -10);
                    t1 -> column("mc_rmx",(mc_pUPS-momentum).mag());
                    
                    t1 -> column("addpi",-1);
                    t1 -> column("addpi0",-1);
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
            
            if (rmx > 0 && rmx < 2.6) 
            {
                int tag = dynamic_cast<UserInfo&>(ALamC.userInfo()).channel(),
                dstch=-1, dch=-1, lcch=0;
                double mD=-1, mDst=-1, mKs=-1, mL=-1, mLc=-1, hl = -10, hlc = -10, cosW = -10, angchi = -10, 
                       pvis=-10, qW=1000, p_2d_from_lamc_cms=-1, p_2d_from_lamc_labs=-1, dphi_lc_lam = -1000, plamlcs = -10, plamclab = -10, plamccms=-10;
        
                if (tag==1)
                {
                    dch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                    mD =  dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                    if( dch==3 || dch==6 ) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                    //mKs = ALamC.child(1).child(0).mass();
                }
                else if (tag==2)
                {
                    dch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                    mD =  dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                     if( dch==2) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                        else if (dch == 3) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(2).userInfo()).mass();
                }
                else if (tag==3)
                {
                    dstch = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).channel();
                    mDst = dynamic_cast<UserInfo&>(ALamC.child(1).userInfo()).mass();
                    dch = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).channel();
                    mD = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).userInfo()).mass();
                    if( dstch==1 && (dch==3 || dch==6)) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(0).userInfo()).mass();
                    else if (dstch==2 && dch==2) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(0).userInfo()).mass();
                    else if (dstch==2 && dch == 3) mKs = dynamic_cast<UserInfo&>(ALamC.child(1).child(0).child(2).userInfo()).mass();
                }
                else if (tag==4)
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
                    lam_c_rec++;
                    
                    int addpi=0, addpi0=0, totcharge=calcuCharge(&(*l))+calcuCharge(&(*a));
                    for(std::vector<Particle>::iterator i=pi0.begin(); i!=pi0.end();++i)
                    {
                        if( !(checkSame(*l,*i)||checkSame(*a,*i)) ) addpi0++;
                    }
                    for(std::vector<Particle>::iterator pi = pions.begin(); pi!=pions.end(); ++pi)
                    {
                        if( !(checkSame(*l,*pi)||checkSame(*a,*pi)) ) addpi++;
                        totcharge += pi->charge();
                    }
                    
                    rm =(pUPS-(momentum+LamC.p())).mag();
                    Mvis = (momentum+LamC.p()).mag();
                    

                    lcch = dynamic_cast<UserInfo&>(LamC.userInfo()).channel(); 
                    mLc = LamC.mass();          
                    pvis = pStar(momentum+LamC.p(),elec,posi).vect().mag();
                    if(lcch==2)
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
                    if (lcch==3 || lcch==4)
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
                    
                    t1 -> column("mc_lcch", mc_lcch);
                    t1 -> column("mc_pnu", pStar(mc_p_nu,elec,posi).vect().mag());
                    t1 -> column("mc_plamc", pStar(mc_LamC.p(),elec,posi).vect().mag());
                    t1 -> column("mc_angnv", pStar(mc_p_nu,elec,posi).vect().angle(-pStar(momentum+LamC.p(),elec,posi).vect()));  
                    t1 -> column("mcanglcx", pStar(mc_LamC.p(),elec,posi).vect().angle(-pStar(momentum,elec,posi).vect()));
                    t1 -> column("mc_rmx",(mc_pUPS-momentum).mag());
                    
                    t1 -> column("addpi",addpi);
                    t1 -> column("addpi0",addpi0);
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
                    
                    t1 -> column("mc_lcch", mc_lcch);
                    t1 -> column("mc_pnu", -1);
                    t1 -> column("mc_plamc", -1);
                    t1 -> column("mc_angnv", -10);
                    t1 -> column("mcanglcx", -10);
                    t1 -> column("mc_rmx",(mc_pUPS-momentum).mag());
                    
                    t1 -> column("addpi",-1);
                    t1 -> column("addpi0",-1);
                    t1 -> column("totcharg",-100);
                    t1->dumpData();
                }
                
            }
        }
        nlamc_rec+=lam_c_rec;
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
