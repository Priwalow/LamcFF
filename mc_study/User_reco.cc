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
        t1 = BASF_Histogram->ntuple ("data_mc","fox ecms lcch ml mlc ch_lamc hl hlc philclam plamclab plamccms coslclab coslccms mc_lcch mc_pnu mcplccms mcclccms mcplclab mcclclab pflidh piflidh lflcidh d2flcidh d3flcidh pflmidh piflmidh lflcmidh d2lcmidh lcmidh ppiflsm l2dlcsm nlda nlcda mcq mchw mcchi mcplslc mcpnulc"); // not ALL momenta in CMS! 	lepton cosTheta in CMS, rholam, rholamcms	
        t2 = BASF_Histogram->ntuple ("gen_mc","fox ecms mc_ecms mlc ch_lamc lcch hlc hl q hw chi lcp2dcm lcp2dlab  philclam plslc pnulc plamclab plamccms coslclab coslccms" ); // not ALL momenta in CMS!
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
    void MysetGenHepEvtInfoCharged(vector<Particle> &plist );
    bool MysetGenHepEvtInfoCharged(Particle & chg );
    void MysetGenHepEvtInfoLambda(vector<Particle> & vlam );
     
    
    
    
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
        double mc_elec, mc_posi;
        mc_elec = mc_pUPS.px()/sin(0.022);
        mc_posi = mc_pUPS.e()-mc_elec;
        
        
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
            //if(pihere && !muhere && !numuhere && !ehere && !nuehere && pi0here && nlamc_daughters==3) 
            //{
            //    mc_lcch=2; 
               // cout << "Lambda pi pi0 mode: m_Lambda = " << mc_lam.mass() << "; m_pi = " << mc_pi.mass() << "; m_pi0 = "<< mc_pi0.mass() << "; Minv = " << (mc_lam.p()+mc_pi.p()+mc_pi0.p()).mag() << endl;
            //}
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
        elec = mc_elec; // because in MC elec and posi are constant 8 and 3.5 (energy_cms = Y(4S) mass 10.58 GeV)
        posi = mc_posi; //
        
        
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
        t2 -> column("pnulc",boostT(mc_p_nu, mc_LamC.p()).vect().mag());        
        
        t2 -> column("philclam",mc_dphi_lc_lam);
        t2 -> column("lcp2dcm",pStar(mc_P_2d_from_LamC,elec,posi).vect().mag());
        t2 -> column("lcp2dlab",mc_P_2d_from_LamC.vect().mag());
        
        t2 -> column("plamclab",mc_LamC.p().vect().mag());
        t2 -> column("plamccms",pStar(mc_LamC.p(),elec,posi).vect().mag() );
        t2 -> column("coslclab",mc_LamC.p().vect().cosTheta());
        t2 -> column("coslccms",pStar(mc_LamC.p(),elec,posi).vect().cosTheta() );
        t2 -> dumpData();
        
        //***************************************************
        
        if(!(nevent%1000)) cout << "total Lambda_c generated in studied modes (1,2,3,4): " << nlamc_gen << "; total Lambda_c reconstructed: " << nlamc_rec << endl;
        
        
        
        //------------------------EXPERIMENTAL DATA ANALYSIS-----------------------------
        //--------------------------------------------------------------------------------
        //------------------------MAKE PARTICLE LISTINGS----------------------------------
        
        //kaons and pions
        std::vector<Particle>  k_p, k_m, pi_p, pi_m, pions;
        makeKPi(k_p, k_m, pi_p, pi_m,1);
        
        
        MysetGenHepEvtInfoCharged(pi_p);
        MysetGenHepEvtInfoCharged(pi_m);
        
        
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
        
        //Lambda hyperons
        std::vector<Particle> lam, lamb;
        makeLam(lam,lamb);
        
        setUserInfo(lam,  11 ); 
        setUserInfo(lamb, 12 ); 
        doMassVertexFit(lam);
        doMassVertexFit(lamb);
        
        MysetGenHepEvtInfoLambda(lam);
        MysetGenHepEvtInfoLambda(lamb);

        
        
        //leptons (e, mu) 
        std::vector<Particle>  e_p,e_m,mu_p,mu_m;
        makeLepton(e_p,e_m,mu_p,mu_m);
        
        MysetGenHepEvtInfoCharged(mu_p);
        MysetGenHepEvtInfoCharged(mu_m);
        MysetGenHepEvtInfoCharged(e_p);
        MysetGenHepEvtInfoCharged(e_m);
        
        
        withMuId(mu_p);
        withMuId(mu_m);
        withEId(e_p);
        withEId(e_m);
        
        
        //######################################    SIGNAL SIDE
        
        std::vector <Particle> Lc, Lcb, allLamc; 
        
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
        
        //combination (Lc,ptype_Lamc, lam, pi_p, pi0, 0.08);
        //combination (Lcb,ptype_Lamc, lamb, pi_m, pi0, 0.08);
        //setUserInfo(Lc,2);
        //setUserInfo(Lcb,2);
        
        //combination (Lc,ptype_Lamc, p_p, k_m, pi_p, 0.05);
        //combination (Lcb,ptype_Lamc, p_m, k_p, pi_m, 0.05);
        //setUserInfo(Lc,5);
        //setUserInfo(Lcb,5);
        
        
        for(std::vector<Particle>::iterator l = Lc.begin(); l!=Lc.end(); ++l)
            allLamc.push_back(*l);
        for(std::vector<Particle>::iterator l = Lcb.begin(); l!=Lcb.end(); ++l)
            allLamc.push_back(*l);    
        
        
        //######################################   FINAL SELECTION
        
        for(std::vector<Particle>::iterator l=allLamc.begin(); l!=allLamc.end();++l)
        {
            Particle &LamC=*l;
            lam_c_rec++;
            

            int lcch;
            double mL, mLc, hlc=-10, hl=-10, dphi_lc_lam=-10;
            
            lcch = dynamic_cast<UserInfo&>(LamC.userInfo()).channel(); 
            mLc = LamC.mass();          
            mL  = dynamic_cast<UserInfo&>(LamC.child(0).userInfo()).mass();
            
            
            // lamc heli
            HepLorentzVector p_lamc;
            if (lcch==1 || lcch==2 || lcch==5) p_lamc=LamC.p();
                else p_lamc = mc_LamC.p();
            
            hlc = -cos(heli(LamC.child(0).p(),pUPS,p_lamc)); 
            
            //lam heli
            Particle proton_from_lam, pi_from_lam;
            HepLorentzVector p_proton_from_lam, p_pi_from_lam; 
            if(lcch!=5)
            {
                
                if (abs(LamC.child(0).child(0).lund())>1000)
                {
                    proton_from_lam = LamC.child(0).child(0);
                    pi_from_lam = LamC.child(0).child(1);
                    p_proton_from_lam=proton_from_lam.p(); 
                    p_pi_from_lam=pi_from_lam.p();
                }
                else
                {
                    proton_from_lam = LamC.child(0).child(1);
                    pi_from_lam = LamC.child(0).child(0);
                    p_proton_from_lam=proton_from_lam.p(); 
                    p_pi_from_lam=pi_from_lam.p();
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
            
            
            //test match
            int  proton_from_lam_idhep=-1000, pion_from_lam_idhep=-1000, lam_from_lamc_idhep=-1000, d2_from_lamc_idhep=-1000, d3_from_lamc_idhep=-1000,
            proton_from_lam_moidhep=-1000, pion_from_lam_moidhep=-1000, lam_from_lamc_moidhep=-1000, d2_from_lamc_moidhep=-1000, lamc_moidhep=-1000, 
            ppi_from_same_mother=-1000, lam_d2lc_same_mother=-1000, lam_d3lc_same_mother=-1000, d2lc_d3lc_same_mother=-1000, n_lamc_daughters= -1000, n_lam_daughters=-1000;
            
            
            //cout << "Matching p from lam" << endl;
            Mdst_sim_trk_Manager &xrefMgr_ = Mdst_sim_trk_Manager::get_manager();
            for(std::vector<Mdst_sim_trk>::iterator i = xrefMgr_.begin(); i != xrefMgr_.end(); ++i) 
            {
                if (proton_from_lam.mdstCharged().trk() != i->trk()) continue;
                proton_from_lam_idhep = i->hepevt().idhep();
                if (!(i->hepevt().mother())) continue;
                proton_from_lam_moidhep = i->hepevt().mother().idhep();
                
            }
            
            //cout << "Matching pi from lam" << endl;
            for(std::vector<Mdst_sim_trk>::iterator i = xrefMgr_.begin(); i != xrefMgr_.end(); ++i) 
            {
                if (pi_from_lam.mdstCharged().trk() != i->trk()) continue;
                pion_from_lam_idhep = i->hepevt().idhep();
                if (!(i->hepevt().mother())) continue;
                pion_from_lam_moidhep = i->hepevt().mother().idhep();
                
            }
            
            //cout << "Matching 2d from lamc" << endl;
            for(std::vector<Mdst_sim_trk>::iterator i = xrefMgr_.begin(); i != xrefMgr_.end(); ++i) 
            {
                if (LamC.child(1).mdstCharged().trk() != i->trk()) continue;
                d2_from_lamc_idhep = i->hepevt().idhep();
                if (!(i->hepevt().mother())) continue;
                d2_from_lamc_moidhep = i->hepevt().mother().idhep();
                
            }
            
            if(proton_from_lam.genHepevt() && pi_from_lam.genHepevt())
            {
                if(proton_from_lam.genHepevt().mother()==pi_from_lam.genHepevt().mother()) 
                    ppi_from_same_mother = 1;
                else ppi_from_same_mother = 0;
            }
            
            //cout << "Matching lam from lamc" << endl;
            if(LamC.child(0).genHepevt())
            {
                //cout << "Lam" << endl;
                lam_from_lamc_idhep = LamC.child(0).genHepevt().idhep();
                if(LamC.child(0).genHepevt().mother())
                {
                    //cout << "Lam mother" << endl;
                    lam_from_lamc_moidhep = LamC.child(0).genHepevt().mother().idhep();
                    if(LamC.child(0).genHepevt() && LamC.child(1).genHepevt())
                    {
                        if(LamC.child(0).genHepevt().mother() == LamC.child(1).genHepevt().mother())
                            lam_d2lc_same_mother=1;
                        else lam_d2lc_same_mother=0;
                    }
                        
                        
                    int first_lamc_da_ID = LamC.child(0).genHepevt().mother().daFirst(), last_lamc_da_ID = LamC.child(0).genHepevt().mother().daLast();
                    n_lamc_daughters = last_lamc_da_ID-first_lamc_da_ID+1;
                    
                    
                    
                    
                    if(n_lamc_daughters == 3 && lam_d2lc_same_mother==1)
                    {
                        int tempid = 0, tempidhep = -1000, sametype=0;
                        
                        for (std::vector<Gen_hepevt>::iterator evt = evt_manager.begin(); evt != evt_manager.end(); ++evt) 
                        {
                            tempid++;
                            if (tempid>last_lamc_da_ID) break;
                            if (tempid<first_lamc_da_ID) continue;
                            tempidhep = evt->idhep();
                            if ((tempidhep==lam_from_lamc_idhep || tempidhep==d2_from_lamc_idhep) && !sametype) //3d particle should be neutral
                            {
                                sametype++;
                                continue;
                            }
                            else 
                            {
                                d3_from_lamc_idhep=tempidhep;
                                break;
                            }
                            
                        }
                    }
                    if(LamC.child(0).genHepevt().mother().mother())
                    {
                        lamc_moidhep = LamC.child(0).genHepevt().mother().mother().idhep();
                    }
                }
                n_lam_daughters = LamC.child(0).genHepevt().daLast()-LamC.child(0).genHepevt().daFirst()+1;
            }
            /*
            cout << "Matching lamc" << endl;
            lamc_idhep = LamC.relation().genHepevt().idhep();
            lamc_ID = LamC.relation().genHepevt().get_ID();
            if(LamC.relation().genHepevt().mother()) lamc_moidhep = LamC.child(0).relation().genHepevt().mother().idhep();
            
            cout << "Matching 3d from lamc" << endl;
            last_lamc_da_ID = LamC.relation().genHepevt().daLast();
            first_lamc_da_ID = LamC.relation().genHepevt().daFirst();
            n_lamc_daughters = last_lamc_da_ID - first_lamc_da_ID + 1;
            if(n_lamc_daughters == 3)
            {
                int tempid = 0, tempidhep = -1000;
                for (std::vector<Gen_hepevt>::iterator evt = evt_manager.begin(); evt != evt_manager.end(); ++evt) 
                {
                    tempid++;
                    if (tempid<first_lamc_da_ID) continue;
                    tempidhep = evt->idhep();
                    if (abs(tempidhep)==3122 || abs(tempidhep)==11 || abs(tempidhep)==13) continue;
                    else 
                    {
                        d3_from_lamc_idhep=tempidhep;
                        break;
                    }
                    
                }
            }

           
            
             cout << "Matching" << endl;
            if (proton_from_lam.relation().genHepevt()) lam_flag++;
            cout << "1" << endl;
            if (pi_from_lam.relation().genHepevt()) lam_flag++;
            cout << "2" << endl;
            if (LamC.child(1).relation().genHepevt()) lamc_flag++;
            cout << "3" << endl;
           if (proton_from_lam.relation().genHepevt().mother() && pi_from_lam.relation().genHepevt().mother()) 
            { // mothers
                if (proton_from_lam.relation().genHepevt().mother()==pi_from_lam.relation().genHepevt().mother()) 
                { // eq mothers
                    if (abs(proton_from_lam.relation().genHepevt().mother().idhep())==3122)
                    {
                        lam_flag++;
                        lamc_flag++;
                    }
                }    
            }
            
 
            if (LamC.child(0).relation().genHepevt().mother() && LamC.child(1).relation().genHepevt().mother()) 
            { // mothers
                if (LamC.child(0).relation().genHepevt().mother()==LamC.child(1).relation().genHepevt().mother()) 
                { // eq mothers
                    if (abs(LamC.child(0).relation().genHepevt().mother().idhep())==4122 && LamC.child(0).relation().genHepevt().mother()==LamC.relation().genHepevt())
                    {
                        lamc_flag++;
                    }
                }
            }*/
            
            
            t1 -> column("fox", fox);  
            t1 -> column("ecms",pUPS.mag());
            
            t1 -> column("lcch",lcch);
            t1 -> column("ml", mL);
            t1 -> column("mlc", mLc);
            t1 -> column("ch_lamc", LamC.charge());
            
            t1->column("hl", hl);
            t1 -> column("hlc", hlc);
            t1 -> column("philclam",dphi_lc_lam);
            
            t1 -> column("plamclab",p_lamc.vect().mag());
            t1 -> column("plamccms",pStar(p_lamc,elec,posi).vect().mag());
            t1 -> column("coslclab",p_lamc.vect().cosTheta());
            t1 -> column("coslccms",pStar(p_lamc,elec,posi).vect().cosTheta() );
            
            t1 -> column("mc_lcch", mc_lcch);
            t1 -> column("mc_pnu", pStar(mc_p_nu,elec,posi).vect().mag());
            t1 -> column("mcplccms", pStar(mc_LamC.p(),elec,posi).vect().mag());
            t1 -> column("mcclccms",pStar(mc_LamC.p(),elec,posi).vect().cosTheta());
            t1 -> column("mcplclab",mc_LamC.p().vect().mag());
            t1 -> column("mcclclab",mc_LamC.p().vect().cosTheta());
            
            t1 -> column("pflidh",proton_from_lam_idhep);
            t1 -> column("piflidh",pion_from_lam_idhep);
            t1 -> column("lflcidh", lam_from_lamc_idhep);
            t1 -> column("d2flcidh",d2_from_lamc_idhep);
            t1 -> column("d3flcidh",d3_from_lamc_idhep);
           
            t1 -> column("pflmidh", proton_from_lam_moidhep);
            t1 -> column("piflmidh",pion_from_lam_moidhep);
            t1 -> column("lflcmidh",lam_from_lamc_moidhep);
            t1 -> column("d2lcmidh",d2_from_lamc_moidhep);
            t1 -> column("lcmidh", lamc_moidhep);
            
            t1 -> column("ppiflsm",ppi_from_same_mother);
            t1 -> column("l2dlcsm",lam_d2lc_same_mother);
            
            t1 -> column("nlcda", n_lamc_daughters);
            t1 -> column("nlda",n_lam_daughters);
            
            t1 -> column("mcq", (mc_LamC.p()-mc_lam.p()).mag());
            t1 -> column("mchw",mc_cosW);
            t1 -> column("mcchi",mc_angchi);
            t1 -> column("mcplslc",boostT(mc_lam.p(), mc_LamC.p()).vect().mag());
            t1 -> column("mcpnulc",boostT(mc_p_nu, mc_LamC.p()).vect().mag()); 

            t1->dumpData();
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
    
    
    
    void  MysetGenHepEvtInfoCharged(vector<Particle> &plist )
    {
        if ( plist.size() == 0 ) { return; }
        
        // loop over the plist
        for ( vector<Particle>::iterator i = plist.begin();
             i != plist.end(); ++i ) {
            MysetGenHepEvtInfoCharged( (Particle &)*i );
             }
    }
    
    bool MysetGenHepEvtInfoCharged(Particle & chg )
    {
        //
        bool bret = false;
        if ( !chg.mdstCharged() ) { return bret; }
        
        const Gen_hepevt & gen = get_hepevt(chg.mdstCharged());
        if ( gen ) {
            chg.relation().genHepevt( gen );
            bret = true;
        }
        
        return bret;
    }
    
    
    
    
    
    
    
    void MysetGenHepEvtInfoLambda(vector<Particle> & vlam )
    { 
        if ( vlam.size() <= 0 ) { return; }
        
        for ( vector<Particle>::iterator i = vlam.begin(); i!= vlam.end(); i++ ) 
        {
            
            // 
            if ( abs((*i).lund()) != 3122 ) { continue; }
            
            // set children pion
            if ( (*i).relation().nChildren() != 2 ) { continue; }
            Particle & pos = (*i).child(0);
        Particle & neg = (*i).child(1); 
        if ( !pos || !neg ) { continue; }
        
        // set Gen_hepevt on children
        bool setGen1 = MysetGenHepEvtInfoCharged( pos );
        bool setGen2 = MysetGenHepEvtInfoCharged( neg );
        if ( !setGen1 || !setGen2 ) { continue; }
        
        // get Gen_hepevt of children
        const Gen_hepevt & gen1 = pos.genHepevt();
        const Gen_hepevt & gen2 = neg.genHepevt();
        if ( !gen1 || !gen2 ) { continue; }
        
        // check children 
        if ((*i).lund() > 0){
            if ( gen1.idhep() != 2212 || gen2.idhep() != -211  ) { continue; }
        }else{ 
            if ( gen1.idhep() != 211  || gen2.idhep() != -2212 ) { continue; }
        }
        
        // get mother Gen_hepevt of chilren
        const Gen_hepevt & mom1 = gen1.mother();
        const Gen_hepevt & mom2 = gen2.mother();
        if ( !mom1 || !mom2 ) { continue; }
        
        // mom1 == mom2 == lambda
        if ( mom1.get_ID() == mom2.get_ID() &&
            abs(mom1.idhep()) == 3122 ) {
            (*i).relation().genHepevt( mom1 );
        break;
            }
            
        }
    }
    
}
