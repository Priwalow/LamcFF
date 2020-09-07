#include "User_reco.h"
#include "userinfo.h"
#include EVTCLS_H //R2 distribution



using namespace Belle;
namespace Belle {
    void User_reco::hist_def( void )
    {
        
        //  extern BelleTupleManager* BASF_Histogram;
        //  t1 = BASF_Histogram->ntuple ("1","ch bestpi0 ml lcch rm rm1 npi mlc npi0 ml1  hlc hl" );   
    };
    //***********************************************************************
    void doMassVertexFit(class vector<Particle> &p_list, double mass=-1);
    void doMassVertexFit(class Particle &P, double mass=-1);
    void doVertexFit(class vector<Particle> &p_list);
    void doVertexFit(class Particle &P);
    void setPi0Error(Particle &p);
    void setPi0Error(std::vector<Particle> &p_list);
    double heli (  HepLorentzVector  * P4A, HepLorentzVector * P4B, HepLorentzVector* P4system);
    double heli (  HepLorentzVector   P4A, HepLorentzVector P4B, HepLorentzVector P4system);
    HepLorentzVector boostT( HepLorentzVector p, HepLorentzVector p_boost);
    HepLorentzVector boostT( HepLorentzVector *p, HepLorentzVector *p_boost);
    void eraseLeptons(std::vector<Particle> &list);
    
    
    void withdRdZcut(class std::vector<Particle> &p_list,double ip_position=0., double drcut = 2., double dzcut = 4.);
    void makeLam(std::vector<Particle> &lam0, std::vector<Particle> &lam0b);
    void makeProton (std::vector<Particle> &p, std::vector<Particle> &pbar);
    
    
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
        
        Evtcls_hadron_info_Manager&  ehimgr =
        Evtcls_hadron_info_Manager::get_manager();
        
        
        
        std::vector<Evtcls_hadron_info>::iterator iti = ehimgr.begin();
        const HepPoint3D&   runIp     = IpProfile::position();
        
        double fox=0;
        int ntrk=0;
        /*if( iti !=  ehimgr.end() && *iti )
         *     {
         *       fox = (*iti).R2();
         *       //       ntrk = (*iti).Ntrk();
    }
    */
        
        
        //------------------------MAKE PARTICLE LISTINGS----------------------------------
        
        //protons
        

        
        std::vector<Particle> p_p, p_m; 
        makeProton(p_p,p_m);
        
        if(p_p.size()+p_m.size()==0) return; 
        
        if (!(nevent%1000))std::cout<<nevent<<" p: "<< p_p.size() << "  anti-p: " << p_m.size() << '\n';
        

        //kaons and pions
        std::vector<Particle>  k_p, k_m, pi_p, pi_m, pions;
        makeKPi(k_p, k_m, pi_p, pi_m,1);
        
        ntrk=k_p.size()+k_m.size();
        
        /*
         *   withdRdZcut(k_p,runIp.z());
         *   withdRdZcut(pi_p,runIp.z());
         *   withdRdZcut(k_m,runIp.z());
         *   withdRdZcut(pi_m,runIp.z());
         */
        
        withKaonId(k_p,0.6,3,1,5,3,4);
        withKaonId(k_p,0.6,3,1,5,3,2);
        withKaonId(k_m,0.6,3,1,5,3,4);
        withKaonId(k_m,0.6,3,1,5,3,2);
        
        withPionId(pi_p,0.6,3,1,5,2,3);
        withPionId(pi_p,0.6,3,1,5,2,4);
        withPionId(pi_m,0.6,3,1,5,2,3);
        withPionId(pi_p,0.6,3,1,5,2,4);
        
        eraseLeptons(k_p);
        eraseLeptons(k_m);
        eraseLeptons(pi_p);
        eraseLeptons(pi_m);
        
        for(std::vector<Particle>::iterator l = pi_m.begin(); l!=pi_m.end(); ++l)
            pions.push_back(*l);
        for(std::vector<Particle>::iterator l = pi_p.begin(); l!=pi_p.end(); ++l)
            pions.push_back(*l);
        
        if (!(nevent%1000))std::cout<<nevent<<" pi_p: "<< pi_p.size() << "  pi_m: " << pi_m.size() << '\n';
        if (!(nevent%1000))std::cout<<nevent<<" k_p: "<< k_p.size() << "  k_m: " << k_m.size() << '\n';
        
        


        //Ks mesons
        std::vector<Particle> k_s;
        makeKs(k_s);
        for(std::vector<Particle>::iterator l = k_s.begin(); l!=k_s.end(); ++l)
        {
            if (abs(l->mass()-0.4976)>0.0075)
            {   
                k_s.erase(l); 
                --l;
                continue;
            }
            
            HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
            Vector3 Pt(l->px(),l->py(),0);
            double Ptot = pStar(l->p(),elec,posi).vect().mag();
            V=V-runIp;
            V.setZ(0.);
            
            if (Ptot<0.5)
            {
                if(V.angle(Pt)>0.3 || l->mdstVee2().z_dist()>0.8)
                {
                    k_s.erase(l); 
                    --l;
                    continue;
                }
            }
                else if (Ptot<1.5)
                {
                    if(V.perp()<0.08 || V.angle(Pt)>0.1 || l->mdstVee2().z_dist()>1.8)
                    {
                        k_s.erase(l); 
                        --l;
                        continue;
                    }
                }
                    else 
                    {
                        if(V.perp()<0.22 || V.angle(Pt)>0.03 || l->mdstVee2().z_dist()>2.4)
                        {
                            k_s.erase(l); 
                            --l;
                            continue;
                        }
                    }
                
        }
        doMassVertexFit(k_s);
        double k_s_chisq, bufchisq=1000000;
        while(k_s.size()>1)
        {
            for(std::vector<Particle>::iterator l = k_s.begin(); l!=k_s.end(); ++l)
            {
                k_s_chisq = dynamic_cast<UserInfo&>(l->userInfo()).vchisq();
                if(k_s_chisq > bufchisq)
                {
                    k_s_chisq = bufchisq;
                    k_s.erase(l); 
                    --l;
                    continue;
                }
                bufchisq = k_s_chisq;
            }
        }

        if (!(nevent%1000))std::cout<<nevent<<" k_s: " << k_s.size() << "; chisq/ndf = " << k_s_chisq << '\n';

        //Pi0 mesons
        std::vector<Particle>  pi0;
        makePi0(pi0);
        
        for(std::vector<Particle>::iterator i=pi0.begin(); i!=pi0.end();++i)
        if(abs(i->mdstPi0().mass()-.135)>0.015)
        {pi0.erase(i); --i;}
            
        //leptons    
       /* std::vector<Particle>  e_p,e_m,mu_p,mu_m;
        makeLepton(e_p,e_m,mu_p,mu_m);
        
        withMuId(mu_p);
        withMuId(mu_m);
        withEId(e_p);
        withEId(e_m);
        */
        if (!(nevent%1000))std::cout<<nevent<<" npi0: " << pi0.size() << '\n';
        


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
   
        
        if (!(nevent%1000))std::cout<<nevent<<" ngamma: " << photons.size() << '\n';
        
        std::vector<Particle> D0, D0_b, D_p, D_m, Dst_p, Dst_m, Dst0, Dst0_b, pi_pm;
        
        
        //######################################    TAG SIDE
        combination (pi_pm, ptype_D0B, pi_p, pi_m);
        
        cout << "Making D0_b" << endl;
        combination (D0_b,ptype_D0B, k_p, pi_m, 0.06);
        setUserInfo(D0_b, 1); 
        combination (D0_b,ptype_D0B, k_p, pi_m, pi_pm, 0.05);
        setUserInfo(D0_b, 2); 
        combination (D0_b,ptype_D0B, k_s, pi_pm, 0.05);
        setUserInfo(D0_b, 3);
        combination (D0_b,ptype_D0B, k_p, pi_m, pi0, 0.06);
        setUserInfo(D0_b, 4);
        
        
        combination (D0_b,ptype_D0B, k_p, pi_m, pi0, pi_pm, 0.06);
        setUserInfo(D0_b, 5);
        combination (D0_b,ptype_D0B, k_s, pi_pm, pi0, 0.06);
        setUserInfo(D0_b, 6);
        
        cout << "Making D_m" << endl;
        combination (D_m,ptype_Dm, k_p, pi_m, pi_m, 0.05);
        setUserInfo(D_m, 1);
        combination (D_m,ptype_Dm, k_s, pi_m, 0.05);
        setUserInfo(D_m, 2);
        combination (D_m,ptype_Dm, k_s, pi_m, pi_pm, 0.05);
        setUserInfo(D_m, 3);
        combination (D_m,ptype_Dm, k_p, k_m, pi_m, 0.05);
        setUserInfo(D_m, 4);
        
        cout << "Selecting D0_b among "<< D0_b.size() << " candidates" << endl;
        doMassVertexFit(D0_b);
        double d0_chisq;
        bufchisq=1000000;
        while(D0_b.size()>1)
        {
            for(std::vector<Particle>::iterator l = D0_b.begin(); l!=D0_b.end(); ++l)
            {
                d0_chisq = dynamic_cast<UserInfo&>(l->userInfo()).vchisq();
                cout << d0_chisq << endl; 
                if((d0_chisq<0) || (d0_chisq > bufchisq))
                {
                    d0_chisq = bufchisq;
                    D0_b.erase(l); 
                    --l;
                    continue;
                }
                bufchisq = d0_chisq;
            }
            cout << D0_b.size() << endl; 
        }
        if (!(nevent%1000))std::cout<<nevent<<" d0_b: " << D0_b.size() << "; chisq/ndf = " << d0_chisq << '\n';
        
        cout << "Selecting D_m among "<< D_m.size() << " candidates" << endl;
        doMassVertexFit(D_m);
        double d_m_chisq;
        bufchisq=1000000;
        while(D_m.size()>1)
        {
            for(std::vector<Particle>::iterator l = D_m.begin(); l!=D_m.end(); ++l)
            {
                d_m_chisq = dynamic_cast<UserInfo&>(l->userInfo()).vchisq();
                if(d_m_chisq > bufchisq)
                {
                    d_m_chisq = bufchisq;
                    D_m.erase(l); 
                    --l;
                    continue;
                }
                bufchisq = d_m_chisq;
            }
        }
        if (!(nevent%1000))std::cout<<nevent<<" d_m: " << D_m.size() << "; chisq/ndf = " << d_m_chisq << '\n';
        
        
                
        cout << "Making Dst0_b" << endl;
        combination (Dst0_b,ptype_Dst0, D0_b, pi0, 0.2);
        setUserInfo(Dst0_b, 1);
        combination (Dst0_b,ptype_Dst0, D0_b, photons, 0.2);
        setUserInfo(Dst0_b, 2);
        
        cout << "Making Dst_m" << endl;
        combination (Dst_m,ptype_Dstm, D0_b, pi_m, 0.2);
        setUserInfo(Dst_m, 1);
        combination (Dst_m,ptype_Dstm, D_m, pi0, 0.2);
        setUserInfo(Dst_m, 2);
        
        cout << "Selecting Dst0_b among "<< Dst0_b.size() << " candidates" << endl;
        for (std::vector<Particle>::iterator i=Dst0_b.begin(); i!=Dst0_b.end();++i)
        {
            Particle dst0b(*i);
            if(abs(dst0b.mass()-dst0b.child(0).mass())>0.025)
            {
                Dst0_b.erase(i); 
                --i;
            }
        }
        
        doMassVertexFit(Dst0_b);
        double dst0_chisq;
        bufchisq=1000000;
        while(Dst0_b.size()>1)
        {
            for(std::vector<Particle>::iterator l = Dst0_b.begin(); l!=Dst0_b.end(); ++l)
            {
                dst0_chisq = dynamic_cast<UserInfo&>(l->userInfo()).vchisq();
                if(dst0_chisq > bufchisq)
                {
                    dst0_chisq = bufchisq;
                    Dst0_b.erase(l); 
                    --l;
                    continue;
                }
                bufchisq = dst0_chisq;
            }
        }
        if (!(nevent%1000))std::cout<<nevent<<" dst0_b: " << Dst0_b.size() << "; chisq/ndf = " << dst0_chisq << '\n';
        
        
        cout << "Selecting Dst_m among "<< Dst_m.size() << " candidates" << endl;
        for (std::vector<Particle>::iterator i=Dst_m.begin(); i!=Dst_m.end();++i)
        {
            Particle dstm(*i);
            if(abs(dstm.mass()-dstm.child(0).mass())>0.025)
            {
                Dst_m.erase(i); 
                --i;
            }
        }
        
        doMassVertexFit(Dst_m);
        double dstm_chisq;
        bufchisq=1000000;
        while(Dst_m.size()>1)
        {
            for(std::vector<Particle>::iterator l = Dst_m.begin(); l!=Dst_m.end(); ++l)
            {
                dstm_chisq = dynamic_cast<UserInfo&>(l->userInfo()).vchisq();
                if(dstm_chisq > bufchisq)
                {
                    dstm_chisq = bufchisq;
                    Dst_m.erase(l); 
                    --l;
                    continue;
                }
                bufchisq = dstm_chisq;
            }
        }
        if (!(nevent%1000))std::cout<<nevent<<" dst_m: " << Dst_m.size() << "; chisq/ndf = " << dstm_chisq << '\n';
        
        cout << "Making recoil" << endl;
        std::vector <Particle> L_, L_b;
        combination (L_b,ptype_Lamc, p_m, D0_b);
        combination (L_b,ptype_Lamc, p_m, D_m, pi_p);
        combination (L_b,ptype_Lamc, p_m, Dst_m, pi_p);
        combination (L_b,ptype_Lamc, p_m, Dst0_b);
        
        doVertexFit(L_b);
        double recoil_chisq;
        bufchisq=1000000;
        while(L_b.size()>1)
        {
            for(std::vector<Particle>::iterator l = L_b.begin(); l!=L_b.end(); ++l)
            {
                recoil_chisq = dynamic_cast<UserInfo&>(l->userInfo()).vchisq();
                if(recoil_chisq > bufchisq)
                {
                    recoil_chisq = bufchisq;
                    L_b.erase(l); 
                    --l;
                    continue;
                }
                bufchisq = recoil_chisq;
            }
        }
        
        if (!(nevent%1000))std::cout<<nevent<<" recoil candidates: " << L_b.size() << "; chisq/ndf = " << recoil_chisq << '\n';
        cout << "Selecting events" << endl;
        for (std::vector<Particle>::iterator a=L_b.begin(); a!=L_b.end();++a)
        {
            Particle &ALamC=*a;
            
            HepLorentzVector momentum=ALamC.p();
            //       std::cout <<"a1\n";
            // final selection 
          
            
            if ( abs((pUPS-momentum).mag()-2.286)<1.3) 
            {*status=1; skimmed++; return;}
        }
        
        if (!(nevent%1000))std::cout<<nevent<<"     Skimmed: "<<skimmed << '\n';
        
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
            if( (prob_kpr > 0.4) || (prob_pipr > 0.4)) continue; //(prob_kpr > 0.4) || (prob_pipr > 0.4)
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
    
    void setPi0Error(vector<Particle> &p_list){
        for(vector<Particle>::iterator i=p_list.begin(); i!=p_list.end(); ++i)
            setPi0Error(*i);
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
    
    
    HepLorentzVector boostT( HepLorentzVector p, HepLorentzVector p_boost){ //p --boost--> p_boost;                                                                             
        HepLorentzVector tmp(p);
        tmp.boost(-p_boost.boostVector());
        
        return tmp;
    }
    
    
    HepLorentzVector boostT( HepLorentzVector *p, HepLorentzVector *p_boost){ //p --boost--> p_boost;                                                                             
        HepLorentzVector tmp(*p);
        tmp.boost(-p_boost->boostVector());
        return tmp;
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
