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
   if( iti !=  ehimgr.end() && *iti )
     {
       fox = (*iti).R2();
       //       ntrk = (*iti).Ntrk();
     }



   //----------------------------------------------------------

   std::vector<Particle> lam, lamb; 
   makeLam(lam,lamb);
   setUserInfo(lam,  11 ); 
   setUserInfo(lamb, 12 ); 
   //   doMassVertexFit(lam);
   //   doMassVertexFit(lamb);
   
   if (lam.size()+lamb.size()==0)
     return;


   std::vector<Particle>  e_p,e_m,mu_p,mu_m;
   makeLepton(e_p,e_m,mu_p,mu_m);
   withMuId(mu_p);
   withMuId(mu_m);

   withEId(e_p);
   withEId(e_m);

   std::vector<Particle>  k_p, k_m, pi_p, pi_m,pions;
   makeKPi(k_p, k_m, pi_p, pi_m,1);
   for(std::vector<Particle>::iterator l = pi_m.begin(); l!=pi_m.end(); ++l)
     pions.push_back(*l);
   for(std::vector<Particle>::iterator l = pi_p.begin(); l!=pi_p.end(); ++l)
     pions.push_back(*l);

   ntrk=k_p.size()+k_m.size();
   
   withdRdZcut(k_p,runIp.z());
   withdRdZcut(pi_p,runIp.z());
   withdRdZcut(k_m,runIp.z());
   withdRdZcut(pi_m,runIp.z());

   withKaonId(k_p,0.6,3,1,5);
   withKaonId(k_m,0.6,3,1,5);


   std::vector<Particle>  pi0;
   makePi0(pi0);
   
   for(std::vector<Particle>::iterator i=pi0.begin(); i!=pi0.end();++i)
     if(i->mdstPi0().gamma(0).ecl().energy()<0.05||
	i->mdstPi0().gamma(1).ecl().energy()<0.05||
  	abs(i->mdstPi0().mass()-.135)>0.015)
       {pi0.erase(i); --i;}
   
   //   setPi0Error(pi0);

   //#################################       SIGNAL SIDE  
   std::vector <Particle> Lc, Lcb; 

   combination (Lc,ptype_Lamc, lam, e_p);
   combination (Lcb,ptype_Lamc, lamb, e_m);
   combination (Lc,ptype_Lamc, lam, e_m);     //fake
   combination (Lcb,ptype_Lamc, lamb, e_p);   //fake
   //   setUserInfo(Lc,2);
   //   setUserInfo(Lcb,2);

   combination (Lc,ptype_Lamc, lam, mu_p);
   combination (Lcb,ptype_Lamc, lamb, mu_m);
   combination (Lc,ptype_Lamc, lam, mu_m);     //fake
   combination (Lcb,ptype_Lamc, lamb, mu_p);   //fake
   ///   setUserInfo(Lc,3);
   //   setUserInfo(Lcb,3);


   //   doVertexFit(Lc);
   //   doVertexFit(Lcb);
  
   for(std::vector<Particle>::iterator l = Lc.begin(); l!=Lc.end(); ++l)
     if (l->mass()>ptype_Lamc.mass()+0.1)
       {Lc.erase(l); --l;}
   for(std::vector<Particle>::iterator l = Lcb.begin(); l!=Lcb.end(); ++l)
     if (l->mass()>ptype_Lamc.mass()+0.1)
       {Lcb.erase(l); --l;}
	 
   combination (Lc,ptype_Lamc, lam, pi_p,0.1);
   combination (Lcb,ptype_Lamc, lamb, pi_m,0.1);
   //   setUserInfo(Lc,1);
   //   setUserInfo(Lcb,1);

   combination (Lc,ptype_Lamc, lam, pi_p,pi0,0.1);
   combination (Lcb,ptype_Lamc, lamb, pi_m,pi0,0.1);
   //   setUserInfo(Lc,1);
   //   setUserInfo(Lcb,1);


   if (Lc.size()+Lcb.size()==0)
     return;

   //######################################    TAG SIDE


   std::vector<Particle> p_p, p_m; 
   makeProton(p_p,p_m);
  
   std::vector<Particle> k_s;
   makeKs(k_s);
   for(std::vector<Particle>::iterator l = k_s.begin(); l!=k_s.end(); ++l)
     {
       HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
       Vector3 P(l->px(),l->py(),0);
       V=V-runIp;
       V.setZ(0.);
       if (abs(l->mass()-0.4977)>0.015 || V.perp()<0.1 ||
  	   V.angle(P)>0.01 || l->mdstVee2().z_dist()>10. ) {
  	 k_s.erase(l); --l;
       }
     }
   //   doMassVertexFit(k_s);   
     

   std::vector <Particle> L_, L_b;
   combination (L_,ptype_Lamc,  p_p, k_m);
   combination (L_b,ptype_Lamc, p_m, k_p);
   // setUserInfo(L_,1);
   //setUserInfo(L_b,1);

   combination (L_,ptype_Lamc,  p_p, k_s);
   combination (L_b,ptype_Lamc, p_m, k_s);
   //setUserInfo(L_,2);
   //setUserInfo(L_b,2);

   for (std::vector<Particle>::iterator i=lam.begin(); i!=lam.end();++i)
     L_.push_back(*i);
   for (std::vector<Particle>::iterator i=lamb.begin(); i!=lamb.end();++i)
     L_b.push_back(*i);
    //setUserInfo(L_,3);
    //setUserInfo(L_b,3);

   std::vector<Particle> A; 
   combination(A,ptype_UPS4,Lc,L_b);
   combination(A,ptype_UPS4,Lcb,L_);

   for (std::vector<Particle>::iterator a=A.begin(); a!=A.end();++a)
     {
       Particle &All=*a;
       Particle &LamC=All.child(0);
       Particle &ALamC=All.child(1);
       
       HepLorentzVector momentum=ALamC.p();
       int charge=ALamC.child(0).charge()+ALamC.child(1).charge()+LamC.child(0).charge()+LamC.child(1).charge();
       int n_pi=0;
       int n_pi0=pi0.size();
       for (std::vector<Particle>::iterator pi=pions.begin(); pi!=pions.end();++pi)
	 {
	   //pion cuts
	   if ( checkSame(*a,*pi) ) continue;
	   //
	   n_pi++;
	   charge+=pi->charge();
	   momentum+=pi->p();
	 }
       //       std::cout <<"a1\n";
       // final selection 
       if (charge!=0) continue;
       if (abs((pUPS-(momentum+LamC.p())).mag()<1.5))
	 {*status=1; skimmed++; return;}
       //       double rm =(pUPS-momentum).mag();
       //       if ((pUPS-(momentum+LamC.p())).mag()<2 && (1<rm && rm <3))
       //	 {
	   /*
	   int ch=dynamic_cast<UserInfo&>(ALamC.userInfo()).channel();
	   t1->column("ch",ch);
	   t1->column("npi0",pi0.size());//without pi0
	   t1->column("bestpi0",0);//without pi0                                                                                                                                                                                         
	   t1->column("ml",dynamic_cast<UserInfo&>(LamC.child(0).userInfo()).mass());// lambda mass                                                                                                                                    
	   if (ch>10)
	     t1->column("ml1",dynamic_cast<UserInfo&>(ALamC.userInfo()).mass());// lambda2 mass                                                                                                                                    
	   else
	     t1->column("ml1",0);
	   t1->column("lcch",dynamic_cast<UserInfo&>(LamC.userInfo()).channel());
	   t1->column("rm",rm);
	   t1->column("rm1",(pUPS-(momentum+LamC.p())).mag());
	   t1->column("npi",n_pi);
	   t1->column("mlc",LamC.mass());// lambdac mass                                                                                                                                                                               
	   
	   // lamc heli
	   HepLorentzVector RestP=pStar(ALamC.p(),elec,posi);
	   HepLorentzVector LamcP=HepLorentzVector(-RestP.vect(),RestP.e()); 
	   t1->column("hlc",cos(heli(pStar(LamC.child(0).p(),elec,posi),RestP,LamcP)));
		      	   
	   //lam heli
	   HepLorentzVector p_proton_from_lam; 
	   if (abs(LamC.child(0).child(0).lund())>1000)
	     p_proton_from_lam=LamC.child(0).child(0).p(); 
	   else
	     p_proton_from_lam=LamC.child(0).child(1).p(); 
	   t1->column("hl",cos(heli (p_proton_from_lam, momentum,  LamC.child(0).p())));
	   
	   t1->dumpData();
	   */
       //	   *status=1; skimmed++; return; 
       //	 }

       //       Particle *bestPi0;
       //       double DeltaM=10000;

       for (std::vector<Particle>::iterator pi=pi0.begin(); pi!=pi0.end();++pi)
	 {
	   if ( checkSame(*a,*pi) ) continue;
	   if (abs((pUPS-(momentum+LamC.p()+pi->p())).mag()<1.5))
	     {*status=1; skimmedPi0++; return;}
	 }
     }
   if (!(nevent%1000))std::cout<<nevent<<"     Skimmed: "<<skimmed<<"    SkimmedPi0: "<<skimmedPi0<<'\n';
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
	if(abs(tmp.child(0).lund()>1000)) 
	  prob_kpr=atc_pid(3, 1, 5, 3, 4).prob(tmp.child(0).mdstCharged()); 
	if(abs(tmp.child(1).lund()>1000)) 
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
	if(prob_kpr > 0.9 ) continue;
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

}
