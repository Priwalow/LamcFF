#include "belle.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "panther/panther.h"
#include <iostream>
#include MDST_H
#include BELLETDF_H
#include HEPEVT_H
#include "particle/utility.h"
#include "particle/combination.h"
#include "kid/atc_pid.h"
#include "mdst/mdst.h"
#include "mdst/Muid_mdst.h"
#include "eid/eid.h"
#include "ip/IpProfile.h"


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


using namespace std;

class User_reco : public Module {
public:
  User_reco ( void );
  ~User_reco ( void ){};
  void init ( int* ){};
  void term ( void ){};
  void disp_stat ( const char* ){};
  void hist_def ( void );
  void event ( BelleEvent*, int* );
  void begin_run ( BelleEvent*, int* );
  void end_run   ( BelleEvent*, int* ){};
  void other ( int*, BelleEvent*, int* ){};

  /*******************************************
    Define things here that you want to use
    inside and outside of the event function
  *******************************************/
  BelleTuple *t1,*t2,*t3,*t4,*t5  ;


private:

  Ptype ptype_D0 ;
  Ptype ptype_D0B ;
  Ptype ptype_Dp ;
  Ptype ptype_Dm ;
  Ptype ptype_Dsp ;
  Ptype ptype_Dsm ;
  Ptype ptype_UPS4 ;
  Ptype ptype_Dstm ;
  Ptype ptype_Dstp ;
  Ptype ptype_Dst0 ;
  Ptype ptype_DstB ;
  Ptype ptype_Lamc ;
  Ptype ptype_ALamc ;
  Ptype ptype_PP;
  Ptype ptype_AP;
  Ptype ptype_Lamc2625; 
  Ptype ptype_ALamc2625;
  Ptype ptype_Lamc2595; 
  Ptype ptype_ALamc2595;
  Ptype ptype_ASigc;
  Ptype ptype_Sigc;
};
extern "C" Module_descr *mdcl_User_reco ()
{ /* main */
  User_reco *module = new User_reco;
  Module_descr *dscr = new Module_descr ( "User_reco", module );
  return dscr;
};

/* Constructor  */
 User_reco::User_reco ( void ):
   ptype_D0("D0") ,
   ptype_D0B("D0B") ,
   ptype_Dp("D+") ,
   ptype_Dm("D-") ,
   ptype_Dsp("Ds+") ,
   ptype_Dsm("Ds-") ,
   ptype_UPS4("UPS4") ,
     ptype_Dstm("D*-"),
     ptype_Dstp("D*+"),
     ptype_Dst0("D*0"),
     ptype_DstB("D*B"),
     ptype_PP("P+"),
     ptype_AP("AP+"),
     ptype_Lamc("LAMC"),
     ptype_ALamc("ALAMC"),
     ptype_Lamc2625("LC2625"),
     ptype_ALamc2625("ALC2625"),
     ptype_Lamc2595("LC2595"),
     ptype_ALamc2595("ALC2595"),
     ptype_ASigc("ASIC+"),  //2.453800
     ptype_Sigc("SIGC+")   // 2.453800
     {
  /*******************************************
   Define things here that you want to use
   inside and outside of the event function
  ******************************************/


};

void User_reco::begin_run ( BelleEvent*, int* )
{
  eid::init_data();
  // Get IP profile data from $BELLE_POSTGRES_SERVER
  IpProfile::begin_run();
  // Dump IP profile data to STDOUT (optional)
  IpProfile::dump();
}

#if defined(BELLE_NAMESPACE)
}
#endif

