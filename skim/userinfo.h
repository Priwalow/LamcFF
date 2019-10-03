#include "belle.h"
#include "particle/Particle.h"
#include "particle/ParticleUserInfo.h"
#include <vector>

// For Interface to Set UserInfo Class

using namespace Belle;
using namespace std;
namespace Belle{

void setUserInfo(Particle &p, unsigned ch=0);
void setUserInfo(std::vector<Particle> &p, unsigned ch=0);

// UserInfo Class

class UserInfo : public ParticleUserInfo
{
public:
  /// Default constructor
  UserInfo();

  UserInfo(unsigned);

  /// Copy constructor
  UserInfo(const UserInfo &);

  /// Destructor
  virtual ~UserInfo();

  /// constructs self object.
  UserInfo * clone(void) const;

  /// Copy operator
  UserInfo & operator = (const UserInfo &);

  
  inline UserInfo * getPointer()
    {return this;}
  
public:
  void vchisq(const double &v) { m_vchisq = v; }
  void p(const HepLorentzVector &v) { m_p = v; }
  void mchisq(const double &v) { m_mchisq = v; }
  void mass(const double &v) { m_mass = v; } 

  void vx(const double &v) { m_vx = v; }
  void vy(const double &v) { m_vy = v; } 
  void vz(const double &v) { m_vz = v; }
  void id(const double &v) { m_id = v; } 
  void K_angl(const double &v) { m_K_angl = v; }
  void Pi_angl(const double &v) { m_Pi_angl = v; }
  const double & vchisq(void) const { return m_vchisq; }
  const HepLorentzVector & p(void) const { return m_p; }
  const double & mchisq(void) const { return m_mchisq; }
  const double & mass(void) const { return m_mass; }
  const double & vx(void) const { return m_vx; }
  const double & vy(void) const { return m_vy; }
  const double & vz(void) const { return m_vz; }
  const double & id(void) const { return m_id; }
  const double & K_angl(void) const { return m_K_angl; }
  const double & Pi_angl(void) const { return m_Pi_angl; }
  
  void vmass(const double &v) { m_vmass = v; }
  const double & vmass(void) const { return m_vmass; }

  void channel(const unsigned &v) { m_channel = v; }
  const unsigned & channel(void) const { return m_channel; }

private:
  double   m_vchisq;
  double   m_mchisq;
  double   m_vmass;
  double   m_mass;

  double   m_vx;
  double   m_vy;
  double   m_vz;  
  double   m_id;
  double   m_K_angl;
  double   m_Pi_angl;
  unsigned m_channel;
  HepLorentzVector m_p;
};
}
