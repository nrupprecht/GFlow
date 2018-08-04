#ifndef __LENNARD_JONES_A_HPP__GFLOW__
#define __LENNARD_JONES_A_HPP__GFLOW__

#include "force.hpp"

namespace GFlowSimulation {

  /*
  *  @class LennardJones where all particles have the same force strength and 
  *  radius.
  *
  *  Lennard Jones force. The particle sigma will represent the force cutoff,
  *  generally 2.5*sig, where sig is the inter-particle distance where V=0.
  *  We use cutoff=2.5 by default, but it can be changed. Strength is the
  *  "epsilon" parameter in LJ.
  *
  *  In addition, all particles are assumed to have the same radius, [sigma].
  *
  *  It turns out, this does not seem to be faster than the regular Lennard 
  *  Jones force class.
  *
  *  @see LennardJones
  */
  class LennardJonesA : public Force {
  public:
    //! Constructor
    LennardJonesA(GFlow *);

    //! Calculate all the forces between atoms in the verlet lists
    virtual void calculateForces() const final;

    //! Set the force strength.
    void setStrength(RealType);

    //! Set the particle cutoff
    void setSigma(RealType);

  private:
    //! Calculate force strength
    void forceStrength(RealType*, const RealType*, const RealType, const int, const int) const;

    //! The assumed cutoff radius of the particles.
    RealType sigma;

    //! LJ strength.
    RealType strength;
    //! LJ cutoff - this is 2.5 by default.
    RealType cutoff;
  };

}
#endif // __LENNARD_JONES_A_HPP__GFLOW__