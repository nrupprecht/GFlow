#ifndef __HARD_SPHERE_A_HPP__GFLOW__
#define __HARD_SPHERE_A_HPP__GFLOW__

#include "force.hpp"

namespace GFlowSimulation {

  /** 
  *  @brief A hard sphere force where all particles have the same repulsion and radius.
  *
  *  The force between particles is proportional to their overlap, (r1 + r2 - distance),
  *  with a constant of proportionality [repulsion] which is a parameter of this class.
  *  In other words, all hard spheres have the same constant of repulsion.
  *
  *  In addition, all particles are assumed to have the same radius, [sigma].
  *
  *  It turns out, this does not seem to be faster than the regular Hard Sphere force class.
  *
  *  @see HardSphere
  */
  class HardSphereA : public Force {
    public:
    //! Constructor
    HardSphereA(GFlow *);

    //! Calculate all the forces between atoms in the verlet lists.
    virtual void calculateForces() const final;

    //! Set the repulsion parameter.
    void setRepulsion(RealType);

    //! Set the hard sphere radius parameter.
    void setSigma(RealType);

  private:
    //! Calculate force strength
    inline void forceStrength(RealType*, const RealType*, const RealType, const int, const int) const;

    //! The (assumed) particle radius.
    RealType sigma;

    //! The repulsion of the hard spheres. This is assumed to be the same for all hard spheres.
    RealType repulsion;
  };

}
#endif // __HARD_SPHERE_A_HPP__GFLOW__