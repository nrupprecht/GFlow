#ifndef __HARD_SPHERE__GFLOW__
#define __HARD_SPHERE__GFLOW__

#include "force.hpp"

namespace GFlowSimulation {

  /** 
  *  @brief A hard sphere force where all particles have the same repulsion.
  *
  *  The force between particles is proportional to their overlap, (r1 + r2 - distance),
  *  with a constant of proportionality [repulsion] which is a parameter of this class.
  *  In other words, all hard spheres have the same constant of repulsion.
  *
  *  @see HardSphereA
  */
  class HardSphere : public Force {
  public:
    //! Constructor
    HardSphere(GFlow *);
    //! Destructor
    ~HardSphere();

    //! Calculate all the forces between atoms in the verlet lists.
    virtual void calculateForces() const final;

    //! Set the repulsion parameter.
    void setRepulsion(RealType r);

  private:
    //! Calculate force strength
    inline void forceStrength(RealType*, const RealType*, const RealType, const int, const int) const;

    //! The repulsion of the hard spheres. This is assumed to be the same for all hard spheres.
    RealType repulsion;

    // --- TESTS
    RealType *_x1, *_x2, *_disp, *_dist, *_f;
    int BLOCK_SIZE;
  };

}
#endif // __HARD_SPHERE__GFLOW__