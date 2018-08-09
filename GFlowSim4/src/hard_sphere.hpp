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

    //! Calculate all the forces between atoms in the verlet lists.
    virtual void calculateForces() const final;

    //! Set the repulsion parameter.
    void setRepulsion(RealType r);

    //! @param[in] normal
    //! @param[in] distance
    //! @param[in] id1
    //! @param[in] id2
    //! @param[in] simData
    //! @param[in] param_pack A parameter pack, passed in from force. Contains characteristic 
    //! constants of the force, and extra data the force needs.
    //! @param[in,out] virial The virial, to be updated by this functiton.
    static void force(RealType*, const RealType, const int, const int, const SimData*, const RealType*, RealType&);

  private:
    //! Calculate force strength
    inline void forceStrength(RealType*, const RealType, const int, const int) const;

    //! The repulsion of the hard spheres. This is assumed to be the same for all hard spheres.
    RealType repulsion;
  };

}
#endif // __HARD_SPHERE__GFLOW__