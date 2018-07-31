#ifndef __GENERIC_FORCE_HPP__GFLOW__
#define __GENERIC_FORCE_HPP__GFLOW__

#include "force.hpp"

namespace GFlowSimulation {

  template<int NTerms> class GenericForce : public Force {
     // Constructor
    GenericForce(GFlow *);

    //! Calculate all the forces between atoms in the verlet lists
    virtual void calculateForces() const final;

  private:
    // Calculate force strength
    inline void forceStrength(RealType*, const RealType*, const RealType, const int, const int) const;

    RealType coefficients[NTerms];
  };

}
#endif // __GENERIC_FORCE_HPP__GFLOW__