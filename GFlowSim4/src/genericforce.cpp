#include "genericforce.hpp"

namespace GFlowSimulation {

  template<int NTerms> GenericForce<NTerms>::GenericForce(GFlow *gflow) : Force(gflow) {};

  template<int NTerms> void GenericForce<NTerms>::calculateForces() const {

  }

  template<int NTerms> inline void GenericForce<NTerms>::forceStrength(RealType *F, const RealType 
    *normal, const RealType distance, const int id1, const int id2) const {

  }

}