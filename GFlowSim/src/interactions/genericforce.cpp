#include "genericforce.hpp"

namespace GFlowSimulation {

  template<int NTerms> GenericForce<NTerms>::GenericForce(GFlow *gflow) : Interaction(gflow) {};

  template<int NTerms> void GenericForce<NTerms>::initialize() {}

}