#include "genericforce.hpp"

namespace GFlowSimulation {

  template<int NTerms> GenericForce<NTerms>::GenericForce(GFlow *gflow) : Interaction(gflow) {};

}