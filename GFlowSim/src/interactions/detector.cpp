#include "detector.hpp"
// Other files
#include "../interactionhandlers/verletlist-pairs.hpp"

namespace GFlowSimulation {

  Detector::Detector(GFlow *gflow) : Interaction(gflow), ke_threshold(0.001) {};

  Detector::Detector(GFlow *gflow, RealType k) : Interaction(gflow), ke_threshold(k) {};
  
  void Detector::setKEThreshold(RealType k) {
    ke_threshold = k;
  }

}