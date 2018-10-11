#include "lennard_jones.hpp"

namespace GFlowSimulation {

  LennardJones::LennardJones(GFlow *gflow) : Interaction(gflow) {
    // Set up parameters
    parameters = new RealType[2];
    parameters[0] = DEFAULT_LENNARD_JONES_STRENGTH;
    parameters[1] = DEFAULT_LENNARD_JONES_CUTOFF;
    // Set the force function
    simd_kernelPtr = &force<simd_float>;
    serial_kernelPtr = &force<float>;
    // Set data needed - sigmas
    data_needed.push_back(0); // Address of sigma array
  };

  void LennardJones::initialize() {
    parameters[0] = DEFAULT_LENNARD_JONES_STRENGTH;
    parameters[1] = DEFAULT_LENNARD_JONES_CUTOFF;
  }

  void LennardJones::setStrength(RealType s) {
    if (s>=0) parameters[0] = s;
  }

  void LennardJones::setCutoff(RealType cut) {
    if (cut>1.) parameters[1] = cut;
  }
  
}
