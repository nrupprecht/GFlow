#include "hard_sphere.hpp"

namespace GFlowSimulation {

  HardSphere::HardSphere(GFlow *gflow) : Interaction(gflow) {
    parameters = new RealType;
    parameters[0] = DEFAULT_HARD_SPHERE_REPULSION;
    // Set the force function
    simd_kernelPtr = &force<simd_float>;
    serial_kernelPtr = &force<float>;
    // Set data needed - sigmas
    data_needed.push_back(0); // Address of sigma array
  };

  void HardSphere::initialize() {
    parameters[0] = DEFAULT_HARD_SPHERE_REPULSION;
  }

  void HardSphere::setRepulsion(RealType r) { 
    parameters[0] = r; 
  }

}
