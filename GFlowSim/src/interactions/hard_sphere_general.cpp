#include "hard_sphere_general.hpp"

namespace GFlowSimulation {

  HardSphereGeneral::HardSphereGeneral(GFlow *gflow) : Interaction(gflow) {
    parameters = new RealType[2];
    parameters[0] = DEFAULT_HARD_SPHERE_REPULSION;
    parameters[1] = DEFAULT_HARD_SPHERE_DISSIPATION; // Dissipation
    // Set the force function
    simd_kernelPtr = &force<simd_float>;
    serial_kernelPtr = &force<float>;
    // Set data needed - sigmas
    data_needed.push_back(0); // Address of sigma array
    // Set vector data needed - velocities
    vec_data_needed.push_back(1);
  };

  void HardSphereGeneral::initialize() {
    parameters[0] = DEFAULT_HARD_SPHERE_REPULSION;
    parameters[1] = DEFAULT_HARD_SPHERE_DISSIPATION; // Dissipation
  }

  void HardSphereGeneral::setRepulsion(RealType r) { 
    parameters[0] = r; 
  }

}
