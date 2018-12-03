#include "hard_sphere.hpp"

namespace GFlowSimulation {

  HardSphere::HardSphere(GFlow *gflow) : Interaction(gflow) {
    parameters = new RealType;
    parameters[0] = DEFAULT_HARD_SPHERE_REPULSION;
    // Set the force function
    simd_kernelPtr = &force<simd_float>;
    serial_kernelPtr = &force2d;//&force<float>;
    // Set data needed - sigmas
    data_needed.push_back(0); // Address of sigma array
  };

  void HardSphere::initialize() {
    parameters[0] = DEFAULT_HARD_SPHERE_REPULSION;
  }

  void HardSphere::setRepulsion(RealType r) { 
    parameters[0] = r; 
  }

  void HardSphere::force2d(
    RealType*       buffer_out, 
    const RealType* normal, 
    const RealType  mask, 
    const RealType  distance, 
    const RealType* soa_data, 
    const RealType* vec_data, 
    const RealType* param_pack, 
    RealType*       data_pack, 
    int             dimensions) 
  {

    if (mask!=1) return;

    // Expect: Soa data:
    //  soa_data[0] - sigma, 1
    //  soa_data[1] - sigma, 2
    const RealType sg1  = soa_data[0];
    const RealType sg2  = soa_data[1];
    // Expect: Param pack:
    //  param_pack[0] - Repulsion
    const RealType repulsion = param_pack[0];

    // Calculate magnitude
    RealType magnitude = repulsion*(sg1 + sg2 - distance);

    buffer_out[0] = magnitude*normal[0];
    buffer_out[2] = -buffer_out[0];
    buffer_out[1] = magnitude*normal[1];
    buffer_out[3] = -buffer_out[1];
  }

}
