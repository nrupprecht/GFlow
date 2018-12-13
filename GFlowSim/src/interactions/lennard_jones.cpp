#include "lennard_jones.hpp"

namespace GFlowSimulation {

  LennardJones::LennardJones(GFlow *gflow) : Interaction(gflow) {
    // Set up parameters
    parameters = new RealType[2];
    parameters[0] = DEFAULT_LENNARD_JONES_STRENGTH;
    parameters[1] = DEFAULT_LENNARD_JONES_CUTOFF;
    // Set the force function
    simd_kernelPtr   = &force<simd_float>;
    serial_kernelPtr = &force_serial;
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

  void LennardJones::force_serial(RealType *buffer_out, const RealType* normal, const RealType mask, const RealType distance, 
    const RealType *soa_data, const RealType *vec_data, const RealType *param_pack, RealType *data_pack, int dimensions) 
  {

    //if (mask!=1) return;

    // Param Pack: Expect:
    const RealType cutoff = 2.5; // param_pack[0];

    // Expect:
    //  soa_data[0] - sigma, 1
    //  soa_data[1] - sigma, 2
    const RealType sg1  = soa_data[0];
    const RealType sg2  = soa_data[1];

    RealType inv_dist = 1./distance;


    RealType gamma = (sg1+sg2)/cutoff*inv_dist;
    RealType g3 = gamma*gamma*gamma; 
    RealType g6 = g3*g3;
    RealType g12 = g6*g6;
   
    // Calculate magnitude
    RealType magnitude = 24.*DEFAULT_LENNARD_JONES_STRENGTH*(2.*g12 - g6)*inv_dist;
  
    // Out:
    //  buffer_out[0::DIM]       = force, 1
    //  buffer_out[DIM+1::2*DIM] = force, 2
    scalarMultVec(magnitude, normal, buffer_out, dimensions); // F1 += f
    copy_negative(buffer_out, &buffer_out[dimensions], dimensions);
  }
  
}
