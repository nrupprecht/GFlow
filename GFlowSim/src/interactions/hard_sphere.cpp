#include "hard_sphere.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../utility/simd_utility.hpp"

namespace GFlowSimulation {

  HardSphere::HardSphere(GFlow *gflow) : Interaction(gflow) {
    parameters = new RealType;
    *parameters = DEFAULT_HARD_SPHERE_REPULSION;
    // Set the force function
    kernelPtr = &force;
  };

  void HardSphere::initialize() {
    *parameters = DEFAULT_HARD_SPHERE_REPULSION;
  }

  void HardSphere::setRepulsion(RealType r) { 
    parameters[0] = r; 
  }

  //! @param[out] buffer_out For forces, torques, etc.
  //! @param[in] normal
  //! @param[in] distance
  //! @param[in] soa_data
  //! @param[in] param_pack A parameter pack, passed in from force. Should be of the form { repulsion } (length 1).
  //! @param[in,out] data_pack Data to update in the function. Should be of the form  { virial } (length 1). 
  //! Add the f_i \dot r_i to this.
  void HardSphere::force(simd_float *buffer_out, simd_float* normal, const simd_float mask, const simd_float distance, const simd_float *soa_data, 
    const RealType *param_pack, RealType *data_pack) 
  {
    // Expect:
    //  soa_data[0] - sigma, 1
    //  soa_data[1] - repulsion, 1
    //  soa_data[2] - sigma, 2
    //  soa_data[3] - repulsion, 2

    const simd_float &sg1  = soa_data[0];
    const simd_float &sg2  = soa_data[1]; // 2
    /*
    const simd_float &rep1 = soa_data[2];
    const simd_float &rep2 = soa_data[3];
    */
    const simd_float rep1 = simd_set1(DEFAULT_HARD_SPHERE_REPULSION);
    const simd_float rep2 = rep1;

    simd_float acc1 = simd_add(rep1, rep2);
    simd_float acc2 = simd_add(sg1, sg2);
    simd_float acc3 = simd_sub(acc2, distance);

    //cout << "sg + sg - dist: " << simd_to_str(acc3) << endl;

    const simd_float repl = simd_mult(rep1, rep2);

    simd_float magnitude = simd_mult(repl, acc3);

    // Apply force mask
    simd_float masked_magnitude = simd_mask(magnitude, mask);

    /*
    cout << "Masked magnitude: " << simd_to_str(masked_magnitude) << endl;
    cout << "Normal vectors: " << simd_vec_to_str(normal, DIMENSIONS) << endl;
    */

    // Out:
    //  buffer_out[0::DIM]       = force, 1
    //  buffer_out[DIM+1::2*DIM] = force, 2
    simd_scalar_mult_vec( masked_magnitude, normal, &buffer_out[0], DIMENSIONS); // F1 += f
    simd_scalar_mult_vec( minus_one, &buffer_out[0], &buffer_out[DIMENSIONS], DIMENSIONS); // F2 -= f

    //cout << "F1 =  " << simd_vec_to_str( &buffer_out[0] , DIMENSIONS) << endl;
    //cout << "F2 =  " << simd_vec_to_str( &buffer_out[DIMENSIONS] , DIMENSIONS) << endl;

    /*
    // Calculate force strength
    RealType magnitude = param_pack[0]*(simData->sg[id1] + simData->sg[id2] - distance);
    // Update the virial
    data_pack[0] += magnitude;
    // Force strength x Normal vector -> Sets normal to be the vectorial force between the particles.
    scalarMultVec(magnitude, normal);
    // Add the force to the buffers
    plusEqVec (simData->f[id1], normal);
    minusEqVec(simData->f[id2], normal);
    */
  }

}
