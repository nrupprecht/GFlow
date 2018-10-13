#include "verletlist-pairs.hpp"
// Other files
#include "../utility/memory.hpp"
#include "../utility/simd_utility.hpp"

#include "../interactions/hard_sphere.hpp"
#include "../interactions/lennard_jones.hpp"

namespace GFlowSimulation {

  VerletListPairs::VerletListPairs(GFlow *gflow) : InteractionHandler(gflow), min_simd_size(0), use_simd(true) {};

  void VerletListPairs::addPair(const int id1, const int id2) {
    verlet_a.push_back(id1);
    verlet_b.push_back(id2);
  }

  void VerletListPairs::clear() {
    verlet_a.clear();
    verlet_b.clear();
  }

  int VerletListPairs::size() const {
    return 2*verlet_a.size();
  }

  void VerletListPairs::executeKernel(Kernel<simd_float> simd_kernel, Kernel<float> serial_kernel, 
    const RealType *param_pack, RealType *data_pack, 
    const vector<int>& data_needed, const vector<int>& vec_data_needed) const 
  {
    // If the kernel is null, then there is no point looping through everything
    if (simd_kernel==nullptr && serial_kernel==nullptr) return;

    // How much extra data is needed
    int data_size = data_needed.size();
    int vec_data_size = vec_data_needed.size();

    // --- Set up some data arrays
    RealType  **Xrt  = nullptr; // SOA format "scratch work" array
    simd_float *Xsd1 = nullptr; // Simd vector for position
    simd_float *Xsd2 = nullptr; // Simd vector for position
    simd_float *dX   = nullptr; // Simd vector for displacement
    simd_float *norm = nullptr; // Simd vector for normal force
    if (use_simd) { // Only allocate these if neccessary
      Xrt  = alloc_array_2d<RealType>(sim_dimensions, simd_data_size); // SOA format "scratch work" array
      Xsd1 = new simd_float[sim_dimensions]; // Simd vector for position
      Xsd2 = new simd_float[sim_dimensions]; // Simd vector for position
      dX   = new simd_float[sim_dimensions]; // Simd vector for displacement
      norm = new simd_float[sim_dimensions]; // Simd vector for normal force
    }
    RealType *temp   = new RealType[simd_data_size];   // A scratchwork buffer

    // Buffer for force output
    int buffer_size = sim_dimensions; // @todo Make this an input parameter.
    simd_float *buffer_out = new simd_float[2*buffer_size];
    float      *serial_buffer_out = new float[2*buffer_size];
    simd_float *scalar_data = nullptr;
    RealType   *serial_scalar_data = nullptr;
    simd_float *vec_data = nullptr;
    RealType   *serial_vec_data = nullptr;
    RealType   *normal = new RealType[DIMENSIONS];

    // Arrays for getting data
    RealType  **arrays   = nullptr;
    RealType ***vec_arrays = nullptr;
    if (data_size) {
      // Create array
      arrays = new RealType*[data_size];
      // Set pointers in "arrays"
      for (int i=0; i<data_size; ++i) {
        int entry = data_needed[i];
        arrays[i] = Base::simData->data_array[entry];
      }
      // Allocate arrays
      if (use_simd) scalar_data = new simd_float[2*data_size];
      serial_scalar_data = new RealType[2*data_size];
    }
    if (vec_data_size) {
      // Create array
      vec_arrays = new RealType**[vec_data_size];
      // Set pointers in "vec_arrays"
      for (int i=0; i<vec_data_size; ++i) {
        int entry = vec_data_needed[i];
        vec_arrays[i] = Base::simData->vec_data_array[entry];
      }
      // Allocate arrays
      if (use_simd) vec_data = new simd_float[2*vec_data_size*sim_dimensions];
      serial_vec_data = new RealType[2*vec_data_size*sim_dimensions];
    }
    
    // Get the positions
    RealType **x = Base::simData->X();
    RealType **f = Base::simData->F();
    RealType *sg = Base::simData->Sg();

    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags

    RealType TV[8], TV2[8];

    int i=0, total = verlet_a.size();
    if (use_simd) {
      for (; i<total-min_simd_size; i+=simd_data_size) {
        // Number of valid particles
        int size = min(simd_data_size, static_cast<int>(verlet_a.size())-i);
        simd_float valid_mask = simd_mask_length(size);

        // --- Pack up data
        // Pack positions
        load_vector_data_simd(&verlet_a[i], x, Xsd1, size, sim_dimensions, Xrt);
        load_vector_data_simd(&verlet_b[i], x, Xsd2, size, sim_dimensions, Xrt);
        
        // Pack sigmas
        simd_float Sg_sd1, Sg_sd2;
        load_scalar_data_simd(&verlet_a[i], sg, Sg_sd1, size, temp);
        load_scalar_data_simd(&verlet_b[i], sg, Sg_sd2, size, temp);
        
        // Pack scalar data
        for (int ed=0; ed<data_size; ++ed)
          load_scalar_data_simd(&verlet_a[i], arrays[ed], scalar_data[ed], size, temp);
        for (int ed=0; ed<data_size; ++ed)
          load_scalar_data_simd(&verlet_b[i], arrays[ed], scalar_data[data_size+ed], size, temp);
        // Pack vector data
        for (int ed=0; ed<vec_data_size; ++ed)
          load_vector_data_simd(&verlet_a[i], vec_arrays[ed], &vec_data[ed*sim_dimensions], size, sim_dimensions, Xrt);
        for (int ed=0; ed<vec_data_size; ++ed)
          load_vector_data_simd(&verlet_a[i], vec_arrays[ed], &vec_data[(vec_data_size+ed)*sim_dimensions], size, sim_dimensions, Xrt);
        

        // --- Check distances
        simd_get_displacement(Xsd1, Xsd2, dX, bounds, boundaryConditions, sim_dimensions);

        // ---Get distance squared, calculate cutoff radius
        simd_float dX2;
        simd_vector_sqr(dX, dX2, sim_dimensions);
        simd_float cutoff  = simd_add(Sg_sd1, Sg_sd2);
        simd_float cutoff2 = simd_mult(cutoff, cutoff);

        // Check if distance is less than cutoff distance. If so, 0xFFFFFFFF is returned, if not, 0x0 is returned,
        // so we can do a bitwise and of the mask and the force strength. "And" this with valid_mask.
        simd_float mask = simd_mask(valid_mask, simd_less_than(dX2, cutoff2));

        // --- Calculate distance and inverse distance
        simd_float distance    = simd_sqrt(dX2);
        simd_float invDistance = simd_set1(1.)/distance;

        // Get normal vectors
        simd_scalar_mult_vec(invDistance, dX, norm, sim_dimensions);

        // Calculate force
        simd_kernel(buffer_out, norm, mask, distance, scalar_data, vec_data, param_pack, data_pack);

        // Update forces
        update_vector_data_size(&verlet_a[i], Base::simData->F(), &buffer_out[0], size, buffer_size);
        update_vector_data_size(&verlet_b[i], Base::simData->F(), &buffer_out[buffer_size], size, buffer_size);
      }
    }
    
    // Serial part - for left overs, or if we aren't using simd
    for (; i<total; ++i) {
      int id1 = verlet_a[i];
      int id2 = verlet_b[i];
      // Get the displacement between the particles (stored in "normal")
      getDisplacement(x[id1], x[id2], normal, bounds, boundaryConditions);
      // Mast the distance squared with the "particles are real" type mask, c1
      RealType dsqr = sqr(normal);
      // Check if the particles should interact
      if (dsqr < sqr(sg[id1] + sg[id2])) {
        RealType distance = sqrt(dsqr);
        scalarMultVec(1./distance, normal);
        
        // Set scalar data
        for (int dt=0; dt<data_size; ++dt)
          serial_scalar_data[dt] = arrays[dt][id1];
        for (int dt=0; dt<data_size; ++dt)
          serial_scalar_data[dt+data_size] = arrays[dt][id2];
        // Set vector data
        for (int dt=0; dt<vec_data_size; ++dt) {
          copyVec(vec_arrays[dt][id1], &serial_vec_data[sim_dimensions*dt]);
          copyVec(vec_arrays[dt][id2], &serial_vec_data[sim_dimensions*(dt+vec_data_size)]);
        }

        // Calculate force.
        serial_kernel(serial_buffer_out, normal, 1., distance, serial_scalar_data, serial_vec_data, param_pack, data_pack);
        
        // Add the force to the buffers
        // plusEqVec(f[id1], &serial_buffer_out[0]);
        // plusEqVec(f[id2], &serial_buffer_out[buffer_size]); 
      }
    }

    // --- Clean up arrays
    if (use_simd) {
      dealloc_array_2d(Xrt);
      delete [] Xsd1;
      delete [] Xsd2;
      delete [] dX;
      delete [] buffer_out;
      delete [] norm;
      if (scalar_data) delete [] scalar_data;
      if (vec_data)    delete [] vec_data;
    }
    delete [] serial_buffer_out;
    delete [] serial_scalar_data;
    delete [] temp;
    delete [] normal;
  }

}
