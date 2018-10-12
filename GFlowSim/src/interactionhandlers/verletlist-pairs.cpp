#include "verletlist-pairs.hpp"
// Other files
#include "../utility/memory.hpp"
#include "../utility/simd_utility.hpp"

#include "../interactions/hard_sphere.hpp"
#include "../interactions/lennard_jones.hpp"

namespace GFlowSimulation {

  VerletListPairs::VerletListPairs(GFlow *gflow) : InteractionHandler(gflow), min_simd_size(simd_data_size), use_simd(true) {};

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
    RealType   *temp   = new RealType[simd_data_size];   // A scratchwork buffer

    // Extra data
    simd_float *soa_data = nullptr;
    RealType  **arrays   = nullptr;
    if (data_size) {
      // Create array
      arrays = new RealType*[data_size];
      // Set pointers in "arrays"
      for (int i=0; i<data_size; ++i) {
        int entry = data_needed[i];
        arrays[i] = Base::simData->data_array[i];
      }
      // Two halves - first half will be for head particle (will be duplicated), second half will be for neighbors
      if (use_simd) soa_data = new simd_float[2*data_size]; 
    }
    // Buffer for force output
    int buffer_size = sim_dimensions; // @todo Make this an input parameter.
    simd_float *buffer_out = new simd_float[2*buffer_size];
    float *buffer_out_float = new float[2*buffer_size];
    RealType normal[DIMENSIONS];
    RealType *data = new RealType[2*data_size];

    // Get the positions
    RealType **x = Base::simData->x;
    RealType **f = Base::simData->f;
    RealType *sg = Base::simData->Sg();

    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags

    int i=0, total = verlet_a.size();
    if (use_simd) {
      for (; i<total-simd_data_size; i+=simd_data_size) {
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
        // Pack other data
        for (int ed=0; ed<data_size; ++ed)
          load_scalar_data_simd(&verlet_a[i], arrays[ed], soa_data[ed], size, temp);
        for (int ed=0; ed<data_size; ++ed)
          load_scalar_data_simd(&verlet_b[i], arrays[ed], soa_data[data_size+ed], size, temp);

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
        simd_float distance = simd_sqrt(dX2);
        simd_float invDistance = simd_set1(1.)/dX2;

        // Get normal vectors
        simd_scalar_mult_vec(invDistance, dX, norm, sim_dimensions);

        // Calculate force
        simd_kernel( buffer_out, norm, mask, distance, soa_data, nullptr, param_pack, data_pack);

        // Update forces
        update_vector_data_size(&verlet_a[i], Base::simData->f, &buffer_out[0], size, buffer_size);
        update_vector_data_size(&verlet_b[i], Base::simData->f, &buffer_out[buffer_size], size, buffer_size);
      }
    }
    
    // Seriel part - for left overs, or if we aren't using simd
    
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
        
        // Set data
        for (int dt=0; dt<data_size; ++dt) {
          data[dt] = arrays[dt][id1];
          data[dt+data_size] = arrays[dt][id2];
        }

        // Calculate force.
        serial_kernel( buffer_out_float, normal, 1., distance, data, nullptr, param_pack, data_pack);
        
        // Add the force to the buffers
        plusEqVec(f[id1], &buffer_out_float[0]);
        plusEqVec(f[id2], &buffer_out_float[buffer_size]); 
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
      if (soa_data) delete [] soa_data;
    }
    delete [] buffer_out_float;
    delete [] data;
    delete [] temp;
  }

}
