#include "verletlist.hpp"
// Other files
#include "../utility/memory.hpp"
#include "../utility/simd_utility.hpp"

#include "../interactions/hard_sphere.hpp"
#include "../interactions/lennard_jones.hpp"

namespace GFlowSimulation {

  VerletList::VerletList(GFlow *gflow) : InteractionHandler(gflow), lastHead(-1), min_simd_size(0), use_simd(true) {};

  void VerletList::addPair(const int id1, const int id2) {
    if (id1==lastHead) verlet.push_back(id2);
    else {
      heads.push_back(verlet.size());
      verlet.push_back(id1);
      verlet.push_back(id2);
      lastHead = id1;
    }
  }

  void VerletList::close() {
    // Mark the end of the verlet list with a ficticious head.
    heads.push_back(verlet.size());
  }

  void VerletList::clear() {
    verlet.clear();
    heads.clear();
    lastHead = -1;
  }

  int VerletList::size() const {
    return verlet.size();
  }

  void VerletList::executeKernel(Kernel<simd_float> simd_kernel, Kernel<float> serial_kernel, 
    const RealType *param_pack, RealType *data_pack, const vector<int>& data_needed, const vector<int>& vec_data_needed) const 
  {
    // If the kernel is null, then there is no point looping through everything
    if (simd_kernel==nullptr && serial_kernel==nullptr) return;

    // How much extra data is needed
    int data_size = data_needed.size();
    int vec_data_size = vec_data_needed.size();

    //! @todo - Make sure verlet list works with vector data.
    if (vec_data_size>0) throw false;

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
    RealType   *f_temp = new RealType[sim_dimensions]; // An accumulator for the head particle's force

    // Extra data
    simd_float *soa_data = nullptr;
    RealType  **arrays   = nullptr;
    if (data_size) {
      // Create array
      arrays = new RealType*[data_size];
      // Set pointers in "arrays"
      for (int i=0; i<data_size; ++i) {
        int entry = data_needed[i];
        arrays[i] = Base::simData->data_array[entry];
      }
      // Two halves - first half will be for head particle (will be duplicated), second half will be for neighbors
      if (use_simd) soa_data = new simd_float[2*data_size]; 
    }
    // Buffer for force output
    int buffer_size = sim_dimensions; // @todo Make this an input parameter.
    simd_float *buffer_out = new simd_float[2*buffer_size];

    // Get the positions
    RealType **x = Base::simData->X();
    RealType **f = Base::simData->F();
    RealType *sg = Base::simData->Sg();

    float *buffer_out_float = new float[2*buffer_size];
    RealType normal[DIMENSIONS];
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags
    RealType *data = new RealType[2*data_size];

    for (int i=0; i<heads.size()-1; ++i) {
      // Get next head
      int h0 = heads[i], h1 = heads[i+1];
      int id1 = verlet[h0];
      int j = h0+1;
      
      // Clear accumulator for head particle's force
      zeroVec(f_temp);

      // Potentially use simd functions
      if (use_simd && h1-j>=min_simd_size) {
        // Get head particle position data and sigma
        simd_vector_set1(x[id1], Xsd1, sim_dimensions);
        simd_float Sg_sd1 = simd_set1(sg[id1]);

        // Get other head particle data - pack into first half of soa_data array
        for (int ed=0; ed<data_size; ++ed)
          soa_data[ed] = simd_set1(arrays[ed][id1]);

        // --- Interact with neighbors using simd
        for (; min_simd_size<h1-j; j+=simd_data_size) {
          // Number of valid particles
          int size = min(simd_data_size, h1-j);
          simd_float valid_mask = simd_mask_length(size);

          // --- Pack up data
          // Pack positions
          load_vector_data_simd(&verlet[j], x, Xsd2, size, sim_dimensions, Xrt);
          // Pack sigmas
          simd_float Sg_sd2;
          load_scalar_data_simd(&verlet[j], sg, Sg_sd2, size, temp);
          // Pack other data - into second half of soa_data array
          for (int ed=0; ed<data_size; ++ed)
            load_scalar_data_simd(&verlet[j], arrays[ed], soa_data[data_size+ed], size, temp);

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
          simd_float invDistance = simd_set1(1.)/distance;

          // Get normal vectors
          simd_scalar_mult_vec(invDistance, dX, norm, sim_dimensions);
          simd_kernel(buffer_out, norm, mask, distance, soa_data, nullptr, param_pack, data_pack);

          // Update forces for other particles
          update_vector_data_size(&verlet[j],  Base::simData->F(), &buffer_out[buffer_size], size, buffer_size);

          // Update head particle force buffer
          simd_consolidate_update(f_temp, &buffer_out[0], buffer_size);
        }
      }

      // Seriel part - for left overs, or if we aren't using simd
      for (; j<h1; ++j) {
        int id2 = verlet[j];
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
          serial_kernel(buffer_out_float, normal, 1., distance, data, nullptr, param_pack, data_pack);
         
          // Add the force to the buffers
          plusEqVec(f_temp, &buffer_out_float[0]);
          plusEqVec(f[id2], &buffer_out_float[buffer_size]); 
        }
      }

      // Done with the head - update the head's force.
      plusEqVec(f[id1], f_temp);
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
    delete [] f_temp;
    delete [] buffer_out_float;
    delete [] data;
    delete [] temp;
  }

}
