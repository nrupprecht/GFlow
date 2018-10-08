#include "verletlist.hpp"
// Other files
#include "../utility/memory.hpp"
#include "../utility/simd_utility.hpp"

namespace GFlowSimulation {

  VerletList::VerletList(GFlow *gflow) : InteractionHandler(gflow), lastHead(-1) {};

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

  void VerletList::executeKernel(Kernel kernel, const RealType *param_pack, RealType *data_pack, const vector<int>& data_needed) const 
  {
    // If the kernel is null, then there is no point looping through everything
    if (kernel==nullptr) return;

    // How much extra data is needed
    int data_size = data_needed.size();

    // ****
    data_size = 1;

    // --- Set up some data arrays
    // Position data
    RealType  **Xrt  = alloc_array_2d<RealType>(sim_dimensions, simd_data_size); // SOA format "scratch work" array
    simd_float *Xsd1 = new simd_float[sim_dimensions];
    simd_float *Xsd2 = new simd_float[sim_dimensions];
    simd_float *dX   = new simd_float[sim_dimensions];
    simd_float *norm = new simd_float[sim_dimensions];
    RealType   *temp = new RealType[simd_data_size];
    RealType   *f_temp = new RealType[sim_dimensions];
    RealType   *f_head = new RealType[sim_dimensions];
    // Radii float vectors
    simd_float Sg_sd1, Sg_sd2;
    // Extra data
    simd_float *soa_data = nullptr;
    RealType  **arrays   = nullptr;
    if (data_size) {
      /*
      arrays = new RealType*[data_size];
      // Set pointers in "arrays"
      for (int i=0; i<data_size; ++i) {
        int entry = data_needed[i];
        arrays[i] = Base::simData->dataF[entry]; // @todo This needs to be able to access sigma
      }
      */

      /// **** Doing this by hand right now.
      arrays = new RealType*[1];
      arrays[0] = Base::simData->sg;

      // Two halves - first half will be for head particle (will be duplicated), second half will be for neighbors
      soa_data = new simd_float[2*data_size]; 
    }
    // Buffer for force output
    int buffer_size = sim_dimensions; // @todo Make this an input parameter.
    simd_float *buffer_out = new simd_float[2*buffer_size];

    // Get the positions
    RealType **x = Base::simData->x;
    RealType *sg = Base::simData->sg;

    //cout << "Heads: " << heads.size() << endl;

    for (int i=0; i<heads.size()-1; ++i) {
      // Get next head
      int h0 = heads[i], h1 = heads[i+1];
      int id_h = verlet[h0];
      int j = h0+1;
      // Clear accumulator for head particle's force
      zeroVec(f_temp);
      // Get head particle position data and sigma
      simd_vector_set1(x[id_h], Xsd1, sim_dimensions);
      Sg_sd1 = simd_set1(sg[id_h]);
      // Get other head particle data - pack into first half of soa_data array
      for (int ed=0; ed<data_size; ++ed)
        soa_data[ed] = simd_set1(arrays[ed][id_h]);
      // --- Interact with neighbors
      for (int j=h0+1; j<h1; j+=simd_data_size) {
        int size = min(simd_data_size, h1-j);
        simd_float valid_mask = simd_mask_length( size);

        //cout << "Valid mask: " << valid_mask << endl;

        // Pack positions
        load_vector_data_simd(&verlet[j], x, Xsd2, size, sim_dimensions, Xrt);
        // Pack sigmas
        load_scalar_data_simd(&verlet[j], sg, Sg_sd2, size, temp);
        // Pack other data - into second half of soa_data array
        for (int ed=0; ed<data_size; ++ed)
          load_scalar_data_simd(&verlet[j], arrays[ed], soa_data[data_size+ed], size, temp);

        // --- Check distances
        simd_vector_sub(Xsd1, Xsd2, dX, sim_dimensions);

        /*
        cout << "X1: " << simd_vec_to_str( Xsd1, sim_dimensions) << endl;
        cout << "X2: " << simd_vec_to_str(Xsd2, sim_dimensions) << endl;
        */

        simd_float dX2;
        simd_vector_sqr(dX, dX2, sim_dimensions);
        simd_float cutoff  = simd_add(Sg_sd1, Sg_sd2);
        simd_float cutoff2 = simd_mult(cutoff, cutoff);

        // Check if distance is less than cutoff distance. If so, 0xFFFFFFFF is returned, if not, 0x0 is returned,
        // so we can do a bitwise and of the mask and the force strength. "And" this with valid_mask.
        simd_float mask = simd_mask(valid_mask, simd_less_than(dX2, cutoff2));

        /*
        cout << "Mask: " << mask << endl;
        cout << "Dx: " << simd_vec_to_str(dX, sim_dimensions) << endl;
        cout << "Dsqr: " << dX2 << endl;
        */

        simd_float distance = simd_sqrt(dX2);
        simd_float invDistance = simd_inv_sqrt(dX2);

        // Get normal vectors
        simd_scalar_mult_vec(invDistance, dX, norm, sim_dimensions);
        kernel(buffer_out, norm, mask, distance, soa_data, nullptr, nullptr);

        // Update forces for other particles
        // update_vector_data_size(&verlet[j], Base::simData->f, &buffer_out[data_size], size, sim_dimensions);

        // Update head particle force buffer
        // simd_consolidate_update(f_head, buffer_out, data_size);

      }
      // Collect forces for head particle and write them to the force buffer      
      plusEqVec(Base::simData->f[id_h], f_temp);
    }

    // --- Clean up arrays
    dealloc_array_2d(Xrt);
    delete [] Xsd1;
    delete [] Xsd2;
    delete [] dX;
    delete [] buffer_out;
    delete [] f_temp;
    delete [] f_head;
    if (soa_data) delete [] soa_data;
    if (temp) delete [] temp;

    ///*********


    /*
    // Get the data we need
    int id1(0), id2(0); // List length, id pointers
    RealType **x = Base::simData->x, **f = Base::simData->f, *sg = Base::simData->sg;
    RealType displacement[DIMENSIONS]; // To calculate displacement, normal vector
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags

    // --- Go through all particles
    for (int i=0; i<verlet.size(); i+=2) {
      id1 = verlet[i];
      id2 = verlet[i+1];
      // Type mask
      bool mask = (simData->type[id1]>-1 || simData->type[id2]>-1);
      // Get the displacement between the particles
      getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions);
      // Mast the distance squared with the "particles are real" type mask, c1
      RealType dsqr = sqr(displacement);
      // Check if the particles should interact
      if (mask && dsqr < sqr(sg[id1] + sg[id2])) {
        RealType distance = sqrt(dsqr);
        scalarMultVec(1./distance, displacement);
        // Calculate force strength. Normal will hold the force strength after the function is called.
        kernel(displacement, distance, id1, id2, Base::simData, param_pack, data_pack);
      }
    }
    */
  }

}
