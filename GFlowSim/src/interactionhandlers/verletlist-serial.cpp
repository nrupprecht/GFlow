#include "verletlist-serial.hpp"
// Other files
#include "../utility/memory.hpp"
#include "../utility/simd_utility.hpp"

#include "../interactions/hard_sphere.hpp"
#include "../interactions/lennard_jones.hpp"

namespace GFlowSimulation {

  VerletListSerial::VerletListSerial(GFlow *gflow) : InteractionHandler(gflow) {};

  void VerletListSerial::addPair(const int id1, const int id2) {
    verlet.push_back(id1);
    verlet.push_back(id2);
  }

  void VerletListSerial::clear() {
    verlet.clear();
  }

  int VerletListSerial::size() const {
    return verlet.size();
  }

  void VerletListSerial::executeKernel(Kernel<simd_float> simd_kernel, Kernel<float> serial_kernel, 
    const RealType *param_pack, RealType *data_pack, const vector<int>& data_needed) const 
  {
    // If the kernel is null, then there is no point looping through everything
    if (serial_kernel==nullptr) return;

    // How much extra data is needed
    int data_size = data_needed.size();

    RealType  **arrays   = nullptr;
    if (data_size) {
      // Create array
      arrays = new RealType*[data_size];
      // Set pointers in "arrays"
      for (int i=0; i<data_size; ++i) {
        int entry = data_needed[i];
        //arrays[i] = Base::simData->dataF[entry]; // @todo This needs to be able to access sigma
        if (i==0) arrays[i] = Base::simData->sg;
      }
      // Two halves - first half will be for head particle (will be duplicated), second half will be for neighbors
    }
    // Buffer for force output
    int buffer_size = sim_dimensions; // @todo Make this an input parameter.
    float *buffer_out_float = new float[2*buffer_size];
    RealType normal[DIMENSIONS];
    RealType *data = new RealType[2*data_size];

    // Get the positions
    RealType **x = Base::simData->x;
    RealType **f = Base::simData->f;
    RealType *sg = Base::simData->sg;

    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags

    // --- Go through all particles
    for (int i=0; i<verlet.size(); i+=2) {
      int id1 = verlet[i];
      int id2 = verlet[i+1];
      // Get the displacement between the particles
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
        serial_kernel( buffer_out_float, normal, 1., distance, data, param_pack, data_pack);
        
        // Add the force to the buffers
        plusEqVec(f[id1], &buffer_out_float[0]);
        plusEqVec(f[id2], &buffer_out_float[buffer_size]); 
      }
    }


    if (arrays) delete [] arrays;
  }

}
