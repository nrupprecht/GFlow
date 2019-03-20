#include "verletlist-pairs.hpp"
// Other files
#include "../base/interaction.hpp"

namespace GFlowSimulation {

  VerletListPairs::VerletListPairs(GFlow *gflow) : InteractionHandler(gflow) {};

  void VerletListPairs::addPair(const int id1, const int id2) {
    verlet.push_back(id1);
    verlet.push_back(id2);
  }

  void VerletListPairs::clear() {
    verlet.clear();
  }

  int VerletListPairs::size() const {
    return verlet.size();
  }

  void VerletListPairs::execute(const Kernel kernel, RealType *param_pack) const {
    // If the kernel is null, then there is no point looping through everything
    if (kernel==nullptr) return;

    // Get the positions
    RealType **x = Base::simData->X();
    RealType **f = Base::simData->F();
    RealType *sg = Base::simData->Sg();
    int *type = Base::simData->Type();

    RealType *displacement = new RealType[sim_dimensions];
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag *boundaryConditions = new BCFlag[sim_dimensions]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions, sim_dimensions); // Keep a local copy of the wrap frags

    RealType sg1, sg2, dx, dy, rsqr, r, invr, magnitude, repulsion = param_pack[0];

    // --- Go through all particles
    for (int i=0; i<verlet.size(); i+=2) {
      int id1 = verlet[i];
      int id2 = verlet[i+1];

      // Check if the types are good
      if (type[id1]<0 || type[id2]<0) continue;

      // Get the displacement
      getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions, sim_dimensions);
      // Get the distance squared
      RealType dsqr = sqr(displacement, sim_dimensions);
      // Check if the particles should interact
      if (dsqr < sqr(sg[id1] + sg[id2])) {
        RealType distance = sqrt(dsqr);
        // Compute the interaction
        kernel(simData, id1, id2, displacement, distance, param_pack, sim_dimensions);
      }
    }

    // Clean up
    delete [] displacement;
    delete [] boundaryConditions;
  }

}
