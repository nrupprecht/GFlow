#include "verletlist.hpp"
// Other files
#include "../base/interaction.hpp"

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

  void VerletList::execute(const Interaction *interaction) const {
    // If the kernel is null, then there is no point looping through everything
    if (interaction==nullptr) return;

    // Get the positions
    RealType **x = Base::simData->X();
    RealType *sg = Base::simData->Sg();

    RealType *displacement = new RealType[sim_dimensions];
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag *boundaryConditions = new BCFlag[sim_dimensions]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions, sim_dimensions); // Keep a local copy of the wrap frags

    // Go through verlet list
    for (int i=0; i<heads.size()-1; ++i) {
      // Get next head
      int h0 = heads[i], h1 = heads[i+1];
      int id1 = verlet[h0];
      int j = h0+1;
      
      // Seriel part - for left overs, or if we aren't using simd
      for (; j<h1; ++j) {
        int id2 = verlet[j];
        getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions, sim_dimensions);
        // Mast the distance squared with the "particles are real" type mask, c1
        RealType dsqr = sqr(displacement, sim_dimensions);
        // Check if the particles should interact
        if (dsqr < sqr(sg[id1] + sg[id2])) {
          RealType distance = sqrt(dsqr);
          // Compute the interaction
          interaction->compute(id1, id2, displacement, distance);
        }
      }
    }
  }

}
