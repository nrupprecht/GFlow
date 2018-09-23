#include "verletlist.hpp"

namespace GFlowSimulation {

  VerletList::VerletList(GFlow *gflow) : InteractionHandler(gflow) {};

  void VerletList::addPair(const int id1, const int id2) {
    verlet.push_back(id1);
    verlet.push_back(id2);
  }

  void VerletList::clear() {
    verlet.clear();
  }

  int VerletList::size() const {
    return verlet.size();
  }

  void VerletList::executeKernel(Kernel kernel, const RealType *param_pack, RealType *data_pack) const 
  {
    // If the kernel is null, then there is no point looping through everything
    if (kernel==nullptr) return;
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
  }

}
