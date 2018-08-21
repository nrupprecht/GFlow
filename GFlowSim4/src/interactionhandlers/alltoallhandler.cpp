#include "alltoallhandler.hpp" 

namespace GFlowSimulation {

  AllToAllHandler::AllToAllHandler(GFlow *gflow) : InteractionHandler(gflow) {};

  int AllToAllHandler::size() const {
    return Base::simData->number * (Base::simData->number - 1) / 2;
  }

  void AllToAllHandler::executeKernel(Kernel kernel, const RealType *param_pack, RealType *data_pack) const {
    int number = Base::simData->number;
    RealType **x = Base::simData->x, **f = Base::simData->f, *sg = Base::simData->sg;

    RealType displacement[DIMENSIONS], normal[DIMENSIONS]; // To calculate displacement, normal vector
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags

    // Loop through all pairs of particles
    for (int id1=0; id1<number; ++id1) {
      // Check that id1's type is valid
      if (simData->type[id1]<0) continue;
      // Loop through all particles of id > id1
      for (int id2=id1+1; id2<number; ++id2) {
        // Type mask
        if (simData->type[id2]<0) continue;
        // Get the displacement between the particles
        getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions);
        // Mast the distance squared with the "particles are real" type mask, c1
        RealType dsqr = sqr(displacement);
        // Check if the particles should interact
        if (dsqr < sqr(sg[id1] + sg[id2])) {
          RealType distance = sqrt(dsqr);
          scalarMultVec(1./distance, displacement, normal);
          // Calculate force strength. Normal will hold the force strength after the function is called.
          kernel(normal, distance, id1, id2, Base::simData, param_pack, data_pack);
        }
      }
    }
  }

}