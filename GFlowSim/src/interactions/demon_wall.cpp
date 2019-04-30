#include "demon_wall.hpp"
// Other files
#include "../interactionhandlers/verletlist-pairs.hpp"

namespace GFlowSimulation {

  DemonWall::DemonWall(GFlow *gflow) : Interaction(gflow, new VerletListPairs(gflow)) {};

  void DemonWall::interact() const {
    // Do dimensional check.
    // \todo Should probably have some sort of global error message system.
    if (sim_dimensions!=2) return;

    // Get the data pointers.
    RealType **x = Base::simData->X();
    RealType **v = Base::simData->V();
    RealType *sg = Base::simData->Sg();
    RealType *im = Base::simData->Im();
    int    *type = Base::simData->Type();
    int    entry = simData->getIntegerData("Side");
    int    *side = nullptr;
    if (0<=entry) side = Base::simData->IntegerData(entry); // Left: 0, Right: 1

    // Make sure all needed pointers are non null.
    // \todo Should probably have some sort of global error message system.
    if (x==nullptr || v==nullptr || sg==nullptr || im==nullptr || type==nullptr || side==nullptr) return;

    // Get the bounds and boundary conditions
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[2];
    copyVec(Base::gflow->getBCs(), boundaryConditions, 2); // Keep a local copy of the bcs
    // Extract bounds related data
    RealType bnd_x = bounds.wd(0);
    RealType bnd_y = bounds.wd(1);

    // Needed constants
    RealType sg1, sg2, dx, dy, rsqr, r, invr, magnitude;
    // Point to the actual list from the verlet list object. Since we set the handler at initialization to 
    // be of type VerletListPairs, this cast should always succeed.
    vector<int> &verlet = dynamic_cast<VerletListPairs*>(handler)->verlet;

    // --- Go through all particles
    for (int i=0; i<verlet.size(); i+=2) {
      int id1 = verlet[i];
      int id2 = verlet[i+1];
      // Check if the types are good
      if (type[id1]<0 || type[id2]<0) continue;
      // The object of larger type will be the wall, make this id2
      if (type[id1]>type[id2]) std::swap(id1, id2);
      // Calculate displacement.
      dx = x[id1][0] - x[id2][0];
      dy = x[id1][1] - x[id2][1];
      // Harmonic corrections to distance.
      if (boundaryConditions[0]==BCFlag::WRAP) {
        RealType dX = bnd_x - fabs(dx);
        if (dX<fabs(dx)) dx = dx>0 ? -dX : dX;
      }  
      if (boundaryConditions[1]==BCFlag::WRAP) {
        RealType dY = bnd_y - fabs(dy);
        if (dY<fabs(dy)) dy = dy>0 ? -dY : dY;
      }  
      // Calculate squared distance
      rsqr = dx*dx + dy*dy;
      // Get radii
      sg1 = sg[id1];
      sg2 = sg[id2];
      // If close, interact.
      if (0<rsqr && rsqr < sqr(sg1 + sg2)) {

        // Left side particle trying to go right
        if (side[id1]==0 && -0.05<=dx+sg1 && 0<v[id1][0]) {
          v[id1][0] = -v[id1][0];
        }
        // Right side particle trying to go left
        else if (side[id1]==1 && dx-sg1<=0.05 && v[id1][0]<0)
          v[id1][0] = -v[id1][0];
      }
    }
  }

}