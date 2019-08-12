#include "interaction__3d.hpp"

namespace GFlowSimulation {

  Interaction3d::Interaction3d(GFlow *gflow) : Interaction(gflow) {};

  void Interaction3d::interact() const {
    // Common tasks.
    Interaction::interact();

    // Do dimensional check.
    // \todo Should probably have some sort of global error message system.
    if (sim_dimensions!=3) return;

    // Get the data pointers.
    RealType **x = Base::simData->X();
    RealType **f = Base::simData->F();
    RealType *sg = Base::simData->Sg();
    int    *type = Base::simData->Type();

    // Make sure all needed pointers are non null.
    // \todo Should probably have some sort of global error message system.
    if (x==nullptr || f==nullptr || sg==nullptr || type==nullptr) return;

    // Needed constants
    RealType sg1, sg2, dr[3], rsqr;

    // --- Go through all particles in verlet.
    for (int i=0; i<verlet.size(); i+=2) {
      int id1 = verlet[i];
      int id2 = verlet[i+1];
      // Check if the types are good
      if (type[id1]<0 || type[id2]<0) continue;
      // Calculate displacement.
      dr[0] = x[id1][0] - x[id2][0];
      dr[1] = x[id1][1] - x[id2][1];
      dr[2] = x[id1][2] - x[id2][2];
      // Calculate squared distance
      rsqr = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      // Get radii
      sg1 = sg[id1];
      sg2 = sg[id2];
      // If close, interact.
      if (rsqr < sqr((sg1 + sg2)*cutoff)) kernel(id1, id2, sg1, sg2, rsqr, dr, f);
    }

    // --- Do verlet wrap part.
    if (verlet_wrap.empty()) return;
    
    // Get the bounds and boundary conditions
    RealType bnd_x = gflow->getBounds().wd(0);
    RealType bnd_y = gflow->getBounds().wd(1);
    RealType bnd_z = gflow->getBounds().wd(2);
    bool wrapX = Base::gflow->getBC(0)==BCFlag::WRAP;
    bool wrapY = Base::gflow->getBC(1)==BCFlag::WRAP;
    bool wrapZ = Base::gflow->getBC(1)==BCFlag::WRAP;

    // --- Go through all particles in verlet_wrap.
    for (int i=0; i<verlet_wrap.size(); i+=2) {
      int id1 = verlet_wrap[i];
      int id2 = verlet_wrap[i+1];
      // Check if the types are good
      if (type[id1]<0 || type[id2]<0) continue;
      dr[0] = x[id1][0] - x[id2][0];
      dr[1] = x[id1][1] - x[id2][1];
      dr[2] = x[id1][2] - x[id2][2];
      // Harmonic corrections to distance.
      if (wrapX) {
        RealType dX = bnd_x - fabs(dr[0]);
        if (dX<fabs(dr[0])) dr[0] = dr[0]>0 ? -dX : dX;
      }  
      if (wrapY) {
        RealType dY = bnd_y - fabs(dr[1]);
        if (dY<fabs(dr[1])) dr[1] = dr[1]>0 ? -dY : dY;
      } 
      if (wrapZ) {
        RealType dZ = bnd_z - fabs(dr[2]);
        if (dZ<fabs(dr[2])) dr[2] = dr[2]>0 ? -dZ : dZ;
      } 
      // Calculate squared distance
      rsqr = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      // Get radii
      sg1 = sg[id1];
      sg2 = sg[id2];
      // If close, interact.
      if (rsqr < sqr((sg1 + sg2)*cutoff)) kernel(id1, id2, sg1, sg2, rsqr, dr, f);
    }
  }

}