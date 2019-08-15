#include "detector__3d.hpp"

namespace GFlowSimulation {

  Detector_3d::Detector_3d(GFlow *gflow) : Detector(gflow) {};

  Detector_3d::Detector_3d(GFlow *gflow, RealType k) : Detector(gflow, k) {};

  void Detector_3d::interact() const {
    // Common tasks
    Interaction::interact();

    // Do dimensional check.
    // \todo Should probably have some sort of global error message system.
    if (sim_dimensions!=3) return;

    // Get the data pointers.
    RealType **x = Base::simData->X();
    RealType **v = Base::simData->V();
    RealType *sg = Base::simData->Sg();
    RealType *im = Base::simData->Im();
    int    *type = Base::simData->Type();

    // Make sure all needed pointers are non null.
    // \todo Should probably have some sort of global error message system.
    if (x==nullptr || v==nullptr || sg==nullptr || im==nullptr || type==nullptr) return;

    // Needed constants
    RealType sg1, sg2, dx, dy, dz, rsqr;

    // --- Go through all particles
    for (int i=0; i<verlet.size(); i+=2) {
      // Get next pair of interacting particles.
      int id1 = verlet[i];
      int id2 = verlet[i+1];
      // Check if the types are good
      if (type[id1]<0 || type[id2]<0) continue;
      // Calculate displacement
      dx = x[id1][0] - x[id2][0];
      dy = x[id1][1] - x[id2][1];
      dz = x[id1][2] - x[id2][2];
      // Calculate squared distance.
      rsqr = dx*dx + dy*dy + dz*dz;
      // Get radii
      sg1 = sg[id1];
      sg2 = sg[id2];
      // If close, interact.
      if (rsqr < sqr((sg1 + sg2)*cutoff)) {
        // Calulate KE
        RealType ke = 0.5*sqr(v[id2], 3)/im[id2];
        // If too large 
        if (ke_threshold<ke) {
          gflow->terminate();
          return;
        }
      }
    }

    // --- Do verlet wrap part.
    if (verlet_wrap.size()==0) return;

    // Get the bounds and boundary conditions
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[3];
    copyVec(Base::gflow->getBCs(), boundaryConditions, 3); // Keep a local copy of the bcs
    // Extract bounds related data
    RealType bnd_x = bounds.wd(0);
    RealType bnd_y = bounds.wd(1);
    RealType bnd_z = bounds.wd(2);

    // --- Go through all particles
    for (int i=0; i<verlet_wrap.size(); i+=2) {
      // Get next pair of interacting particles.
      int id1 = verlet_wrap[i];
      int id2 = verlet_wrap[i+1];
      // Check if the types are good
      if (type[id1]<0 || type[id2]<0) continue;
      // Calculate displacement
      dx = x[id1][0] - x[id2][0];
      dy = x[id1][1] - x[id2][1];
      dz = x[id1][2] - x[id2][2];
      // Harmonic corrections to distance.
      if (boundaryConditions[0]==BCFlag::WRAP) {
        RealType dX = bnd_x - fabs(dx);
        if (dX<fabs(dx)) dx = dx>0 ? -dX : dX;
      }  
      if (boundaryConditions[1]==BCFlag::WRAP) {
        RealType dY = bnd_y - fabs(dy);
        if (dY<fabs(dy)) dy = dy>0 ? -dY : dY;
      } 
      if (boundaryConditions[2]==BCFlag::WRAP) {
        RealType dZ = bnd_z - fabs(dz);
        if (dZ<fabs(dz)) dz = dz>0 ? -dZ : dZ;
      }   
      // Calculate squared distance.
      rsqr = dx*dx + dy*dy + dz*dz;
      // Get radii
      sg1 = sg[id1];
      sg2 = sg[id2];
      // If close, interact.
      if (rsqr < sqr((sg1 + sg2)*cutoff)) {
        // Calulate KE
        RealType ke = 0.5*sqr(v[id2], 3)/im[id2];
        // If too large 
        if (ke_threshold<ke) {
          gflow->terminate();
          return;
        }
      }
    }
  }

}