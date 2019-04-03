#include "harmonicbond_2d.hpp"

namespace GFlowSimulation {

  HarmonicBond_2d::HarmonicBond_2d(GFlow *gflow) : HarmonicBond(gflow) {};

  HarmonicBond_2d::HarmonicBond_2d(GFlow *gflow, RealType K) : HarmonicBond(gflow, K) {};

  void HarmonicBond_2d::interact() const {
    // Call parent class.
    HarmonicBond::interact();
    // This object only acts in two dimensions.
    if (sim_dimensions!=2) return;
    // Get the number of bonds. left and right will have the same size - we checked this last time
    // updateLocalIDs was called.
    int nbonds = left.size();

    // Check if local ids need updating.
    if (simData->getNeedsRemake()) updateLocalIDs();

    // Get simdata, check if the local ids need updating
    RealType **x = simData->X();
    RealType **f = simData->F();
    int    *type = simData->Type();
    // Displacement vector
    RealType dx, dy, r, invr, dr, scale;
    // Get the bounds and boundary conditions
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[2];
    copyVec(Base::gflow->getBCs(), boundaryConditions, 2); // Keep a local copy of the bcs
    // Extract bounds related data
    RealType bnd_x = bounds.wd(0);
    RealType bnd_y = bounds.wd(1);

    // Go through all bonds, calculating forces.
    for (int i=0; i<nbonds; ++i) {
      // Get the global, then local ids of the particles.
      int id1 = left[i], id2 = right[i];
      // Check if the types are good
      if (type[id1]<0 || type[id2]<0) continue;
      // Calculate displacement
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
      // Get the distance, inverse distance
      r = sqrtf(dx*dx + dy*dy);
      invr = 1./r;
      // Calculate displacement from equilibrium
      dr = r - distance[i];
      // Makes dX into the force. The springConstant*dr comes from the force, the invr comes from 
      // normalizing {dx, dy}
      scale = springConstant*dr*invr;
      dx *= scale;
      dy *= scale;
      // Add forces to particles
      f[id1][0] -= dx;
      f[id1][1] -= dy;
      f[id2][0] += dx;
      f[id2][1] += dy;
    }
  }

}