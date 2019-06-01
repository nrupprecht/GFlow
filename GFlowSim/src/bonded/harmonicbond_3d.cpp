#include "harmonicbond_3d.hpp"

namespace GFlowSimulation {

  HarmonicBond_3d::HarmonicBond_3d(GFlow *gflow) : HarmonicBond(gflow) {};

  HarmonicBond_3d::HarmonicBond_3d(GFlow *gflow, RealType K) : HarmonicBond(gflow, K) {};

  void HarmonicBond_3d::interact() const {
    // Call parent class.
    HarmonicBond::interact();
    // This object only acts in two dimensions.
    if (sim_dimensions!=3) return;
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
    RealType dx, dy, dz, r, invr, dr, scale;
    // Get the bounds and boundary conditions
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[3];
    copyVec(Base::gflow->getBCs(), boundaryConditions, 3); // Keep a local copy of the bcs
    // Extract bounds related data
    RealType bnd_x = bounds.wd(0);
    RealType bnd_y = bounds.wd(1);
    RealType bnd_z = bounds.wd(2);

    // Go through all bonds, calculating forces.
    for (int i=0; i<nbonds; ++i) {
      // Get the global, then local ids of the particles.
      int id1 = left[i], id2 = right[i];
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
      // Get the distance, inverse distance
      r = sqrtf(dx*dx + dy*dy + dz*dz);
      invr = 1./r;
      // Calculate displacement from equilibrium
      dr = r - distance[i];
      // Makes dX into the force. The springConstant*dr comes from the force, the invr comes from 
      // normalizing {dx, dy, dz}
      scale = springConstant*dr*invr;
      dx *= scale;
      dy *= scale;
      dz *= scale;
      // Add forces to particles
      f[id1][0] -= dx;
      f[id1][1] -= dy;
      f[id1][2] -= dz;
      f[id2][0] += dx;
      f[id2][1] += dy;
      f[id2][2] += dz;

      // Potential energy
      if (do_potential) {
        potential += 0.5*springConstant*sqr(dr);
      }
    }
  }

}