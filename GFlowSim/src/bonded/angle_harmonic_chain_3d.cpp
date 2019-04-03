#include "angle_harmonic_chain_3d.hpp"

namespace GFlowSimulation {

  AngleHarmonicChain_3d::AngleHarmonicChain_3d(GFlow *gflow) : AngleHarmonicChain(gflow) {};

  AngleHarmonicChain_3d::AngleHarmonicChain_3d(GFlow *gflow, RealType K) : AngleHarmonicChain(gflow, K) {};

  void AngleHarmonicChain_3d::interact() const {
    // Call parent class.
    AngleHarmonicChain::interact();
    // Get simdata, check if the local ids need updating
    SimData *sd = Base::simData;
    RealType **x = sd->X();
    RealType **f = sd->F();
    if (sd->getNeedsRemake()) updateLocalIDs();
    // If there are not enough particles in the buffer
    if (local_ids.size()<3) return;
    // Variables and buffers
    RealType dx1, dy1, dz1, dx2, dy2, dz2, rsqr1, rsqr2, r1, r2, angle_cos, a11, a12, a22;
    RealType f1[3], f2[3];
    // Get the bounds and boundary conditions
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[3];
    copyVec(Base::gflow->getBCs(), boundaryConditions, 3); // Keep a local copy of the bcs
    // Extract bounds related data
    RealType bnd_x = bounds.wd(0);
    RealType bnd_y = bounds.wd(1);
    RealType bnd_z = bounds.wd(2);

    // If there are fewer than three particles, this loop will not execute
    for (int i=0; i+2<local_ids.size(); ++i) {
      // Get the global, then local ids of the particles.
      int id1 = local_ids[i], id2 = local_ids[i+1], id3=local_ids[i+2];
      // First bond
      dx1 = x[id1][0] - x[id2][0];
      dy1 = x[id1][1] - x[id2][1];
      dz1 = x[id1][2] - x[id2][2];
      // Second bond
      dx2 = x[id3][0] - x[id2][0];
      dy2 = x[id3][1] - x[id2][1];
      dz2 = x[id3][2] - x[id2][2];
      // Harmonic corrections to distance.
      if (boundaryConditions[0]==BCFlag::WRAP) {
        RealType dX = bnd_x - fabs(dx1);
        if (dX<fabs(dx1)) dx1 = dx1>0 ? -dX : dX;
        dX = bnd_x - fabs(dx2);
        if (dX<fabs(dx2)) dx2 = dx2>0 ? -dX : dX;
      }  
      if (boundaryConditions[1]==BCFlag::WRAP) {
        RealType dY = bnd_y - fabs(dy1);
        if (dY<fabs(dy1)) dy1 = dy1>0 ? -dY : dY;
        dY = bnd_y - fabs(dy2);
        if (dY<fabs(dy2)) dy2 = dy1>0 ? -dY : dY;
      } 
      if (boundaryConditions[2]==BCFlag::WRAP) {
        RealType dZ = bnd_z - fabs(dz1);
        if (dZ<fabs(dz1)) dz1 = dz1>0 ? -dZ : dZ;
        dZ = bnd_z - fabs(dy2);
        if (dZ<fabs(dz2)) dz2 = dz1>0 ? -dZ : dZ;
      } 
      // Get distance squared, distance
      rsqr1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
      r1 = sqrt(rsqr1);
      rsqr2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
      r2 = sqrt(rsqr2);
      // Find the cosine of the angle between the atoms
      angle_cos = dx1*dx2 + dy1*dy2 + dz1*dz2;
      angle_cos /= (r1*r2);
      // Calculate the force
      a11 = angleConstant*angle_cos / rsqr1;
      a12 = -angleConstant / (r1*r2);
      a22 = angleConstant*angle_cos / rsqr2;
      // Force buffers
      f1[0] = a11*dx1 + a12*dx2;
      f1[1] = a11*dy1 + a12*dy2;
      f1[2] = a11*dz1 + a12*dz2;
      f2[0] = a22*dx2 + a12*dx1;
      f2[1] = a22*dy2 + a12*dy1;
      f2[2] = a22*dz2 + a12*dz1;

      // First atom
      f[id1][0] += f1[0];
      f[id1][1] += f1[1];
      f[id1][2] += f1[2];
      // Second atom
      f[id2][0] -= (f1[0] + f2[0]);
      f[id2][1] -= (f1[1] + f2[1]);
      f[id2][2] -= (f1[2] + f2[2]);
      // Third atom
      f[id3][0] += f2[0];
      f[id3][1] += f2[1];
      f[id3][2] += f2[2];
    }
    // For two atoms, there is no angle. Just compute the bond force.
  }

}