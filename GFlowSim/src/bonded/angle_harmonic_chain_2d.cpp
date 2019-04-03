#include "angle_harmonic_chain_2d.hpp"

namespace GFlowSimulation {

  AngleHarmonicChain_2d::AngleHarmonicChain_2d(GFlow *gflow) : AngleHarmonicChain(gflow) {};

  AngleHarmonicChain_2d::AngleHarmonicChain_2d(GFlow *gflow, RealType K) : AngleHarmonicChain(gflow, K) {};

  void AngleHarmonicChain_2d::interact() const {
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
    RealType dx1, dy1, dx2, dy2, rsqr1, rsqr2, r1, r2, invr1, invr2, angle_cos, a11, a12, a22;
    RealType f1[2], f3[2];
    // Get the bounds and boundary conditions
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[2];
    copyVec(Base::gflow->getBCs(), boundaryConditions, 2); // Keep a local copy of the bcs
    // Extract bounds related data
    RealType bnd_x = bounds.wd(0);
    RealType bnd_y = bounds.wd(1);

    // If there are fewer than three particles, this loop will not execute
    for (int i=0; i+2<local_ids.size(); ++i) {
      // Get the global, then local ids of the particles.
      int id1 = local_ids[i], id2 = local_ids[i+1], id3=local_ids[i+2];
      // First bond
      dx1 = x[id1][0] - x[id2][0];
      dy1 = x[id1][1] - x[id2][1];
      // Second bond
      dx2 = x[id3][0] - x[id2][0];
      dy2 = x[id3][1] - x[id2][1];
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
      // Get distance squared, distance
      rsqr1 = dx1*dx1 + dy1*dy1;
      r1 = sqrt(rsqr1);
      rsqr2 = dx2*dx2 + dy2*dy2;
      r2 = sqrt(rsqr2);
      // Find the cosine of the angle between the atoms
      angle_cos = dx1*dx2 + dy1*dy2;
      angle_cos /= (r1*r2);
      // Calculate the force
      a11 = angleConstant * angle_cos / rsqr1;
      a12 = -angleConstant / (r1*r2);
      a22 = angleConstant * angle_cos / rsqr2;
      // Force buffers
      f1[0] = a11*dx1 + a12*dx2;
      f1[1] = a11*dy1 + a12*dy2;
      f3[0] = a22*dx2 + a12*dx1;
      f3[1] = a22*dy2 + a12*dy1;
      // First atom
      f[id1][0] += f1[0];
      f[id1][1] += f1[1];
      // Second atom
      f[id2][0] -= (f1[0] + f3[0]);
      f[id2][1] -= (f1[1] + f3[1]);
      // Third atom
      f[id3][0] += f3[0];
      f[id3][1] += f3[1];
    }
    // For two atoms, there is no angle. Just compute the bond force.
  }

}