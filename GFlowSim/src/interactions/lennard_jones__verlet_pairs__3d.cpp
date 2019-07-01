#include "lennard_jones__verlet_pairs__3d.hpp"

namespace GFlowSimulation {

  LennardJones_VerletPairs_3d::LennardJones_VerletPairs_3d(GFlow *gflow) : LennardJones(gflow) {};

  void LennardJones_VerletPairs_3d::interact() const {
    // Common tasks
    LennardJones::interact();

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

    // Get the bounds and boundary conditions
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[3];
    copyVec(Base::gflow->getBCs(), boundaryConditions, 3); // Keep a local copy of the bcs
    // Extract bounds related data
    RealType bnd_x = bounds.wd(0);
    RealType inv_bnd_x = 1./bnd_x;
    RealType bnd_y = bounds.wd(1);
    RealType inv_bnd_y = 1./bnd_y;
    RealType bnd_z = bounds.wd(2);
    RealType inv_bnd_z = 1./bnd_z;

    // Needed constants
    RealType sg1, sg2, dx, dy, dz, rsqr, r, invr, gamma, g3, g6, g12, magnitude;

    // --- Go through all particles
    for (int i=0; i<verlet.size(); i+=2) {
      int id1 = verlet[i];
      int id2 = verlet[i+1];
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
      // Calculate squared distance
      rsqr = dx*dx + dy*dy + dz*dz;
      // Get radii
      sg1 = sg[id1];
      sg2 = sg[id2];
      // If close, interact. 
      if (rsqr < sqr((sg1 + sg2)*cutoff)) {
        // Calculate distance, inverse distance.
        r = sqrt(rsqr);
        invr = 1./r;
        // Create a normal vector
        dx *= invr;
        dy *= invr;
        dz *= invr;
        // Calculate the magnitude of the force
        gamma = (sg1+sg2)*invr;
        g3  = gamma*gamma*gamma; 
        g6  = g3*g3;
        g12 = g6*g6;
        // Calculate magnitude
        magnitude = 24.*strength*(2.*g12 - g6)*invr;
        // Update forces
        f[id1][0] += magnitude * dx;
        f[id1][1] += magnitude * dy;
        f[id1][2] += magnitude * dz;
        f[id2][0] -= magnitude * dx;
        f[id2][1] -= magnitude * dy;
        f[id2][2] -= magnitude * dz;

        // Calculate potential
        if (do_potential) {
          potential += 4.*strength*(g12 - g6) - potential_energy_shift;
        }
        // Calculate virial
        if (do_virial) {
          virial += magnitude*r;
        }
      }
    }
  }

}