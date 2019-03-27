#include "hard_sphere__verlet_pairs__3d.hpp"
// Other files
#include "../interactionhandlers/verletlist-pairs.hpp"

namespace GFlowSimulation {

  HardSphere_VerletPairs_3d::HardSphere_VerletPairs_3d(GFlow *gflow) : HardSphere(gflow, new VerletListPairs(gflow)) {};

  void HardSphere_VerletPairs_3d::interact() const {
    // Common tasks
    HardSphere::interact();

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
    RealType sg1, sg2, dx, dy, dz, rsqr, r, invr, Fn;
    // Point to the actual list from the verlet list object. Since we set the handler at initialization to 
    // be of type VerletListPairs, this cast should always succeed.
    vector<int> &verlet = dynamic_cast<VerletListPairs*>(handler)->verlet;

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
      if (rsqr < sqr(sg1 + sg2)) {
        // Calculate distance, inverse distance.
        r = sqrt(rsqr);
        invr = 1./r;
        // Create a normal vector.
        dx *= invr;
        dy *= invr;
        dz *= invr;
        // Calculate the magnitude of the force.
        Fn = repulsion*(sg1 + sg2 - r);
        // Update forces
        f[id1][0] += Fn * dx;
        f[id2][0] -= Fn * dx;
        f[id1][1] += Fn * dy;
        f[id2][1] -= Fn * dy;
        f[id1][2] += Fn * dz;
        f[id2][2] -= Fn * dz;
        
        // Calculate potential
        if (do_potential) {
          potential += 0.5*repulsion*sqr(r - sg1 - sg2);
        }
        // Calculate virial
        if (do_virial) {
          virial += Fn*r;
        }
      }
    }
  }

}