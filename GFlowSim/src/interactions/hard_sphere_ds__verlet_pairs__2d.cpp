#include "hard_sphere_ds__verlet_pairs__2d.hpp"
// Other files
#include "../interactionhandlers/verletlist-pairs.hpp"
#include "../utility/simd_generic.hpp" // For un_clamp

namespace GFlowSimulation {

  HardSphereDs_VerletPairs_2d::HardSphereDs_VerletPairs_2d(GFlow *gflow) : HardSphereDs(gflow, new VerletListPairs(gflow)) {};

  void HardSphereDs_VerletPairs_2d::interact() const {
    // Common tasks
    HardSphereDs::interact();
    
    // Do dimensional check.
    // \todo Should probably have some sort of global error message system.
    if (sim_dimensions!=2) return;

    // Get the data pointers.
    RealType **x = Base::simData->X();
    RealType **v = Base::simData->V();
    RealType **f = Base::simData->F();
    RealType *sg = Base::simData->Sg();
    int    *type = Base::simData->Type();

    // Make sure all needed pointers are non null.
    // \todo Should probably have some sort of global error message system.
    if (x==nullptr || v==nullptr || f==nullptr || sg==nullptr || type==nullptr) return;

    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[2];
    copyVec(Base::gflow->getBCs(), boundaryConditions, 2); // Keep a local copy of the bcs

    // Needed constants
    RealType sg1, sg2, dx, dy, rsqr, r, invr, Fn, dvx, dvy, vn;
    // Point to the actual list from the verlet list object. Since we set the handler at initialization to 
    // be of type VerletListPairs, this cast should always succeed.
    vector<int> &verlet = dynamic_cast<VerletListPairs*>(handler)->verlet;

    // --- Go through all particles
    for (int i=0; i<verlet.size(); i+=2) {
      // Get next pair of interacting particles.
      int id1 = verlet[i];
      int id2 = verlet[i+1];

      // Check if the types are good.
      if (type[id1]<0 || type[id2]<0) continue;

      // Calculate displacement.
      dx = x[id1][0] - x[id2][0];
      dy = x[id1][1] - x[id2][1];
      // Calculate squared distance.
      rsqr = dx*dx + dy*dy;
      // Get radii
      sg1 = sg[id1];
      sg2 = sg[id2];

      if (rsqr < sqr(sg1 + sg2)) {
        // Calculate distance, inverse distance.
        r = sqrt(rsqr);
        invr = 1./r;

        // Create a normal vector.
        dx *= invr;
        dy *= invr;
        // Calculate the magnitude of the force.
        Fn = repulsion*(sg1 + sg2 - r);

        // Calculate relative velocity.
        dvx = v[id2][0] - v[id1][0];
        dvy = v[id2][1] - v[id1][1];
        // Caluclate normal velocity.
        vn = dvx*dx + dvy*dy;
        // Dissipation only occurs on loading the spring.
        Fn += dissipation * un_clamp(vn);

        // Update forces.
        f[id1][0] += Fn * dx;
        f[id2][0] -= Fn * dx;
        f[id1][1] += Fn * dy;
        f[id2][1] -= Fn * dy;
      }
    }
  }

}