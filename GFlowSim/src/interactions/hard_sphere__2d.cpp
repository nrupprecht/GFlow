#include "hard_sphere__2d.hpp"

namespace GFlowSimulation {

  HardSphere_2d::HardSphere_2d(GFlow *gflow) : HardSphere(gflow) {

    setRepulsion(5*DEFAULT_HARD_SPHERE_REPULSION);

  };

  void HardSphere_2d::interact() const {
    // Common tasks
    HardSphere::interact();

    // Do dimensional check.
    // \todo Should probably have some sort of global error message system.
    if (sim_dimensions!=2) return;

    // Get the data pointers.
    RealType **x = Base::simData->X();
    RealType **f = Base::simData->F();
    RealType *sg = Base::simData->Sg();
    int    *type = Base::simData->Type();

    // Make sure all needed pointers are non null.
    // \todo Should probably have some sort of global error message system.
    if (x==nullptr || f==nullptr || sg==nullptr || type==nullptr) return;

    // Needed constants
    RealType sg1, sg2, dx, dy, rsqr, r, invr, magnitude;

    // --- Go through all particles in verlet.
    for (int i=0; i<verlet.size(); i+=2) {
      int id1 = verlet[i];
      int id2 = verlet[i+1];
      // Check if the types are good
      if (type[id1]<0 || type[id2]<0) continue;
      // Calculate displacement.
      dx = x[id1][0] - x[id2][0];
      dy = x[id1][1] - x[id2][1];
      // Calculate squared distance
      rsqr = dx*dx + dy*dy;
      // Get radii
      sg1 = sg[id1];
      sg2 = sg[id2];
      // If close, interact.
      if (rsqr < sqr(sg1 + sg2)) {
        // Calculate distance, inverse distance.
        r = sqrt(rsqr);
        invr = 1./r;
        // Create a normal vector
        dx *= invr;
        dy *= invr;
        // Calculate the magnitude of the force
        magnitude = repulsion*(sg1 + sg2 - r);
        // Update forces
        f[id1][0] += magnitude * dx;
        f[id2][0] -= magnitude * dx;
        f[id1][1] += magnitude * dy;
        f[id2][1] -= magnitude * dy;
        // Calculate potential
        if (do_potential) {
          potential += 0.5*repulsion*sqr(r - sg1 - sg2);
        }
        // Calculate virial
        if (do_virial) {
          virial += magnitude*r;
        }
      }
    }

    // --- Do verlet wrap part.
    if (verlet_wrap.empty()) return;
    
    // Get the bounds and boundary conditions
    BCFlag boundaryConditions[2];
    copyVec(Base::gflow->getBCs(), boundaryConditions, 2); // Keep a local copy of the bcs
    // Extract bounds related data
    RealType bnd_x = gflow->getBounds().wd(0);
    RealType bnd_y = gflow->getBounds().wd(1);
    // --- Go through all particles in verlet_wrap.
    for (int i=0; i<verlet_wrap.size(); i+=2) {
      int id1 = verlet_wrap[i];
      int id2 = verlet_wrap[i+1];
      // Check if the types are good
      if (type[id1]<0 || type[id2]<0) continue;
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
      if (rsqr < sqr(sg1 + sg2)) {
        // Calculate distance, inverse distance.
        r = sqrt(rsqr);
        invr = 1./r;
        // Create a normal vector
        dx *= invr;
        dy *= invr;
        // Calculate the magnitude of the force
        magnitude = repulsion*(sg1 + sg2 - r);
        // Update forces
        f[id1][0] += magnitude * dx;
        f[id2][0] -= magnitude * dx;
        f[id1][1] += magnitude * dy;
        f[id2][1] -= magnitude * dy;
        // Calculate potential
        if (do_potential) {
          potential += 0.5*repulsion*sqr(r - sg1 - sg2);
        }
        // Calculate virial
        if (do_virial) {
          virial += magnitude*r;
        }
      }
    }
  }

  void HardSphere_2d::kernel(int id1, int id2, RealType R1, RealType R2, RealType rsqr, RealType *dr, RealType **f) const {
    RealType r = sqrt(rsqr);
    RealType invr = 1./r;
    // Create a normal vector
    dr[0] *= invr;
    dr[1] *= invr;
    // Calculate the magnitude of the force
    RealType magnitude = repulsion*(R1 + R2 - r);

    // Update forces
    f[id1][0] += magnitude * dr[0];
    f[id2][0] -= magnitude * dr[0];
    f[id1][1] += magnitude * dr[1];
    f[id2][1] -= magnitude * dr[1];

    // Calculate potential
    if (do_potential) {
      potential += 0.5*repulsion*sqr(r - R1 - R2);
    }
    // Calculate virial
    if (do_virial) {
      virial += magnitude*r;
    }
  }

}