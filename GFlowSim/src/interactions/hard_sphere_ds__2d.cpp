#include "hard_sphere_ds__2d.hpp"
// Other files
#include "../utility/simd_generic.hpp"

namespace GFlowSimulation {

  HardSphereDs_2d::HardSphereDs_2d(GFlow *gflow) : HardSphereDs(gflow) {};

  bool HardSphereDs_2d::checks() {
    RealType **x = Base::simData->X();
    RealType **v = Base::simData->V();
    RealType **f = Base::simData->F();
    RealType *sg = Base::simData->Sg();
    int    *type = Base::simData->Type();
    // Check the condition
    return (x!=nullptr && v!=nullptr && f!=nullptr && sg!=nullptr && type!=nullptr && sim_dimensions==2);
  }

  void HardSphereDs_2d::interact() const {
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

    // Needed constants
    RealType sg1, sg2, dx, dy, rsqr, r, invr, Fn, dvx, dvy, vn;

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
      // If close, interact.
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

        // Calculate potential
        if (do_potential) {
          potential += 0.5*repulsion*sqr(r - sg1 - sg2);
        }
        // Calculate virial
        if (do_virial) {
          virial += repulsion*(sg1 + sg2 - r)*r;
        }
      }
    }

    // --- Do verlet wrap part.
    if (verlet_wrap.empty()) return;

    // Get the bounds and boundary conditions
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    bool wrapX = gflow->getBC(0)==BCFlag::WRAP;
    bool wrapY = gflow->getBC(1)==BCFlag::WRAP;
    // Extract bounds related data
    RealType bnd_x = bounds.wd(0);
    RealType bnd_y = bounds.wd(1);

    // --- Go through all particles
    for (int i=0; i<verlet_wrap.size(); i+=2) {
      // Get next pair of interacting particles.
      int id1 = verlet_wrap[i];
      int id2 = verlet_wrap[i+1];
      // Check if the types are good.
      if (type[id1]<0 || type[id2]<0) continue;
      // Calculate displacement.
      dx = x[id1][0] - x[id2][0];
      dy = x[id1][1] - x[id2][1];
      // Harmonic corrections to distance.
      if (wrapX) {
        RealType dX = bnd_x - fabs(dx);
        if (dX<fabs(dx)) dx = dx>0 ? -dX : dX;
      }  
      if (wrapY) {
        RealType dY = bnd_y - fabs(dy);
        if (dY<fabs(dy)) dy = dy>0 ? -dY : dY;
      }
      // Calculate squared distance.
      rsqr = dx*dx + dy*dy;
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

        // Calculate potential
        if (do_potential) {
          potential += 0.5*repulsion*sqr(r - sg1 - sg2);
        }
        // Calculate virial
        if (do_virial) {
          virial += repulsion*(sg1 + sg2 - r)*r;
        }
      }
    }
  }

}