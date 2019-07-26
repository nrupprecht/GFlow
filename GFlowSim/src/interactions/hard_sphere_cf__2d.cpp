#include "hard_sphere_cf__2d.hpp"
// Other files
#include "../utility/simd_generic.hpp" // For un_clamp
#include "../integrators/angularvelocityverlet2d.hpp"

namespace GFlowSimulation {

  HardSphereCf_2d::HardSphereCf_2d(GFlow *gflow) : HardSphereCf(gflow) {
    simData->requestScalarData("Om");
    simData->requestScalarData("Tq");
    // Add an angular integrator.
    gflow->addIntegrator(new AngularVelocityVerlet2d(gflow));
  };


  void HardSphereCf_2d::interact() const {
    // Common tasks
    HardSphereCf::interact();

    // Do dimensional check.
    // \todo Should probably have some sort of global error message system.
    if (sim_dimensions!=2) return;

    // Get the data pointers.
    RealType **x = simData->X();
    RealType **v = Base::simData->V();
    RealType **f = simData->F();
    RealType *sg = simData->Sg();
    int om_add = simData->getScalarData("Om");
    int tq_add = simData->getScalarData("Tq");
    RealType *om = simData->ScalarData(om_add);
    RealType *tq = simData->ScalarData(tq_add);
    int    *type = simData->Type();

    // Make sure all needed pointers are non null.
    // \todo Should probably have some sort of global error message system.
    if (x==nullptr || v==nullptr || f==nullptr || sg==nullptr || 
      type==nullptr || om==nullptr || tq==nullptr) return;

    // Needed constants.
    RealType sg1, sg2, dx, dy, rsqr, r, invr, Fn, dvx, dvy, vn, dOm;

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
        // Calculate the positive tangential direction.
        RealType dtx = dy, dty = -dx;
        // Calculate angular velocities
        RealType dvs = om[id1]*sg[id1] + om[id2]*sg[id2]; // Relative surface velocity due to angular velocity.
        RealType vt = dvx*dtx + dvy*dty;
        RealType Ft = mu*Fn*sign(dvs + vt);
        // Update forces.
        RealType Fx = Fn * dx + Ft * dtx;
        RealType Fy = Fn * dy + Ft * dty;
        f[id1][0] += Fx;
        f[id2][0] -= Fx;
        f[id1][1] += Fy;
        f[id2][1] -= Fy;
        // Update torques
        tq[id1] -= Ft*sg[id1];
        tq[id2] -= Ft*sg[id2];

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
        // Calculate the positive tangential direction.
        RealType dtx = dy, dty = -dx;
        // Calculate angular velocities
        RealType dvs = om[id1]*sg[id1] + om[id2]*sg[id2]; // Relative surface velocity due to angular velocity.
        RealType vt = dvx*dtx + dvy*dty;
        RealType Ft = mu*Fn*sign(dvs + vt);
        // Update forces.
        RealType Fx = Fn * dx + Ft * dtx;
        RealType Fy = Fn * dy + Ft * dty;
        f[id1][0] += Fx;
        f[id2][0] -= Fx;
        f[id1][1] += Fy;
        f[id2][1] -= Fy;
        // Update torques
        tq[id1] -= Ft*sg[id1];
        tq[id2] -= Ft*sg[id2];

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