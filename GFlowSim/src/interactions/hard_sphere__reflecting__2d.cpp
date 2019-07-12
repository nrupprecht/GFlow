#include "hard_sphere__reflecting__2d.hpp"

namespace GFlowSimulation {

  HardSphere_Reflecting_2d::HardSphere_Reflecting_2d(GFlow *gflow) : Interaction(gflow) {};

  void HardSphere_Reflecting_2d::interact() const {
    // Do dimensional check.
    // \todo Should probably have some sort of global error message system.
    if (sim_dimensions!=2) return;

    // Get the data pointers.
    RealType **x = Base::simData->X();
    RealType **v = Base::simData->V();
    RealType *sg = Base::simData->Sg();
    RealType *im = Base::simData->Im();
    int    *type = Base::simData->Type();

    // Make sure all needed pointers are non null.
    // \todo Should probably have some sort of global error message system.
    if (x==nullptr || v==nullptr || sg==nullptr || im==nullptr || type==nullptr) return;

    // Needed constants
    RealType sg1, sg2, dx, dy, rsqr, r, invr, magnitude;

    // --- Go through all particles
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

        // Initial projected velocities
        RealType v1 = v[id1][0]*dx + v[id1][1]*dy;
        RealType v2 = v[id2][0]*dx + v[id2][1]*dy;

        // Relative velocities
        RealType dvx = v[id1][0] - v[id2][0];
        RealType dvy = v[id1][1] - v[id2][1];

        // Only interact if particles are going towards each other.
        if (dx*dvx+dy*dvy>0) continue;

        // Calculate final velocities based on momentum.
        if (im[id1]>0 && im[id2]>0) {
          // Mass related constants
          RealType m1 = 1./im[id1], m2 = 1./im[id2];
          RealType mT = m1+m2, mD = m1 - m2;
          RealType invMT = 1./mT;
          // Final velocities (in the projected direction)
          RealType vf1 = mD*invMT * v1 + 2*m2*invMT * v2;
          RealType vf2 = 2*m1*invMT * v1 - mD*invMT * v2;
          // Adjust velocities
          v[id1][0] -= (v1 - vf1)*dx;
          v[id1][1] -= (v1 - vf1)*dy;
          v[id2][0] -= (v2 - vf2)*dx;
          v[id2][1] -= (v2 - vf2)*dy;
        }
        // Both particles have infinite mass. No force.
        else if (im[id1]==0 && im[id2]==0);
        // One particle has infinite mass. 
        else {
          // Elastic collision results in normal momenta reflecting
          v[id1][0] -= 2*v1*dx;
          v[id1][1] -= 2*v1*dy;
          v[id2][0] -= 2*v2*dx;
          v[id2][1] -= 2*v2*dy;
        }
      }
    }

    // --- Do verlet wrap part.
    if (verlet_wrap.empty()) return;

    // Get the bounds and boundary conditions
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[2];
    copyVec(Base::gflow->getBCs(), boundaryConditions, 2); // Keep a local copy of the bcs
    // Extract bounds related data
    RealType bnd_x = bounds.wd(0);
    RealType bnd_y = bounds.wd(1);

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

        // Initial projected velocities
        RealType v1 = v[id1][0]*dx + v[id1][1]*dy;
        RealType v2 = v[id2][0]*dx + v[id2][1]*dy;

        // Relative velocities
        RealType dvx = v[id1][0] - v[id2][0];
        RealType dvy = v[id1][1] - v[id2][1];

        // Only interact if particles are going towards each other.
        if (dx*dvx+dy*dvy>0) continue;

        // Calculate final velocities based on momentum.
        if (im[id1]>0 && im[id2]>0) {
          // Mass related constants
          RealType m1 = 1./im[id1], m2 = 1./im[id2];
          RealType mT = m1+m2, mD = m1 - m2;
          RealType invMT = 1./mT;
          // Final velocities (in the projected direction)
          RealType vf1 = mD*invMT * v1 + 2*m2*invMT * v2;
          RealType vf2 = 2*m1*invMT * v1 - mD*invMT * v2;
          // Adjust velocities
          v[id1][0] -= (v1 - vf1)*dx;
          v[id1][1] -= (v1 - vf1)*dy;
          v[id2][0] -= (v2 - vf2)*dx;
          v[id2][1] -= (v2 - vf2)*dy;
        }
        // Both particles have infinite mass. No force.
        else if (im[id1]==0 && im[id2]==0);
        // One particle has infinite mass. 
        else {
          // Elastic collision results in normal momenta reflecting
          v[id1][0] -= 2*v1*dx;
          v[id1][1] -= 2*v1*dy;
          v[id2][0] -= 2*v2*dx;
          v[id2][1] -= 2*v2*dy;
        }
      }
    }
  }

}