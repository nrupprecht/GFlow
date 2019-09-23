#include "rigidbody_2d.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  RigidBody_2d::RigidBody_2d(GFlow *gflow) : Body(gflow), Group(gflow) {};

  void RigidBody_2d::post_forces() {
    // If size is zero, this body is empty.
    if (size()==0) return;
    // If needed, remake local ids.
    if (simData->getNeedsRemake()) update_local_ids();

    // Get the position, force, and inverse mass arrays.
    auto x = simData->X();
    auto f = simData->F();
    auto im = simData->Im();

    // Get the bounds and boundary conditions
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[2];
    copyVec(Base::gflow->getBCs(), boundaryConditions, 2); // Keep a local copy of the bcs
    // Extract bounds related data
    RealType bnd_x = bounds.wd(0);
    RealType bnd_y = bounds.wd(1);

    // Calculate center of mass position. We have to be careful, since there may be harmonic boundary conditions.
    RealType rcm[2], r0[2], fnet[2], I(0), tau(0);
    zeroVec(rcm,  2);
    zeroVec(fnet, 2);
    // Set the reference point to be the position of the first particle.
    copyVec(x[local_ids[0]], r0, 2);
    // Go through all body particles.
    for (auto id : local_ids) {
      RealType dx = x[id][0] - r0[0];
      RealType dy = x[id][1] - r0[1];
      // Harmonic corrections to distance.
      if (boundaryConditions[0]==BCFlag::WRAP) {
        RealType dX = bnd_x - fabs(dx);
        if (dX<fabs(dx)) dx = dx>0 ? -dX : dX;
      }  
      if (boundaryConditions[1]==BCFlag::WRAP) {
        RealType dY = bnd_y - fabs(dy);
        if (dY<fabs(dy)) dy = dy>0 ? -dY : dY;
      }  
      // Get the mass of the particle - assumes none of the particles have infinite mass.
      RealType m = 1./im[id];
      // Update rcm
      rcm[0] += m*dx;
      rcm[1] += m*dy;
      // Update fnet
      fnet[0] += f[id][0];
      fnet[1] += f[id][1];
      // Update moment of inertia
      I += m*(dx*dx + dy*dy);
      // Update torque
      tau += f[id][0]*dy - f[id][1]*dx;
    }
    // Divide by total mass and add to reference position to get the actual com position.
    rcm[0] /= mass;
    rcm[1] /= mass;
    rcm[0] += r0[0];
    rcm[1] += r0[1];
    // Correct r0 for harmonic boundary conditions
    if (boundaryConditions[0]==BCFlag::WRAP) {
      if (rcm[0]<bounds.min[0]) rcm[0] += bnd_x;
      else if (bounds.max[0]<rcm[0]) rcm[0] -= bnd_x;
    }
    if (boundaryConditions[1]==BCFlag::WRAP) {
      if (rcm[1]<bounds.min[1]) rcm[1] += bnd_y;
      else if (bounds.max[1]<rcm[1]) rcm[1] -= bnd_y;
    }
    // We now have the correct center of mass. Calcuate linear and angular acceleration.
    RealType A_lin[2] = { fnet[0]/mass, fnet[1]/mass };
    RealType alpha = tau/I;

    // Set the forces for each particle.
    for (auto id : local_ids) {
      RealType dx = x[id][0] - r0[0];
      RealType dy = x[id][1] - r0[1];
      // Harmonic corrections to distance.
      if (boundaryConditions[0]==BCFlag::WRAP) {
        RealType dX = bnd_x - fabs(dx);
        if (dX<fabs(dx)) dx = dx>0 ? -dX : dX;
      }  
      if (boundaryConditions[1]==BCFlag::WRAP) {
        RealType dY = bnd_y - fabs(dy);
        if (dY<fabs(dy)) dy = dy>0 ? -dY : dY;
      } 
      // Get the mass of the particle - assumes none of the particles have infinite mass.
      RealType m = 1./im[id];
      // Set the force on the particle.
      f[id][0] = m*(A_lin[0] - alpha*dy);
      f[id][1] = m*(A_lin[1] + alpha*dx);
    }
  }

  void RigidBody_2d::add(int i) {
    Group::add(i);
    // Update mass
    mass += 1./simData->Im(simData->getLocalID(i));
  }

}