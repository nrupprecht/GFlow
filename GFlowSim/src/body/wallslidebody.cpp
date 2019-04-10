#include "wallslidebody.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  WallSlideBody::WallSlideBody(GFlow *gflow, int d) : Body(gflow), slide_dimension(d) {};

  WallSlideBody::WallSlideBody(GFlow *gflow, Group group, int d) : Body(gflow), slide_dimension(d) {
    // Set the group part of this body.
    *dynamic_cast<Group*>(this) = group;
  }

  void WallSlideBody::pre_integrate() {
    // Set all velocities to zero
    RealType **v = simData->V();
    for (auto id : local_ids)
      zeroVec(v[id], sim_dimensions);
  }

  void WallSlideBody::correct() {
    // Make sure slide dimension is valid
    if (slide_dimension<0 || sim_dimensions<=slide_dimension) return;
    // If necessary, update local ids.
    if (simData->getNeedsRemake()) update_local_ids(simData);

    // Get the force and inverse mass arrays.
    RealType **f = simData->F();
    RealType *im = simData->Im();

    // Net force on and total mass of the body
    RealType Fnet = 0;
    RealType M = 0;
    bool has_inf = false;

    // Accumulation loop
    for (const auto id : local_ids) {
      // Accumulate mass
      if (im[id]>0) M += 1./im[id];
      else has_inf = true;
      // Accumulate force
      Fnet += f[id][slide_dimension];
      // Clear the force vector.
      zeroVec(f[id], sim_dimensions);
    }

    // Check if any object had infinite mass
    if (has_inf) return; // There is an infinitely massive object.

    // Net acceleration
    RealType A = Fnet/M;
    // Force setting loop
    for (const auto id : local_ids) {
      // Set the slide dimension component of force.
      f[id][slide_dimension] = 1./im[id] * A;
    }
  }


}