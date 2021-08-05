#include <body/wallslidebody.hpp>
// Other files
#include <utility/vectormath.hpp>
#include <base/simdata.hpp>

using namespace GFlowSimulation;

WallSlideBody::WallSlideBody(GFlow *gflow, int d)
    : Body(gflow), Group(gflow), slide_dimension(d) {};

WallSlideBody::WallSlideBody(GFlow *gflow, const Group &group, int d)
    : Body(gflow), Group(gflow), slide_dimension(d) {
  // Set the group part of this body.
  *dynamic_cast<Group *>(this) = group;
}

void WallSlideBody::pre_integrate() {
  // Return if no particles in the wall
  if (global_ids.empty()) {
    return;
  }

  // Update local ids
  update_local_ids();

  // Set all velocities to zero
  auto v = simData->V();
  for (auto id : local_ids) {
    zeroVec(v(id), sim_dimensions);
  }

  // Check length
  auto x = simData->X();
  auto sg = simData->Sg();
  Vec max(sim_dimensions), min(sim_dimensions);
  copyVec(x(local_ids[0]), min);
  max = min;
  for (auto id : local_ids) {
    for (int d = 0; d < sim_dimensions; ++d) {
      RealType xd = x(id, d);
      RealType r = sg(id);
      if (xd - r < min[d]) {
        min[d] = xd - r;
      }
      else if (max[d] < xd + r) {
        max[d] = xd + r;
      }
    }
  }
  // Set length
  Vec dl = max - min;
  length = sqrt(dl * dl);
}

void WallSlideBody::correct() {
  // Make sure slide dimension is valid
  if (slide_dimension < 0 || sim_dimensions <= slide_dimension) {
    return;
  }
  // If necessary, update local ids.
  if (simData->getNeedsRemake()) {
    update_local_ids();
  }

  // Get the force and inverse mass arrays.
  auto f = simData->F();
  auto im = simData->Im();

  // Net force on and total mass of the body
  Fnet = 0;
  RealType M = 0;
  bool has_inf = false;

  // Accumulation loop
  for (const auto id : local_ids) {
    // Accumulate mass
    if (0 < im(id)) {
      M += 1. / im(id);
    }
    else {
      has_inf = true;
    }
    // Accumulate force
    Fnet += f(id, slide_dimension);
    // Clear the force vector.
    zeroVec(f(id), sim_dimensions);
  }

  // Check if any object had infinite mass
  if (has_inf) {
    return; // There is an infinitely massive object.
  }

  // Net acceleration
  RealType A = Fnet / M;
  // Force setting loop
  for (const auto id : local_ids) {
    // Set the slide dimension component of force.
    f(id, slide_dimension) = 1. / im(id) * A;
  }
}

RealType WallSlideBody::getPosition() {
  return simData->X(simData->getLocalID(global_ids[0]), slide_dimension);
}

RealType WallSlideBody::getVelocity() {
  return simData->V(simData->getLocalID(global_ids[0]), slide_dimension);
}

RealType WallSlideBody::getLength() const {
  return length;
}

RealType WallSlideBody::getFnet() const {
  return Fnet;
}
