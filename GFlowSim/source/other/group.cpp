#include "other/group.hpp"

#include <utility>
// Other files
#include "base/simdata.hpp"

using namespace GFlowSimulation;

Group::Group(shared_ptr<SimData> sim_data) {
  this->sim_data = std::move(sim_data);
}

Group::Group(GFlow *gflow) {
  sim_data = gflow->getSimData();
}

Group::Group(const Group &g) {
  global_ids = g.global_ids;
  local_ids = g.local_ids;
  sim_data = g.sim_data;
}

int Group::size() const {
  return global_ids.size();
}

bool Group::empty() const {
  return global_ids.empty();
}

int Group::at(int i) const {
  return local_ids.at(i);
}

int Group::g_at(int i) const {
  return global_ids.at(i);
}

void Group::findCenterOfMass(RealType *rcm) const {
  int sim_dimensions = sim_data->getSimDimensions();
  zeroVec(rcm, sim_dimensions);
  if (size() == 0) {
    return;
  }

  // Dispacement vector
  Vec dr(sim_dimensions), r0(sim_dimensions), Rcm(sim_dimensions);

  // Get the bounds and boundary conditions
  auto gflow = sim_data->getGFlow();
  Bounds bounds = gflow->getBounds(); // Simulation bounds

  // Get position, mass arrays
  auto x = sim_data->X();
  auto im = sim_data->Im();
  RealType mass = 0;

  // Set reference point
  int mid = floor(local_ids.size() / 2);
  copyVec(x(local_ids[mid]), r0.data, sim_dimensions);
  // Go through all particles
  for (auto id : local_ids) {
    // Get displacement between the particle and the reference point.
    gflow->getDisplacement(x(id), r0.data, dr.data);
    // Get the mass of the particle - assumes none of the particles have infinite mass.
    RealType m = 1. / im(id);
    mass += m;
    // Update Rcm
    plusEqVecScaled(Rcm.data, dr.data, m, sim_dimensions);
  }
  // Divide by total mass and add to reference position to get the actual com position.
  scalarMultVec(1. / mass, Rcm.data, sim_dimensions);
  // Add to reference point, and do wrapping
  Rcm += r0;

  // Harmonic boundary conditions.
  for (int d = 0; d < sim_dimensions; ++d) {
    if (Rcm[d] < bounds.min[d]) {
      Rcm[d] += bounds.wd(d);
    }
    else if (bounds.max[d] < Rcm[d]) {
      Rcm[d] -= bounds.wd(d);
    }
  }

  copyVec(Rcm, rcm);
}

void Group::findCOMVelocity(RealType *v) const {
  int sim_dimensions = sim_data->getSimDimensions();
  zeroVec(v, sim_dimensions);
  if (size() == 0) {
    return;
  }
  // Zero vector
  zeroVec(v, sim_dimensions);
  // Compute new momentum
  RealType mass = 0;
  for (auto id : local_ids) {
    RealType m = 1. / sim_data->Im(id);
    plusEqVecScaled(v, sim_data->V(id), m, sim_dimensions);
    mass += m;
  }
  // Normalize
  scalarMultVec(1. / mass, v, sim_dimensions);
}

void Group::addVelocity(RealType *V) const {
  if (size() == 0) {
    return;
  }
  // Get the dimensionality
  int sim_dimensions = sim_data->getSimDimensions();
  // Get arrays
  auto v = sim_data->V();
  for (auto id : local_ids) {
    plusEqVec(v(id), V, sim_dimensions);
  }
}

void Group::addAcceleration(RealType *A) const {
  if (size() == 0) {
    return;
  }
  // Get the dimensionaliry
  int sim_dimensions = sim_data->getSimDimensions();
  // Get arrays
  auto f = sim_data->F();
  auto im = sim_data->Im();
  for (auto id : local_ids) {
    if (im(id) > 0) {
      RealType mass = 1. / im(id);
      plusEqVecScaled(f(id), A, mass, sim_dimensions);
    }
  }
}

int Group::getIndex(int id) const {
  auto it = std::find(local_ids.begin(), local_ids.end(), id);
  if (it == local_ids.end()) {
    return -1;
  }
  // Otherwise, we found the element
  return std::distance(local_ids.begin(), it);
}

bool Group::contains(int id) const {
  return std::find(local_ids.begin(), local_ids.end(), id) != local_ids.end();
}

RealType Group::getChainedLength() const {
  RealType length = 0;
  auto x = sim_data->X();
  auto gflow = sim_data->getGFlow();
  // Compute distances between adjacent (in the order defined by the id list) particles.
  for (int i = 0; i < local_ids.size() - 1; ++i) {
    int id1 = local_ids[i], id2 = local_ids[i + 1];
    length += gflow->getDistance(x(id1), x(id2));
  }
  // Return the length
  return length;
}

void Group::set(Group &group) {
  *this = group;
}

void Group::add(int id) {
  global_ids.push_back(id);
  local_ids.push_back(id);
}

RealType Group::distance(const RealType *point) {
  if (size() == 0) {
    return -1.;
  }
  // Get positions
  auto x = sim_data->X();
  RealType _distance(1000000.), minimum_distance(1000000.);
  // Find closest object in the group to the point.
  for (auto id : local_ids) {
    _distance = sim_data->getGFlow()->getDistance(x(id), point);
    if (_distance < minimum_distance) {
      minimum_distance = _distance;
    }
  }
  return _distance;
}

void Group::findNetForce(RealType *frc) const {
  int sim_dimensions = sim_data->getSimDimensions();
  zeroVec(frc, sim_dimensions);
  if (size() == 0) {
    return;
  }
  // Zero vector
  zeroVec(frc, sim_dimensions);
  // Force array
  auto f = sim_data->F();
  // Compute net force
  for (auto id : local_ids) {
    plusEqVec(frc, f(id), sim_dimensions);
  }
}

RealType Group::findTotalMass() const {
  if (size() == 0) {
    return 0;
  }
  // Total the mass
  RealType m = 0, im = 0;;
  for (auto id : local_ids) {
    im = sim_data->Im(id);
    if (im == 0) {
      return 0;
    }
    m += 1. / im;
  }
  // Return the mass.
  return m;
}

void Group::findClosestObject(const RealType *point, RealType *displacement) const {
  if (size() == 0) {
    return;
  }
  // Get dimensionality.
  int sim_dimensions = sim_data->getSimDimensions();
  auto x = sim_data->X();
  Vec disp(sim_dimensions), maxDisp(sim_dimensions);
  RealType maxD = 0;
  // Find closest object in the group to the point.
  for (auto id : local_ids) {
    sim_data->getGFlow()->getDisplacement(x(id), point, disp.data);
    RealType d = sqr(disp);
    if (d > maxD) {
      maxD = d;
      maxDisp = disp;
    }
  }
  copyVec(maxDisp, displacement);
}

void Group::update_local_ids() const {
  // Make sure sizes are the same
  int _size = size(), _count = 0;
  // Update local ids
  for (int i = 0; i < _size; ++i) {
    int gid = global_ids[i];
    int lid = sim_data->getLocalID(gid);
    if (-1 < lid) {
      // Potentially shrink global id list by moving up the ids.
      local_ids[_count] = lid;
      global_ids[_count] = global_ids[i];
      ++_count;
    }
  }
  global_ids.resize(_count);
  local_ids.resize(_count);
}

void Group::shift_global_ids(const int shift) const {
  for (auto &id : global_ids) {
    id += shift;
  }
}
