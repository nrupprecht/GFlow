#include "group.hpp"
// Other files
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  Group::Group(const Group &g) {
    global_ids = g.global_ids;
    local_ids = g.local_ids;
  }

  Group& Group::operator=(const Group &g) {
    global_ids = g.global_ids;
    local_ids = g.local_ids;
    // Return
    return *this;
  }

  int Group::size() const {
    return global_ids.size();
  }

  void Group::add(int id) {
    global_ids.push_back(id);
    local_ids.push_back(id);
  }

  int Group::at(int i) const {
    return local_ids.at(i);
  }
    
  int Group::g_at(int i) const {
    return global_ids.at(i);
  }

  void Group::findCenterOfMass(RealType *rcm, SimData *simData) const {
    if (size()==0) return;
    // Get the dimensionaliry
    int sim_dimensions = simData->getSimDimensions();
    // Dispacement vector
    RealType *dr = new RealType[sim_dimensions], *r0 = new RealType[sim_dimensions];
    zeroVec(rcm, sim_dimensions);

    // Get the bounds and boundary conditions
    GFlow *gflow = simData->getGFlow();
    Bounds bounds = gflow->getBounds(); // Simulation bounds
    BCFlag *boundaryConditions = new BCFlag[sim_dimensions];
    copyVec(gflow->getBCs(), boundaryConditions, sim_dimensions); // Keep a local copy of the bcs

    // Get position, mass arrays
    RealType **x = simData->X();
    RealType *im = simData->Im();
    RealType mass = 0;

    // Set reference point
    copyVec(x[local_ids[0]], r0, sim_dimensions);
    // Go through all particles
    for (auto id : local_ids) {
      // Get displacement between the particle and the reference point.
      gflow->getDisplacement(x[id], r0, dr);
      // Get the mass of the particle - assumes none of the particles have infinite mass.
      RealType m = 1./im[id];
      mass += m;
      // Update rcm
      plusEqVecScaled(rcm, dr, m, sim_dimensions);
    }
    // Divide by total mass and add to reference position to get the actual com position.
    scalarMultVec(1./mass, rcm, sim_dimensions);
    // Add to reference point, and do wrapping
    plusEqVec(rcm, r0, sim_dimensions);

    // Harmonic boundary conditions.
    for (int d=0; d<sim_dimensions; ++d) {
      if (rcm[d]<bounds.min[d]) rcm[d] += bounds.wd(d);
      else if (bounds.max[d]<rcm[d]) rcm[d] -= bounds.wd(d);
    }

    // Clean up
    delete [] boundaryConditions;
    delete [] dr;
    delete [] r0;
  }

  void Group::findCOMVelocity(RealType *v, SimData *simData) const {
    if (size()==0) return;
    // Get the dimensionaliry
    int sim_dimensions = simData->getSimDimensions();
    // Zero vector
    zeroVec(v, sim_dimensions);
    // Compute new momentum
    RealType mass = 0;
    for (auto id : local_ids) {
      RealType m = 1./simData->Im(id);
      plusEqVecScaled(v, simData->V(id), m, sim_dimensions);
      mass += m;
    }
    // Normalize
    scalarMultVec(1./mass, v, sim_dimensions);
  }

  void Group::findNetForce(RealType *f, SimData *simData) const {
    if (size()==0) return;
    // Get the dimensionaliry
    int sim_dimensions = simData->getSimDimensions();
    // Zero vector
    zeroVec(f, sim_dimensions);
    // Force array
    RealType **F = simData->F();
    // Compute net force
    for (auto id : local_ids) {
      plusEqVec(f, F[id], sim_dimensions);
    }
  }

  void Group::update_local_ids(SimData *simData) {
    // Make sure sizes are the same
    int _size = size();
    // Update local ids
    for (int i=0; i<_size; ++i) {
      int gid = global_ids[i];;
      int lid = simData->getLocalID(gid);
      local_ids[i] = lid;
    }
  }

}