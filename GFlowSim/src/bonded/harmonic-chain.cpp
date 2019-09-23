#include "harmonic-chain.hpp"

namespace GFlowSimulation {

  HarmonicChain::HarmonicChain(GFlow *gflow) : Bonded(gflow), springConstant(DEFAULT_SPRING_CONSTANT) {};

  HarmonicChain::HarmonicChain(GFlow *gflow, RealType K) : Bonded(gflow), springConstant(K) {};

  void HarmonicChain::addAtom(int gid) {
    if (global_ids.empty()) {
      local_ids.push_back(gid);
      global_ids.push_back(gid);
    }
    else {
      int gid1 = *global_ids.rbegin();
      int id1 = simData->getLocalID(gid1), id = simData->getLocalID(gid);
      // Calculate the distance between the particles
      RealType dX[8]; // <-- Assumes that (sim_dimensions < 9)
      Base::gflow->getDisplacement(simData->X(id1), simData->X(id), dX);
      RealType r = magnitudeVec(dX, sim_dimensions);
      // Add global ids
      global_ids.push_back(gid);
      // Add local ids
      local_ids.push_back(id);
      // Add equilibrium distance
      distance.push_back(r);
    }
  }

  int HarmonicChain::size() const {
    return global_ids.size();
  }

  void HarmonicChain::interact() const {
    // Get simdata, check if the local ids need updating
    auto f = simData->F();
    if (simData->getNeedsRemake()) updateLocalIDs();
    RealType dX[8]; // <-- Assumes that (sim_dimensions < 9)
    
    for (int i=0; i+1<local_ids.size(); ++i) {
      // Get the global, then local ids of the particles.
      int id1 = local_ids[i], id2 = local_ids[i+1];

      // Calculate displacement
      Base::gflow->getDisplacement(simData->X(id1), simData->X(id2), dX);
      RealType r = magnitudeVec(dX, sim_dimensions);

      // Calculate displacement from equilibrium
      RealType dr = r - distance[i];
      // Makes nX the unit vector of dX
      normalizeVec(dX, sim_dimensions);
      // nX is now the force vector
      scalarMultVec(springConstant*dr, dX, sim_dimensions);
      // Add forces to particles
      minusEqVec(f[id1], dX, sim_dimensions);
      plusEqVec (f[id2], dX, sim_dimensions);
    }
  }

  void HarmonicChain::updateLocalIDs() const {
    // Make sure sizes are the same
    int nbonds = global_ids.size();
    // Update local ids
    for (int i=0; i<nbonds; ++i) {
      int gid = global_ids[i];
      int id = simData->getLocalID(gid);
      local_ids[i] = id;
    }
  }

}