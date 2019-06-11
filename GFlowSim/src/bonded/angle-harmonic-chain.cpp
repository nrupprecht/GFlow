#include "angle-harmonic-chain.hpp"

namespace GFlowSimulation {

  AngleHarmonicChain::AngleHarmonicChain(GFlow *gflow) 
    : Bonded(gflow), springConstant(DEFAULT_SPRING_CONSTANT), angleConstant(0.005) {};

  AngleHarmonicChain::AngleHarmonicChain(GFlow *gflow, RealType K) 
    : Bonded(gflow), springConstant(K), angleConstant(0.005) {};

  void AngleHarmonicChain::addAtom(int gid) {
    SimData *sd = Base::simData;
    if (global_ids.empty()) {
      local_ids.push_back(gid);
      global_ids.push_back(gid);
    }
    else {
      int gid1 = *global_ids.rbegin();
      int id1 = sd->getLocalID(gid1), id = sd->getLocalID(gid);
      // Calculate the distance between the particles
      RealType dX[8]; // <-- Assumes that (sim_dimensions < 9)
      Base::gflow->getDisplacement(sd->X(id1), sd->X(id), dX);
      RealType r = magnitudeVec(dX, sim_dimensions);
      // Add global ids
      global_ids.push_back(gid);
      // Add local ids
      local_ids.push_back(id);
      // Add equilibrium distance
      distance.push_back(r);
    }
  }

  int AngleHarmonicChain::size() const {
    return global_ids.size();
  }

  void AngleHarmonicChain::updateLocalIDs() const {
    // Make sure sizes are the same
    int nbonds = global_ids.size();
    // Update local ids
    SimData *sd = Base::simData;
    for (int i=0; i<nbonds; ++i) {
      int gid = global_ids[i];
      int id = sd->getLocalID(gid);
      local_ids[i] = id;
    }
  }

}