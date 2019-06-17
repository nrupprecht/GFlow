#include "angle-harmonic-chain.hpp"

namespace GFlowSimulation {

  AngleHarmonicChain::AngleHarmonicChain(GFlow *gflow) 
    : Bonded(gflow), springConstant(DEFAULT_SPRING_CONSTANT), angleConstant(0.005) {};

  AngleHarmonicChain::AngleHarmonicChain(GFlow *gflow, RealType K) 
    : Bonded(gflow), springConstant(K), angleConstant(0.005) {};

  AngleHarmonicChain::AngleHarmonicChain(GFlow *gflow, RealType K, RealType A)
    : Bonded(gflow), springConstant(K), angleConstant(A) {};

  void AngleHarmonicChain::addAtom(int gid) {
    SimData *sd = Base::simData;
    if (global_ids.empty()) add(gid);
    else {
      int gid1 = *global_ids.rbegin(); // Last particle in the chain
      int id1 = sd->getLocalID(gid1), id = sd->getLocalID(gid);
      // Calculate the distance between the particles
      RealType dX[8]; // <-- Assumes that (sim_dimensions < 9)
      Base::gflow->getDisplacement(sd->X(id1), sd->X(id), dX);
      RealType r = magnitudeVec(dX, sim_dimensions);
      // Add global ids
      add(gid);
      // Add equilibrium distance
      distance.push_back(r);
    }
  }

  int AngleHarmonicChain::size() const {
    return Group::size();
  }

  void AngleHarmonicChain::setSpringConstant(RealType s) {
    springConstant = s;
  }

  void AngleHarmonicChain::setAngleConstant(RealType s) {
    angleConstant = s;
  }

}