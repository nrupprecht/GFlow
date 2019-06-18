#include "angle-harmonic-chain.hpp"

namespace GFlowSimulation {

  AngleHarmonicChain::AngleHarmonicChain(GFlow *gflow) 
    : Bonded(gflow), springConstant(DEFAULT_SPRING_CONSTANT), angleConstant(0.005) {
    use_correspondence = true;
  };

  AngleHarmonicChain::AngleHarmonicChain(GFlow *gflow, RealType K) 
    : Bonded(gflow), springConstant(K), angleConstant(0.005) {
    use_correspondence = true;
  };

  AngleHarmonicChain::AngleHarmonicChain(GFlow *gflow, RealType K, RealType A)
    : Bonded(gflow), springConstant(K), angleConstant(A) {
    use_correspondence = true;
  };

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
    forces.push_back(Vec(sim_dimensions));
  }

  int AngleHarmonicChain::size() const {
    return Group::size();
  }

  void AngleHarmonicChain::interact() const {
    // Sets potential, virials to zero
    Bonded::interact();
    // Zero force buffers.
    for (auto &v : forces) v.zero();
    // Possibly update local ids
    if (simData->getNeedsRemake()) update_local_ids(simData);
  }

  void AngleHarmonicChain::setSpringConstant(RealType s) {
    springConstant = s;
  }

  void AngleHarmonicChain::setAngleConstant(RealType s) {
    angleConstant = s;
  }

  const Vec& AngleHarmonicChain::getForce(int i) {
    if (i<0 || forces.size()<=i) throw Exception("Angle harmonic chain get force out of range.");
    return forces[i];
  }

  const Vec& AngleHarmonicChain::getForceByID(int id) {
    if (!use_correspondence) throw Exception("Angle harmonic chain must have use_correspondence = true.");
    auto it = correspondence.find(id);
    if (it==correspondence.end()) throw Exception("Angle harmonic chain does not contain that local id.");
    return forces[it->second];
  }

}