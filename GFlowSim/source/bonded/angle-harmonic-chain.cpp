#include <bonded/angle-harmonic-chain.hpp>

using namespace GFlowSimulation;

AngleHarmonicChain::AngleHarmonicChain(GFlow *gflow)
    : Bonded(gflow), Group(gflow),
      springConstant(1.5 * DEFAULT_SPRING_CONSTANT),
      angleConstant(0.05) {};

AngleHarmonicChain::AngleHarmonicChain(GFlow *gflow, RealType K)
    : Bonded(gflow), Group(gflow), springConstant(K),
      angleConstant(0.05) {};

AngleHarmonicChain::AngleHarmonicChain(GFlow *gflow, RealType K, RealType A)
    : Bonded(gflow), Group(gflow),
      springConstant(K), angleConstant(A) {};

void AngleHarmonicChain::addAtom(int gid) {
  if (global_ids.empty()) {
    add(gid);
  }
  else {
    int gid1 = *global_ids.rbegin(); // Last particle in the chain
    int id1 = simData->getLocalID(gid1), id = simData->getLocalID(gid);
    // Calculate the distance between the particles
    RealType dX[8]; // <-- Assumes that (sim_dimensions < 9)
    Base::gflow->getDisplacement(simData->X(id1), simData->X(id), dX);
    RealType r = magnitudeVec(dX, sim_dimensions);
    // Add global ids
    add(gid);
    // Add equilibrium distance
    distance.push_back(r);
  }
  forces.emplace_back(sim_dimensions);
}

int AngleHarmonicChain::size() const {
  return Group::size();
}

void AngleHarmonicChain::interact() const {
  // Sets potential, virials to zero
  Bonded::interact();
  // Zero force buffers.
  for (auto &v : forces) {
    v.zero();
  }
  // Possibly update local ids
  if (simData->getNeedsRemake()) {
    update_local_ids();
  }
}

void AngleHarmonicChain::setSpringConstant(RealType s) {
  springConstant = s;
}

void AngleHarmonicChain::setAngleConstant(RealType s) {
  angleConstant = s;
}

const Vec &AngleHarmonicChain::getForce(int i) {
  if (i < 0 || forces.size() <= i) {
    throw Exception("Angle harmonic chain get force out of range.");
  }
  return forces[i];
}

const Vec &AngleHarmonicChain::getForceByID(int id) {
  int index = getIndex(id);
  if (index < 0) {
    throw Exception("Chain does not contain that particle.");
  }
  return forces[index];
}
