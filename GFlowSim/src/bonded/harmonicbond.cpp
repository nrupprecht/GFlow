#include "harmonicbond.hpp"
// Other files
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  HarmonicBond::HarmonicBond(GFlow *gflow) : Bond(gflow), springConstant(DEFAULT_SPRING_CONSTANT) {};

  HarmonicBond::HarmonicBond(GFlow *gflow, RealType K) : Bond(gflow), springConstant(K) {};

  void HarmonicBond::addBond(int gid1, int gid2) {
    SimData *sd = Base::simData;
    int id1 = sd->getLocalID(gid1), id2 = sd->getLocalID(gid2);
    // Calculate the distance between the particles
    RealType dX[8]; // <-- Assumes that (sim_dimensions < 9)
    Base::gflow->getDisplacement(sd->X(id1), sd->X(id2), dX);
    RealType r = magnitudeVec(dX, sim_dimensions);
    // Add global ids
    gleft.push_back(gid1);
    gright.push_back(gid2);
    // Add local ids
    left.push_back(id1);
    right.push_back(id2);
    // Add equilibrium distance
    distance.push_back(r);
  }

  void HarmonicBond::post_forces() {
    // Make sure sizes are the same
    if (!checkBondVectors()) throw UnequalBondVectors();
    int nbonds = left.size();

    // Check if local ids need updating.
    if (simData->getNeedsRemake()) updateLocalIDs();

    // Get simdata, check if the local ids need updating
    RealType **x = simData->X();
    RealType **f = simData->F();
    RealType dX[8]; // <-- Assumes that (sim_dimensions < 9)

    for (int i=0; i<nbonds; ++i) {
      // Get the global, then local ids of the particles.
      int id1 = left[i], id2 = right[i];

      // Calculate displacement
      Base::gflow->getDisplacement(x[id1], x[id2], dX);
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

  void HarmonicBond::updateLocalIDs() {
    // Make sure sizes are the same
    if (!checkBondVectors()) throw UnequalBondVectors();
    int nbonds = left.size();

    SimData *sd = Base::simData;
    for (int i=0; i<nbonds; ++i) {
      int gid1 = gleft[i], gid2 = gright[i];
      int id1 = sd->getLocalID(gid1), id2 = sd->getLocalID(gid2);
      left[i] = id1; right[i] = id2;
    }
  }

}