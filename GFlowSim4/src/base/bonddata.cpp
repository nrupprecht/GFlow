#include "bonddata.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  BondData::BondData(GFlow* gflow) : Base(gflow) {};

  void BondData::addBond(int id1, int id2, RealType str, RealType len) {
    pairs.push_back(id1);
    pairs.push_back(id2);
    strength.push_back(str);
    length.push_back(len);
  }

  void BondData::addBond(int id1, int id2, RealType str) {
    // Use the current length as the relaxed length
    RealType disp[DIMENSIONS];
    getDisplacement(simData->X(id1), simData->X(id2), disp, Base::gflow->getBounds(), Base::gflow->getBCs());
    RealType len = magnitudeVec(disp);
    addBond(id1, id2, str, len);
  }

  void BondData::post_forces() {
    // Calculate all the bond forces
    const Bounds bounds = Base::gflow->getBounds();
    const BCFlag *bcs = Base::gflow->getBCs();
    RealType force[DIMENSIONS];
    for (int i=0; i<pairs.size(); i+=2) {
      int id1 = pairs[i], id2 = pairs[i+1];
      RealType *x1 = Base::simData->x[id1], *x2 = Base::simData->x[id2];
      // Displacement
      getDisplacement(x1, x2, force, bounds, bcs); 
      // Distance
      RealType distance = magnitudeVec(force);
      // Find x-x_eq
      RealType diffX = distance - length[i/2];
      // Normalize [force]
      normalizeVec(force);
      // Find force
      scalarMultVec(strength[i/2]*diffX, force); // F = - K \vec{x-x_eq}
      // Add the forces
      minusEqVec(Base::simData->f[id1], force);
      plusEqVec (Base::simData->f[id2], force);
    }
  }

}