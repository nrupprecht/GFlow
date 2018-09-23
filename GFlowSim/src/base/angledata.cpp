#include "angledata.hpp"
// Other files
#include "simdata.hpp"
#include "../utility/vectormath.hpp"
#include <cmath>

namespace GFlowSimulation {

  AngleData::AngleData(GFlow *gflow) : Base(gflow) {};

  void AngleData::addAngle(int id1, int id2, int id3, RealType str, RealType ang) {
    triples.push_back(id1);
    triples.push_back(id2);
    triples.push_back(id3);
    strength.push_back(str);
    angle.push_back(ang);
  }

  void AngleData::addAngle(int id1, int id2, int id3, RealType str) {
    // Use the current angle as the relaxd angle
    RealType x12[DIMENSIONS], x32[DIMENSIONS];
    // Get the displacements
    getDisplacement(simData->X(id1), simData->X(id2), x12, Base::gflow->getBounds(), Base::gflow->getBCs());
    getDisplacement(simData->X(id3), simData->X(id2), x12, Base::gflow->getBounds(), Base::gflow->getBCs());
    // Normalize
    normalizeVec(x12);
    normalizeVec(x32);
    // Calculate angle
    RealType dot = dotVec(x12, x32);
    RealType ang = acos(dot);
    // Add bond
    addAngle(id1, id2, str, ang);
  }

  void AngleData::post_forces() {
    // Calculate all the angle forces
    const Bounds bounds = Base::gflow->getBounds();
    const BCFlag *bcs = Base::gflow->getBCs();
    RealType x12[DIMENSIONS], x32[DIMENSIONS];
    for (int i=0; i<triples.size(); i+=3) {
      // Get indices
      int id1 = triples[i], id2 = triples[i+1], id3 = triples[i+2];
      // Get the displacements
      getDisplacement(simData->X(id1), simData->X(id2), x12, Base::gflow->getBounds(), Base::gflow->getBCs());
      getDisplacement(simData->X(id3), simData->X(id2), x12, Base::gflow->getBounds(), Base::gflow->getBCs());
      // Normalize
      normalizeVec(x12);
      normalizeVec(x32);
      // Calculate angle
      RealType dot = dotVec(x12, x32);
      RealType ang = acos(dot);
      // Compute torque
      RealType torque = strength[i/3]*(angle[i/3] - ang);
      
    }
  }

}