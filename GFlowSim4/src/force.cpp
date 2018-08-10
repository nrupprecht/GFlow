#include "force.hpp"
// Other files
#include "simdata.hpp"
#include "vectormath.hpp"

namespace GFlowSimulation {

  Force::Force(GFlow *gflow) : Base(gflow), virial(0), verletList(new VerletList(gflow)), parameters(nullptr), forcePtr(nullptr) {};

  Force::~Force() {
    if (verletList) delete verletList;
    if (parameters) delete parameters;
  }

  void Force::calculateForces() const {
    // Check if there are forces to calculate. Virial is reset.
    if (!Force::initCalculation() || forcePtr==nullptr) return; 

    //RealType param_pack[] = { strength, cutoff };
    verletList->forceKernel(forcePtr, parameters, &virial);
  }

  bool Force::initCalculation() const {
    // Reset the virial
    virial = 0;
    // Return whether we should calculate forces
    return verletList->vlSize()>0;
  }

  int Force::vlSize() const {
    return verletList->vlSize();
  }

  const VerletList& Force::getVerletList() const {
    return *verletList;
  }

  int Force::getVirial() const {
    return virial;
  }

  void Force::clearVerletList() {
    verletList->clear();
  }

  void Force::addVerletPair(int id1, int id2) {
    // Add the head if it is new
    verletList->addPair(id1, id2);
  }

}