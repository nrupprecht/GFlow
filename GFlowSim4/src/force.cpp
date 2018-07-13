#include "force.hpp"
// Other files
#include "simdata.hpp"
#include "vectormath.hpp"

namespace GFlowSimulation {

  Force::Force(GFlow *gflow) : Base(gflow) {};

  Force::~Force() {}

  int Force::lastHead() const {
    return verletList.lastHead();
  }

  int Force::vlSize() const {
    return verletList.vlSize();
  }

  int Force::vlHSize() const {
    return verletList.vlHSize();
  }

  const VerletList& Force::getVerletList() const {
    return verletList;
  }

  void Force::clearVerletList() {
    verletList.clear();
  }

  void Force::addVerletPair(int id1, int id2) {
    // Add the head if it is new
    if (id1!=lastHead()) verletList.addHead(id1);
    // Add id2 to the head
    verletList.addToHead(id2);
  }

}