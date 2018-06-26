#include "force.hpp"

namespace GFlowSimulation {

  Force::Force(GFlow *gflow) : Base(gflow), typeMap(nullptr) {};

  Force::~Force() {
    if (typeMap)    delete [] typeMap;
  }

  int Force::lastHead() {
    return verletList.lastHead();
  }

  int Force::vlSize() {
    return verletList.vlSize();
  }

  int Force::vlHSize() {
    return verletList.vlHSize();
  }

  const VerletList& Force::getVerletList() {
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