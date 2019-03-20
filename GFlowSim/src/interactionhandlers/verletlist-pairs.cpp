#include "verletlist-pairs.hpp"
// Other files
#include "../base/interaction.hpp"

namespace GFlowSimulation {

  VerletListPairs::VerletListPairs(GFlow *gflow) : InteractionHandler(gflow) {};

  void VerletListPairs::addPair(const int id1, const int id2) {
    verlet.push_back(id1);
    verlet.push_back(id2);
  }

  void VerletListPairs::clear() {
    verlet.clear();
  }

  int VerletListPairs::size() const {
    return verlet.size();
  }

}
