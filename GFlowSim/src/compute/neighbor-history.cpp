#include "neighbor-history.hpp"

namespace GFlowSimulation {

  bool NeighborHistory::in_contact(int id0, int id1) const {
    if (id1<id0) std::swap(id0, id1);

    return false;
  }

}