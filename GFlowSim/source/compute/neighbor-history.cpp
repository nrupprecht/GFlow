#include <compute/neighbor-history.hpp>

using namespace GFlowSimulation;

int &NeighborHistory::in_contact(const int id0, const int j) {
  // See if id1 is in id0's contact vector.
  return contact_vector[id0][j];
}

real &NeighborHistory::value(const int id0, const int j) {
  // Return the history value. For the sake of speed, we don't check bounds.
  return value_vector[id0][j];
}
