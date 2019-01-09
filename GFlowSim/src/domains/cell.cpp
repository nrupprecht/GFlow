#include "cell.hpp"

namespace GFlowSimulation {

  // Dimension setting constructor
  Cell::Cell(int d) {}

  // Copy constructor
  Cell::Cell(const Cell& cell) {
    particle_ids = cell.particle_ids;
    adjacent = cell.adjacent;
  };

  Cell Cell::operator=(const Cell& cell) {
    // Set data
    particle_ids = cell.particle_ids;
    adjacent = cell.adjacent;
    // Return
    return *this;
  }

  int Cell::size() const {
    return particle_ids.size();
  }

}