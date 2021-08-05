#include <domains/cell.hpp>

using namespace GFlowSimulation;

// Copy constructor
Cell::Cell(const Cell &cell) {
  particle_ids = cell.particle_ids;
  adjacent = cell.adjacent;
};

int Cell::size() const {
  return particle_ids.size();
}
