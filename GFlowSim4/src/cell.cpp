#include "cell.hpp"

#include "cmath"

namespace GFlowSimulation {

  Cell::Cell() : cellType(CellType::Unassigned), is_boundary_cell(false), x(nullptr), dx(nullptr), f(nullptr), sg(nullptr) {};

  Cell::~Cell() {}

  void Cell::clear() {
    id_list.clear() ;
  }

  void Cell::add(int id) {
    id_list.push_back(id) ;
  }
}