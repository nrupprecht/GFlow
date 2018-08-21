#include "cell.hpp"

namespace GFlowSimulation {

  Cell::Cell() : cellType(CellType::Unassigned), is_boundary_cell(false), x(nullptr), f(nullptr), sg(nullptr), 
    mask(nullptr), capacity(0), loaded(false) {};

  Cell::~Cell() {}

  void Cell::clear() {
    id_list.clear() ;
  }

  void Cell::add(int id) {
    id_list.push_back(id) ;
  }

  int Cell::size() {
    return id_list.size();
  }

}