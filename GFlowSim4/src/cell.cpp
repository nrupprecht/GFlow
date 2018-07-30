#include "cell.hpp"

#include "cmath"

namespace GFlowSimulation {

  Cell::Cell() : adjacent_cell_id(nullptr), id_list(nullptr), adjacent_cell_id_size(0), id_list_size(0), 
    cellType(CellType::Unassigned), is_boundary_cell(false), id_list_capacity(0) {};

  Cell::~Cell() {
    if (adjacent_cell_id) delete [] adjacent_cell_id;
    if (id_list)          delete [] id_list;
  }

  void Cell::clear() {
    id_list_size = 0;
  }

  void Cell::add(int id) {
    if (id_list_size+1>=id_list_capacity) {
      // Resize
      int new_capacity = ceil( 1.5*(id_list_size+1) );
      int *new_id_list = new int[new_capacity];
      // Copy elements
      for (int i=0; i<id_list_size; ++i) new_id_list[i] = id_list[i];
      // Delete old
      if (id_list) delete [] id_list;
      // Set new
      id_list = new_id_list;
      id_list_capacity = new_capacity;
    }
    // Add new element (we don't check to make sure that this element is 
    // not already here - this should never occur)
    id_list[id_list_size++] = id;
  }
}