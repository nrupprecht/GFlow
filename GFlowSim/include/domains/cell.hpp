#ifndef __CELL_HPP__GFLOW__
#define __CELL_HPP__GFLOW__

#include "../utility/utility.hpp"

namespace GFlowSimulation {

  inline RealType min(RealType *array, int length) {
    if (length==0) return 0;
    RealType m = array[0];
    for (int i=0; i<length; ++i)
      if (array[i]<m) m = array[i];
    return m;
  }

  struct Cell {
    //! \brief Dimension setting constructor
    Cell() {};

    //! \brief Copy constructor
    Cell(const Cell& cell);

    //! \brief Operator equals.
    Cell& operator=(const Cell& cell) = default;

    int size() const;

    //! \brief The ids of the particles contained in the cell.
    vector<int> particle_ids;

    //! \brief Pointers to cells adjacent to the cell.
    vector<Cell*> adjacent;
  };

}
#endif // __CELL_HPP__GFLOW__