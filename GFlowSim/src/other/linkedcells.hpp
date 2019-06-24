#ifndef __LINKED_CELLS_HPP__GFLOW__
#define __LINKED_CELLS_HPP__GFLOW__

namespace GFlowSimulation {

  class LinkedCells {
  public:

    void add();

    //! \brief Reserve some space size.
    void resize(int);

  private:
    //! \brief Array for cell-#, local id. We then sort this.
    vector<pair<int, int> > array;

    //! \brief The total number of cells.
    int total_cells = 0;

    //! \brief The i-th entry in pointers points to the particle in particles that is in the i-th cell.
    int *pointers = nullptr;
    int *particles = nullptr;
  };

}
#endif // __LINKED_CELLS_HPP__GFLOW__