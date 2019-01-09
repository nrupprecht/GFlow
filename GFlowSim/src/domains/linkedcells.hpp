#ifndef __LINKED_CELLS_HPP__GFLOW__
#define __LINKED_CELLS_HPP__GFLOW__

#include "../base/domainbase.hpp"

namespace GFlowSimulation {

  class LinkedCells : public DomainBase {
  public:
    //! \brief Default constructor.
    LinkedCells(GFlow*);

    virtual void initialize() override;

    //! \brief Exchange particles between processors
    virtual void exchange_particles() override;

    //! \brief Get all the particles within a radius of another particle
    virtual void getAllWithin(int, RealType, vector<int>&) override;

    virtual void removeOverlapping(RealType) override;

    virtual void construct() override;

    virtual void setCellSize(RealType) override;

  private:

    //! \brief Creates the linked cell arrays.
    inline void create_cells();

    inline void update_cells();

    //! \brief Check a cell for the neighbors of a particle.
    inline void check_cell(int, int);

    inline bool correct_index(int&, int&);

    //! \brief The head part of the linked cells.
    int *heads;
    //! \brief The number of cells. The size of the head array.
    int ncells;

    //! \brief The list part of the linked cells.
    int *list;
    //! \brief The length of the list array.
    int nlist;

  };

}
#endif // __LINKED_CELLS_HPP__GFLOW__