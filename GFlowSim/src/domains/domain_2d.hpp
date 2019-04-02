#ifndef __DOMAIN_2D_HPP__GFLOW__
#define __DOMAIN_2D_HPP__GFLOW__

#include "../base/domainbase.hpp"

namespace GFlowSimulation {

  class Domain2D : public DomainBase {
  public:
    //! \brief Default constructor.
    Domain2D(GFlow*);

    //! \brief Set up the domain.
    virtual void initialize() override;

    //! \brief Get all the particles within a radius of another particle
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.\n
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, RealType, vector<int>&) override;

    //! \brief Remove particles that are overlapping by more than some fraction.
    virtual void removeOverlapping(RealType) override;

    // --- Mutators

    //! \brief Set the skin depth. This function is virtual, as the inheriting class
    //! may need to remake itself after doing this.
    virtual void setSkinDepth(RealType) override;

    //! \brief Set the cell size. 
    //!
    //! Really, this suggests a cell size. It must be larger than the minimum possible cell size, 
    //! and must evenly divide the size of the domain.
    virtual void setCellSize(RealType) override;

    //! \brief Remakes interactionhandlers (if neccessary).
    //!
    //! This function should be overridden by each child to remake the interaction handlers of the forces as they
    //! see fit. It should, however, be called first in each child function (DomainBase::construct()). This function
    //! will only be called if the handlers need to be constructed.
    //!
    //! This function:
    //! (1) Wraps particle positions to their cannonical form.
    //! (2) Sets the lastUpdate timer
    //! (3) Increments the number_of_remakes counter
    //! (5) Calls the ForceMaster clear() function
    //! (6) Calls fillXVL() to record the positions of the particles
    virtual void construct() override;
 
  private:

    //! \brief Calculates the dimensions of the grid, and the cell widths and inverse widths.
    inline void calculate_domain_cell_dimensions();

    //! \brief Fill the cells with particle ids.
    inline void make_cells();

    //! \brief Adjustments for halo or ghost cells on the min side.
    int min_side_edge_cells[2];

    //! \brief Adjustments for halo or ghost cells on the max side.
    int max_side_edge_cells[2];

    //! \brief Whether there are is halo duplication in the 0, 1 directions
    bool halo_cells[2];

    //! \brief The i-th entry points to the start of the i-th cell in the linked cells array, or to -1 if the cell is empty.
    int *cell_pointers = nullptr;

    //! \brief The number of cells
    int cells_size = 0;

    //! \brief The linked cells - the i-th entry points to the next particle in the same cell as particle with local id i,
    //! or to -1 if it is the last particle in the cell.
    int *linked_cells = nullptr;

    //! \brief The size of the linked cells array, i.e. the number of particles stored in the linked cells.
    int list_size = 0;

  };

}
#endif // __DOMAIN_2D_HPP__GFLOW__