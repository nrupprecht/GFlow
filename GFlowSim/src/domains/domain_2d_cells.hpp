#ifndef __DOMAIN_2D_CELLS_HPP__GFLOW__
#define __DOMAIN_2D_CELLS_HPP__GFLOW__

#include "../base/domainbase.hpp"
#include "cell.hpp"

namespace GFlowSimulation {

  class Domain2DCells : public DomainBase {
  public:
    //! \brief Default constructor.
    Domain2DCells(GFlow*);

    //! \brief Set up the domain.
    virtual void initialize() override;

    //! \brief Get all the particles within a radius of another particle
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.\n
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, vector<int>&, RealType=-1.) override;

    //! \brief Get all the particles withing a radius of some position.
    virtual void getAllWithin(Vec, vector<int>&, RealType=-1.) override {};

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

    //! \brief Construct interactions for a single particle
    //!
    //! If the insert flag is set to true, the particle will first be added to the relevent domain or data structure
    //! (if applicable), if not, it will not be added.
    //! This function should be called after construct.
    virtual void constructFor(int, bool=false) override {
      cout << "Warning: constructFor not implemented for this class." << endl;
      exit(0);
    }
 
  private:
    //! \brief Make halo particles.
    //!
    //! This function assumes that there are currently no halo or ghost particles stored in simdata.
    void construct_halo_particles();

    //! \brief Calculates the dimensions of the grid, and the cell widths and inverse widths.
    virtual void calculate_domain_cell_dimensions() override;

    //! \brief Create the cells.
    inline void create_cells();

    //! \brief Fill the cells with particle ids.
    inline void update_cells();

    //! \brief Go through a cell, starting with a certain particle, checking all the particles for interactions with a specified particle.
    inline void check_cell(int, const Cell*);

    //! \brief Do a cell check for a large particle.
    inline void check_cell_large(int, int);

    //! \brief Adjustments for halo or ghost cells on the min side.
    int min_side_edge_cells[2];

    //! \brief Adjustments for halo or ghost cells on the max side.
    int max_side_edge_cells[2];

    //! \brief A vector holding all the cells in the domain.
    vector<Cell> cells;
  };

}
#endif // __DOMAIN_2D_CELLS_HPP__GFLOW__