#ifndef __DOMAIN_2D_CELLS_HPP__GFLOW__
#define __DOMAIN_2D_CELLS_HPP__GFLOW__

#include "../base/domainbase.hpp"
#include "cell.hpp"

namespace GFlowSimulation {

  class Domain2DCells : public DomainBase {
  public:
    //! \brief Default constructor.
    Domain2DCells(GFlow*);

    //! \brief Get all the particles within a radius of another particle
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.\n
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, vector<int>&, RealType=-1.) override { cout << "Error: unimplemented.\n"; }

    //! \brief Get all the particles withing a radius of some position.
    virtual void getAllWithin(Vec, vector<int>&, RealType=-1.) override { cout << "Error: unimplemented.\n"; }

    // --- Mutators

    //! \brief Construct interactions for a single particle
    //!
    //! If the insert flag is set to true, the particle will first be added to the relevent domain or data structure
    //! (if applicable), if not, it will not be added.
    //! This function should be called after construct.
    virtual void constructFor(int, bool=false) override { cout << "Error: unimplemented.\n"; }

    //! \brief Function that traverses the internal data structures of interaction handlers and calls a function on
    //! all pairs of particles that are within cutoff + skin depth distance of one another.
    //!
    //! The function (a PairFunction) that is passed in should expect to receive particles' id1, id2, wrapping type (0 - no wrapping
    //! required, 1 - wrapping required), radius of particle 1, radius of particle 2, and distance between particles.
    virtual void traversePairs(PairFunction) override;

    //! \brief Function that traverses ghost particles, calling a function on all ghost-real pairs of particls that are within
    //! cutoff + skin depth distance of one another.
    //!
    //! The function (a PairFunction) that is passed in should expect to receive particles' id1, id2, wrapping type (0 - no wrapping
    //! required, 1 - wrapping required), radius of particle 1, radius of particle 2, and distance between particles.
    virtual void traverseGhostPairs(PairFunction) override {};
 
  private:
    //! \brief Make halo particles.
    //!
    //! This function assumes that there are currently no halo or ghost particles stored in simdata.
    void construct_halo_particles();

    //! \brief Update the cells structure.
    virtual void structure_updates() override;

    //! \brief Calculates the dimensions of the grid, and the cell widths and inverse widths.
    virtual void calculate_domain_cell_dimensions() override;

    //! \brief Create the cells.
    virtual void create_cells() override;

    //! \brief Fill the cells with particle ids.
    inline void update_cells();

    //! \brief Go through a cell, starting with a certain particle, checking all the particles for interactions with a specified particle.
    inline void check_cell(int, const Cell*, PairFunction);

    //! \brief Do a cell check for a large particle.
    inline void check_cell_large(int, int, PairFunction);

    //! \brief A vector holding all the cells in the domain.
    vector<Cell> cells;
  };

}
#endif // __DOMAIN_2D_CELLS_HPP__GFLOW__