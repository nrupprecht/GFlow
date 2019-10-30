#ifndef __DOMAIN_LIST_HPP__GFLOW__
#define __DOMAIN_LIST_HPP__GFLOW__

#include "../base/domainbase.hpp" 
#include "../utility/generic-dimension.hpp"

namespace GFlowSimulation {

  /**
  * \brief A true linked cells domain object, where the cells are implemented via a linked list.
  * 
  */
  class DomainList : public DomainBase {
  public:
    //! \brief Default constructor.
    DomainList(GFlow*);

    //! \brief Get all the particles within a radius of another particle
    //!
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, vector<int>&, RealType=-1.) override;

    //! \brief Get all the particles withing a radius of some position.
    virtual void getAllWithin(Vec, vector<int>&, RealType=-1.) override;

    //! \brief Construct interactions for a single particle.
    //!
    //! If the insert flag is set to true, the particle will first be added to the relevent domain or data structure
    //! (if applicable), if not, it will not be added.
    //! This function should be called after construct.
    virtual void constructFor(int, bool=false) override;

    //! \brief Traverse the cells structure, calling a function on pairs of particles that are close to one another.
    //!
    //! The function that is passed in should expect to receive particles' id1, id2, wrapping type (0 - no wrapping
    //! required, 1 - wrapping required), radius of particle 1, radius of particle 2, and distance between particles.
    virtual void traversePairs(PairFunction) override;

    //! \brief Function that traverses ghost particles, calling a function on all ghost-real pairs of particls that are within
    //! cutoff + skin depth distance of one another.
    //!
    //! The function (a PairFunction) that is passed in should expect to receive particles' id1, id2, wrapping type (0 - no wrapping
    //! required, 1 - wrapping required), radius of particle 1, radius of particle 2, and distance between particles.
    virtual void traverseGhostPairs(PairFunction) override {};

  private:

    // --- Overloaded functions

    //! \brief Update the linked cells structure.
    virtual void structure_updates() override;

    //! \brief Calculates the domain cell dimensions, widths, and inverse widths given 
    //! that the cutoff has been calculated. 
    virtual void calculate_domain_cell_dimensions() override;

    //! \brief Create the cells.
    virtual void create_cells() override;

    // --- Data members

    //! \brief The N-th entry points to the first particle in the N-th cell, or to -1 if the cell is empty.
    vector<int> cell_list;
    //! \brief The particle points to the next particle in the same cell, or to -1 if it is the last particle in a cell.
    vector<int> particle_list;

  };

}
#endif // __DOMAIN_LIST_HPP__GFLOW__