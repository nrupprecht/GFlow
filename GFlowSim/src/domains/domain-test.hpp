#ifndef __DOMAIN_TEST_HPP__GFLOW__
#define __DOMAIN_TEST_HPP__GFLOW__

#include "../base/domainbase.hpp"

namespace GFlowSimulation {

  class Domain_Test : public DomainBase {
  public:
    Domain_Test(GFlow*);

    //! \brief Construct interactions for a single particle.
    //!
    //! If the insert flag is set to true, the particle will first be added to the relevent domain or data structure
    //! (if applicable), if not, it will not be added.
    //! This function should be called after construct.
    virtual void constructFor(int, bool=false) override;

    // --- Locator functions

    //! \brief Get all the particles within a radius of another particle
    //!
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.\n
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, vector<int>&, RealType=-1.) override;

    //! \brief Get all the particles withing a radius of some position.
    //!
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(Vec, vector<int>&, RealType=-1.) override;

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
    virtual void traverseGhostPairs(PairFunction) override;

  private:

    //! \brief Calculates the domain cell dimensions, widths, and inverse widths given 
    //! that the cutoff has been calculated.
    virtual void calculate_domain_cell_dimensions() override;

    //! \brief Create the cells.
    virtual void create_cells() override;

    vector<int> cell_heads;
    vector<int> next_in_cell;

  };

}
#endif // __DOMAIN_TEST_HPP__GFLOW__