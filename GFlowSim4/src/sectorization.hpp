#ifndef __SECTORIZATION_HPP__GFLOW__
#define __SECTORIZATION_HPP__GFLOW__

#include "domainbase.hpp"
#include "array.hpp"

namespace GFlowSimulation {

  class Sectorization : public DomainBase {
  public:
    //! Constructor
    Sectorization(GFlow *);

    // Pre-integrate calls sectorize
    virtual void pre_integrate();

    //! @brief Doesn't do anything, but must be overloaded.
    //!
    //! In theory, this exchanges particles between processors. Sectorization is 
    //! designed for single processor jobs, so this function does not do anything.
    //! Since it is purely virtual for DomainBase, it must be overloaded.
    virtual void exchange_particles() {};

    // --- Accessors

    // Get a sector
    const vector<int>& getSector(int *) const;

    // --- Locator functions

    //! @brief Get all the particles within a radius of another particle
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.\n
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, RealType, vector<int>&);

    // --- Mutators

    //! Set the skin depth, and remake the sectors
    virtual void setSkinDepth(RealType);

    //! Set the cutoff factor, and remake the sectors
    virtual void setCutoffFactor(RealType);
    
    // GFlow is a friend class
    friend class GFlow;
    friend class SectorizationData;
    
  private:
    // --- Helper functions

    virtual void remake_verlet();

    // Create sectors
    inline void makeSectors();

    // Create verlet lists for all forces
    inline void makeVerletLists();
    
    // The id's of particles in each sector - these are the actual sectors
    Array< vector<int> > sectors;

    // Current head - for verlet list creation
    int currentHead;
  };

};

#endif // __SECTORIZATION_HPP__GFLOW__