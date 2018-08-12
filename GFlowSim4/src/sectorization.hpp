#ifndef __SECTORIZATION_HPP__GFLOW__
#define __SECTORIZATION_HPP__GFLOW__

#include "domainbase.hpp"
#include "array.hpp"

namespace GFlowSimulation {

  class Sectorization : public DomainBase {
  public:
    //! Constructor
    Sectorization(GFlow *);

    //! @brief Pre-integrate calls sectorize
    virtual void pre_integrate();

    //! @brief Doesn't do anything, but must be overloaded.
    //!
    //! In theory, this exchanges particles between processors. Sectorization is 
    //! designed for single processor jobs, so this function does not do anything.
    //! Since it is purely virtual for DomainBase, it must be overloaded.
    virtual void exchange_particles() {};

    // --- Accessors

    //! @brief Get a sector
    const vector<int>& getSector(int *) const;

    // --- Locator functions

    //! @brief Get all the particles within a radius of another particle
    //!
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.\n
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, RealType, vector<int>&);

    // --- Mutators

    //! @brief Set the skin depth, and remake the sectors
    virtual void setSkinDepth(RealType);

    //! @brief Set the cell size. 
    //!
    //! Really, this suggests a cell size. It must be larger than the minimum possible cell size, 
    //! and must evenly divide the size of the domain.
    virtual void setCellSize(RealType);

    //! @brief Set the cutoff factor, and remake the sectors
    virtual void setCutoffFactor(RealType);
    
    // GFlow is a friend class
    friend class GFlow;
    friend class SectorizationData;
    
  private:
    // --- Helper functions

    //! @brief Remake the verlet lists
    virtual void construct();

    //! @brief Create sectors
    inline void makeSectors();

    //! @brief Create the sectors
    inline void create_cells();

    //! @brief Create verlet lists for all forces
    inline void makeVerletLists();
    
    //! @brief The id's of particles in each sector - these are the actual sectors
    Array< vector<int> > sectors;

    //! @brief Current head - for verlet list creation
    int currentHead;
  };

};

#endif // __SECTORIZATION_HPP__GFLOW__