#ifndef __SECTORIZATION_HPP__GFLOW__
#define __SECTORIZATION_HPP__GFLOW__

#include "gflow.hpp"
#include "array.hpp"

namespace GFlowSimulation {

  class Sectorization : public Base {
  public:
    // Constructor
    Sectorization(GFlow *);

    // Destructor
    ~Sectorization();

    // Pre-integrate calls sectorize
    virtual void pre_integrate();

    // Pre-forces is where resectorization may happen
    virtual void pre_forces();

    // Sectorize the particles
    virtual void sectorize();

    // --- Accessors

    // Get the total number of sectors
    int getNumSectors();

    // Get the skin depth
    RealType getSkinDepth();

    // Get the cutoff
    RealType getCutoff();

    // --- Mutators

    // Set the skin depth
    void setSkinDepth(RealType);
    
    // GFlow is a friend class
    friend class GFlow;
    
  private:
    // --- Helper functions

    // Create sectors
    inline void makeSectors();

    // Check whether we need to resectorize, if so, then do it
    inline void checkSectors();

    // Delete and set as null the xVL array
    inline void nullXVL();

    // Set up xVL array
    inline void setupXVL(int);

    // Bounds - get these from gflow
    Bounds bounds;

    // Number of sectors in each dimension
    int dims[DIMENSIONS];
    
    // The id's of particles in each sector
    Array< vector<int> > sectors;

    // The widths of sectors in each dimension
    RealType widths[DIMENSIONS];
    // The inverse widths of sectors in each dimension
    RealType inverseW[DIMENSIONS];

    // The positions of the particles at the last verlet list creation
    RealType **xVL;
    int sizeXVL; // The size of the xVL array

    // Sectorization constants
    RealType skinDepth, cutoff;
  };

};

#endif // __SECTORIZATION_HPP__GFLOW__