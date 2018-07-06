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

    // Get a sector
    const vector<int>& getSector(int *) const;

    // Get the dims
    const int* getDims() const;

    // Get the widths
    const RealType* getWidths() const;

    // Get the total number of sectors
    int getNumSectors() const;

    // Get the skin depth
    RealType getSkinDepth() const;

    // Get the cutoff
    RealType getCutoff() const;

    // Get the number of remakes
    int getNumberOfRemakes() const;

    // --- Mutators

    // Set the skin depth
    void setSkinDepth(RealType);

    void setCutoffFactor(RealType);
    
    // GFlow is a friend class
    friend class GFlow;
    friend class SectorizationData;
    
  private:
    // --- Helper functions

    // Create sectors
    inline void makeSectors();

    inline RealType maxMotion();

    // Check whether we need to resectorize, if so, then do it
    inline void checkSectors();

    // Create verlet lists for all forces
    inline void makeVerletLists();

    // Given two particles, take care of whether they should exert forces on one another
    inline void pairInteraction(int, int);

    // Delete and set as null the xVL array
    inline void nullXVL();

    // Set up xVL array
    inline void setupXVL(int);

    // Bounds - get these from gflow
    Bounds bounds;

    // Number of sectors in each dimension
    int dims[DIMENSIONS];
    
    // The id's of particles in each sector - these are the actual sectors
    Array< vector<int> > sectors;

    // The widths of sectors in each dimension
    RealType widths[DIMENSIONS];
    // The inverse widths of sectors in each dimension
    RealType inverseW[DIMENSIONS];

    // The positions of the particles at the last verlet list creation
    RealType **xVL;
    int sizeXVL; // The size of the xVL array

    // Sectorization constants
    RealType skinDepth, cutoff, minCutoff, cutoffFactor;

    // Current head - for verlet list creation
    int currentHead;

    // The number of times we have remade the sectors
    int number_of_remakes;

    // Update timers and related
    RealType lastCheck, lastUpdate, updateDelay, max_update_delay;

    // The target move ratio for remake
    RealType mvRatioTollerance;
  };

};

#endif // __SECTORIZATION_HPP__GFLOW__