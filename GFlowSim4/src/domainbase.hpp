#ifndef __DOMAIN_BASE_HPP__GFLOW__
#define __DOMAIN_BASE_HPP__GFLOW__

#include "gflow.hpp"

namespace GFlowSimulation {

  /** @brief The base class for domain decomposition and sectorization classes.
  *
  *  DomainBase classes are responsible for a subvolume (possibly the whole volume)
  *  of the simulation, keeping track of where the particles are and where they need 
  *  to be transfered to (if using MPI parallelization).\n
  *
  *  They are responsible for creating verlet lists for all the forces (if deemed 
  *  necessary) during the pre-forces step of the simulation.\n
  *  
  *  The natural way to do this is to divide the domain into sectors, or cells, and 
  *  slot particles into these cells so we can easily lookup potential neighbors of 
  *  each particle.\n
  *
  */
  class DomainBase : public Base {
  public:
    //! Constructor
    DomainBase(GFlow*);

    //! Destructor
    ~DomainBase();

    virtual void pre_forces();

    // --- Accessors

    //! Get the array of the number of cells in each dimension
    const int* getDims() const;

    //! Get the array of the width of the cells in each dimension
    const RealType* getWidths() const;

    //! Get the total number of cells in the domain
    int getNumCells() const;

    //! Get the skin depth the domain is using
    RealType getSkinDepth() const;

    //! Get the cutoff
    RealType getCutoff() const;

    //! Get the number of times the domain has made verlet lists
    int getNumberOfRemakes() const;

    // --- Mutators

    //! Set the skin depth. This function is virtual, as the inheriting class
    //! may need to remake itself after doing this.
    virtual void setSkinDepth(RealType);

    //! Set the cutoff factor. This function is virtual, as the inheriting class
    //! may need to remake itself after doing this.
    virtual void setCutoffFactor(RealType);

    // GFlow is a friend class
    friend class GFlow;

  protected:
    // --- Helper functions

    //! Remake the verlet lists
    virtual void remake_verlet()=0;
    
    //! Delete and set as null the xVL array
    void nullXVL();

    //! Set up xVL array
    void setupXVL(int);

    //! Fill the xVL array with the positions
    void fillXVL();

    //! Find the maximum amount two particles might have moved closer to one another
    virtual RealType maxMotion();

    //! Check whether particles might have move far enough to warrant verlet list remake.
    virtual bool check_needs_remake();

    // --- Data

    //! The bounds of the domain
    Bounds domain_bounds;
    //! The bounds of the entire simulation
    Bounds bounds;

    //! The number of times we have remade the sectors
    int number_of_remakes;

    //! Number of cells in each dimension
    int dims[DIMENSIONS];

    //! The widths of a cell in each dimension
    RealType widths[DIMENSIONS];

    //! The inverse widths of a cell in each dimension
    RealType inverseW[DIMENSIONS];

    // Sectorization constants
    RealType skinDepth, cutoff, minCutoff, cutoffFactor;

    // Update timers and related
    RealType lastCheck, lastUpdate, updateDelay, max_update_delay;

    //! The target move ratio for remake
    RealType mvRatioTollerance;

    //! An array storing the positions of the particles at the last verlet list creation
    RealType **xVL;
    //! The size of the xVL array
    int sizeXVL;

  };

}
#endif // __DOMAIN_BASE_HPP__GFLOW__