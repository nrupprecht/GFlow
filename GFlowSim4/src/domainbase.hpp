#ifndef __DOMAIN_BASE_HPP__GFLOW__
#define __DOMAIN_BASE_HPP__GFLOW__

#include "gflow.hpp"

namespace GFlowSimulation {

  class DomainBase : public Base {
  public:
    //! Constructor
    DomainBase(GFlow*);

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

  protected:
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

    //! Update timers and related
    RealType lastCheck, lastUpdate, updateDelay, max_update_delay;

    //! The target move ratio for remake
    RealType mvRatioTollerance;

  };

}
#endif // __DOMAIN_BASE_HPP__GFLOW__