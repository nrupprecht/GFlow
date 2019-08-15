#ifndef __DOMAIN_BASE_HPP__GFLOW__
#define __DOMAIN_BASE_HPP__GFLOW__

#include "interactionhandler.hpp"

namespace GFlowSimulation {

  /** 
  *  \brief The base class for domain decomposition classes.
  *
  *  MPI tutorials: <http://mpitutorial.com/>
  *
  */
  class DomainBase : public InteractionHandler {
  public:
    //! Constructor
    DomainBase(GFlow*);

    //! Destructor
    virtual ~DomainBase();

    // --- Accessors

    //! Get the array of the number of cells in each dimension
    const int* getDims() const;

    //! Get the array of the width of the cells in each dimension
    const RealType* getWidths() const;

    //! Get the total number of cells in the domain
    int getNumCells() const;

    //! \brief Get the min small cutoff.
    RealType getCutoff() const;

    // --- Mutators

    //! \brief Set the skin depth. This function is virtual, as the inheriting class
    //! may need to remake itself after doing this.
    virtual void setSkinDepth(RealType) override;

    //! \brief Set the cell size. 
    //!
    //! Really, this suggests a cell size. It must be larger than the minimum possible cell size, 
    //! and must evenly divide the size of the domain.
    virtual void setCellSize(RealType)=0;

    // GFlow is a friend class
    friend class GFlow;

  protected:

    //! \brief Calculates the domain cell dimensions, widths, and inverse widths given 
    //! that the cutoff has been calculated.
    virtual void calculate_domain_cell_dimensions()=0;

    // --- Data

    //! \brief Number of cells in each dimension
    int *dims;

    //! \brief The widths of a cell in each dimension
    RealType *widths;

    //! \brief The inverse widths of a cell in each dimension
    RealType *inverseW;

    // --- Sectorization constants

    //! \brief The target cell size. 
    //!
    //! Cells will be at least this wide in each dimension, but since an integral number of them have
    //! to fit in the domain in each dimension, the actual widths will be different. They will be at
    //! least this large though.
    RealType target_cell_size = 0.;

    //! \brief The target number of particles in each particle's verlet list.
    //!
    //! Despite the fact that particles come in integer numbers, the target list size is a real, since it
    //! determines the skin depth via multiplication, and the *average* list size can be a real.
    RealType target_list_size = 4.;

    //! \brief The minimum allowable cutoff for small particles, 2*max_small_sigma + skin_depth
    RealType min_small_cutoff = 0.; 
  };

}
#endif // __DOMAIN_BASE_HPP__GFLOW__