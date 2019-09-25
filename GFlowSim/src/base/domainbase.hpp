#ifndef __DOMAIN_BASE_HPP__GFLOW__
#define __DOMAIN_BASE_HPP__GFLOW__

#include "interactionhandler.hpp"

namespace GFlowSimulation {

  /** 
  *  \brief The base class for domain decomposition classes.
  *
  */
  class DomainBase : public InteractionHandler {
  public:
    //! Constructor
    DomainBase(GFlow*);

    virtual void initialize() override;

    // --- Accessors

    //! Get the array of the number of cells in each dimension
    const vector<int>& getDims() const;

    //! Get the array of the width of the cells in each dimension
    const vector<RealType>& getWidths() const;

    //! Get the total number of cells in the domain
    int getNumCells() const;

    //! \brief Get the min small cutoff.
    RealType getCutoff() const;

    // GFlow is a friend class
    //friend class GFlow;

  protected:

    //! \brief Calculates the domain cell dimensions, widths, and inverse widths given 
    //! that the cutoff has been calculated.
    virtual void calculate_domain_cell_dimensions()=0;

    //! \brief Create the cells.
    virtual void create_cells()=0;

    // --- Data

    //! \brief Number of cells in each dimension
    vector<int> dims;

    //! \brief The widths of a cell in each dimension
    vector<RealType> widths;

    //! \brief The inverse widths of a cell in each dimension
    vector<RealType> inverseW;

    // --- Sectorization constants

    //! \brief The target cell size. 
    //!
    //! Cells will be at least this wide in each dimension, but since an integral number of them have
    //! to fit in the domain in each dimension, the actual widths will be different. They will be at
    //! least this large though.
    RealType target_cell_size = 0.;

    //! \brief The minimum allowable cutoff for small particles, 2*max_small_sigma + skin_depth
    RealType min_small_cutoff = 0.; 
  };

}
#endif // __DOMAIN_BASE_HPP__GFLOW__