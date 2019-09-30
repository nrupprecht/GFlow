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

    //! \brief Initialize a domain type interaction handler, in part by calling polymorphic functions that child classes overload.
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

  protected:

    //! \brief Calculates the domain cell dimensions, widths, and inverse widths given 
    //! that the cutoff has been calculated.
    virtual void calculate_domain_cell_dimensions()=0;

    //! \brief Create the cells.
    virtual void create_cells()=0;

    // --- Helper functions

    //! \brief Get the tuple index of a cell that a position lies in.
    void get_cell_index_tuple(const RealType*, vector<int>&);

    //! \brief Get the linear index of the cell that a position lies within.
    int get_cell_index(const RealType*);

    //! \brief Turns a linear cell index into a (DIMENSIONS)-dimensional index
    void linear_to_tuple(const int, vector<int>&);
    void linear_to_tuple(const int, int*);

    //! \brief Turns a (DIMENSIONS)-dimensional index into a linear cell index.
    void tuple_to_linear(int&, const vector<int>&);
    void tuple_to_linear(int&, const int*);

    // --- Data

    //! \brief Number of cells in each dimension
    vector<int> dims;

    //! \brief The widths of a cell in each dimension
    vector<RealType> widths;

    //! \brief The inverse widths of a cell in each dimension
    vector<RealType> inverseW;

    //! \brief The number of halo or ghost sectors added below.
    //!
    //! We don't actually need this number (as of now), we only need dim_shift_down, but we keep it for completeness.
    vector<int> dim_shift_up;

    //! \brief The number of halo or ghost sectors added below.
    vector<int> dim_shift_down;

    //! \brief Helper array for doing indexing.
    vector<int> products;

    //! \brief The target cell size. 
    //!
    //! Cells will be at least this wide in each dimension, but since an integral number of them have
    //! to fit in the domain in each dimension, the actual widths will be different. They will be at
    //! least this large though.
    RealType target_cell_size = 0.;
  };

}
#endif // __DOMAIN_BASE_HPP__GFLOW__