
#ifndef __DIMS_DOMAIN_TEST_HPP__GFLOW__
#define __DIMS_DOMAIN_TEST_HPP__GFLOW__

#include "../base/domainbase.hpp"
#include "cell.hpp"

#include "../test/d-vec.hpp"

namespace GFlowSimulation {

  template<int dimensions> 
  class DimsDomain : public InteractionHandler {
  public:
    //! \brief Default constructor.
    DimsDomain(GFlow*);

    //! \brief Get all the particles within a radius of another particle
    //!
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, vector<int>&, RealType=-1.) override;

    //! \brief Get all the particles withing a radius of some position.
    virtual void getAllWithin(Vec, vector<int>&, RealType=-1.) override;

    //! \brief Construct interactions for a single particle.
    //!
    //! If the insert flag is set to true, the particle will first be added to the relevent domain or data structure
    //! (if applicable), if not, it will not be added.
    //! This function should be called after construct.
    virtual void constructFor(int, bool=false) override;

    //! \brief Traverse the cells structure, calling a function on pairs of particles that are close to one another.
    //!
    //! The function that is passed in should expect to receive particles' id1, id2, wrapping type (0 - no wrapping
    //! required, 1 - wrapping required), radius of particle 1, radius of particle 2, and distance between particles.
    virtual void traversePairs(PairFunction) override;

    //! \brief Function that traverses ghost particles, calling a function on all ghost-real pairs of particls that are within
    //! cutoff + skin depth distance of one another.
    //!
    //! The function (a PairFunction) that is passed in should expect to receive particles' id1, id2, wrapping type (0 - no wrapping
    //! required, 1 - wrapping required), radius of particle 1, radius of particle 2, and distance between particles.
    virtual void traverseGhostPairs(PairFunction) override;

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

    //! \brief Get the linear index of the cell that a position lies within.
    int get_halo_cell_index(const RealType*);

    //! \brief Turns a linear cell index into a (DIMENSIONS)-dimensional index
    void linear_to_tuple(const int, vector<int>&);
    void linear_to_tuple(const int, int*);

    //! \brief Turns a (DIMENSIONS)-dimensional index into a linear cell index.
    void tuple_to_linear(int&, const vector<int>&);
    void tuple_to_linear(int&, const int*);

    // --- Data

    //! \brief Number of cells in each dimension
    int grid[dimensions];

    //! \brief The widths of a cell in each dimension
    real widths[dimensions];

    //! \brief The inverse widths of a cell in each dimension
    real inverseW[dimensions];

    //! \brief The number of halo or ghost sectors added below.
    //!
    //! We don't actually need this number (as of now), we only need dim_shift_down, but we keep it for completeness.
    int dim_shift_up[dimensions];

    //! \brief The number of halo or ghost sectors added below.
    int dim_shift_down[dimensions];

    //! \brief Helper array for doing indexing.
    int products[dimensions];


    //! \brief The first particle in a cell. If the cell is empty, heads[cell_number] is -1.
    //!
    //! The size of this array is the number of cells in the domain.
    int *heads = nullptr;

    //! \brief The number of cells in the structure.
    int num_cells = 0;

    //! \brief Keep track of the next particle in a linked list. The next particle after id in the same list is link[id], which
    //! is -1 if there is no next particle.
    //!
    //! The size of this array is the number of particles on the processor.
    int *link = nullptr;

    //! \brief The number of particles that can be stored in the link array.
    int capacity = 0;

    //! \brief The target cell size. 
    //!
    //! Cells will be at least this wide in each dimension, but since an integral number of them have
    //! to fit in the domain in each dimension, the actual widths will be different. They will be at
    //! least this large though.
    real target_cell_size = 0.;
  };

  #include "dims-domain.tpp"

}
#endif // __DOMAIN_TEST_HPP__GFLOW__
