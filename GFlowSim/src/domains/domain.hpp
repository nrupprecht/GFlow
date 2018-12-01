#ifndef __DOMAIN_HPP__GFLOW__
#define __DOMAIN_HPP__GFLOW__

#include "../base/domainbase.hpp"
#include "cell.hpp"

namespace GFlowSimulation {

  /** @brief A (hyper) rectangular patch of space in the simulation.
  *  
  *  Creates and maintains a sectorization composed of many cells. Uses this to 
  *  create verlet lists for the forces. 
  *
  *  Responsible for communicating with neighboring domains.
  */
  class Domain : public DomainBase {
  public:
    //! @brief Default constructor - pass in the Gflow object
    Domain(GFlow *);

    //! @brief Destuctor.
    ~Domain();

    //! @brief Create cells, assign their neighbors, etc.
    virtual void initialize() override;

    //! @brief Pre-integrate calls sectorize
    virtual void pre_integrate() override;

    //! @brief Exchange particles between processors
    virtual void exchange_particles() override;

    // --- Locator functions

    //! @brief Get all the particles within a radius of another particle
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.\n
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, RealType, vector<int>&) override;

    virtual void removeOverlapping(RealType) override { throw false; }

    // --- Mutators

    virtual void setSkinDepth(RealType) override;

    //! @brief Set the cell size. 
    //!
    //! Really, this suggests a cell size. It must be larger than the minimum possible cell size, 
    //! and must evenly divide the size of the domain. Inherited from DomainBase.
    virtual void setCellSize(RealType) override;

  private:
    // --- Helper functions

    //! @brief Adds particles to the correct position(s) in the sectorization cells.
    //!
    //! We do not add halo particles when we inset the particle.
    void addToCell(int);

    //! @brief Get the (linear) index of the cell the position falls in
    int getCellIndex(RealType *);

    #if USE_MPI==1
    //! @brief Calculate the data neccessary to run domains in parallel e.g. domain index, which domain this is, etc.
    void parallel_assignments();
    #endif

    //! @brief Remake the verlet lists for all the forces.
    //!
    //! Resectorizes the particles into cells and calculates verlet lists from the cell decomposition.
    virtual void construct();

    //! @brief Turns a linear cell index into a (DIMENSIONS)-dimensional index
    inline void linear_to_tuple(const int, int*);

    //! @brief Turns a (DIMENSIONS)-dimensional index into a linear cell index.
    inline void tuple_to_linear(int&, const int*);

    //! @brief Find the linear indices of the cells adjacent to a given cell.
    inline void find_adjacent_cells(int*, bool, vector<Cell*>&);

    //! @brief Based on the cutoff distance, reate new cells, assign them their types, etc.
    inline void create_cells();

    //! @brief Fill the cells with particle ids.
    inline void fill_cells();

    inline bool correct_index(int*, bool);

    // --- Data

    //! All the cells in this decomposition
    vector<Cell> cells;

    // --- Parallel decomposition data

    //! @brief The adjacent domains in the positive dimension. For none, we have -1. Otherwise, we have the processor id.
    int *domains_up;
    //! @brief The adjacent domains in the negative dimension. For none, we have -1. Otherwise, we have the processor id.
    int *domains_down;

    //! @brief How many domains exist along each dimension (the decomposition into domains is rectangular).
    int *domain_dims;
    //! @brief The (DIMENSION)-tuple domain index of this domain
    int *domain_index;

    //! @brief The bounds of the domain, including ghost and harmonic cells
    Bounds extended_domain_bounds;

    //! @brief If we have harmonic cells on the positive boundary.
    bool *halo_up;
    //! @brief If we have harmonic cells on the negative boundary.
    bool *halo_down;
    //! @brief If we have ghost cells on the positive boundary
    bool *ghost_up;
    //! @brief If we have ghost cells on the negative boundary
    bool *ghost_down;

    //! @brief Whether to use halo cells to reflect the periodic boundary conditions
    bool use_halo_cells;

    // --- MPI related varibles
    #if USE_MPI==1
    //! @brief The rank of this processor.
    int rank;
    //! @brief The total number of MPI processors.
    int numProc;
    //! @brief This is true when we have done the parallel initialization.
    bool parallel_init;
    #endif 

    // --- Linked cell force related
    inline void load_cell(Cell&);
    inline void release_cell(Cell&);
  };

}
#endif // __DOMAIN_HPP__GFLOW__
