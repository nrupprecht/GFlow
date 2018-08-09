#ifndef __DOMAIN_HPP__GFLOW__
#define __DOMAIN_HPP__GFLOW__

#include "domainbase.hpp"
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
    //! Default constructor - pass in the Gflow object
    Domain(GFlow *);

    //! Create cells, assign their neighbors, etc.
    virtual void initialize();

    // Pre-integrate calls sectorize
    virtual void pre_integrate();

    //! Exchange particles between processors
    virtual void exchange_particles();

    // --- Locator functions

    //! @brief Get all the particles within a radius of another particle
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.\n
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, RealType, vector<int>&);

    // --- Mutators

    virtual void setSkinDepth(RealType);

    //! @brief Set the cell size. 
    //!
    //! Really, this suggests a cell size. It must be larger than the minimum possible cell size, 
    //! and must evenly divide the size of the domain. Inherited from DomainBase.
    virtual void setCellSize(RealType);

  private:
    // --- Helper functions

    //! @brief Adds particles to the correct position(s) in the sectorization cells.
    //!
    //! We do not add halo particles when we inset the particle.
    void addToCell(int);

    //! @brief Get the (linear) index of the cell the position falls in
    int getCellIndex(RealType *);

    //! @brief Calculate the data neccessary to run domains in parallel e.g. domain index, which domain this is, etc.
    void parallel_assignments();

    //! @brief Remake the verlet lists for all the forces.
    //!
    //! Resectorizes the particles into cells and calculates verlet lists from the cell decomposition.
    virtual void remake_verlet();

    //! @brief Turns a linear cell index into a (DIMENSIONS)-dimensional index
    inline void linear_to_tuple(const int, int*);

    //! @brief Turns a (DIMENSIONS)-dimensional index into a linear cell index.
    inline void tuple_to_linear(int&, const int*);

    //! @brief Find the linear indices of the cells adjacent to a given cell.
    inline void find_adjacent_cells(int[DIMENSIONS], bool, vector<Cell*>&);

    //! @brief Based on the cutoff distance, reate new cells, assign them their types, etc.
    inline void create_cells();

    //! @brief Fill the cells with particle ids.
    inline void fill_cells();

    inline bool correct_index(int[DIMENSIONS], bool);

    // --- Data

    //! All the cells in this decomposition
    vector<Cell> cells;

    // --- Parallel decomposition data

    //! The adjacent domains in the positive dimension. For none, we have -1. Otherwise, we have the processor id.
    int domains_up[DIMENSIONS];
    //! The adjacent domains in the negative dimension. For none, we have -1. Otherwise, we have the processor id.
    int domains_down[DIMENSIONS];

    //! How many domains exist along each dimension (the decomposition into domains is rectangular).
    int domain_dims[DIMENSIONS];
    //! The (DIMENSION)-tuple domain index of this domain
    int domain_index[DIMENSIONS];

    //! The bounds of the domain, including ghost and harmonic cells
    Bounds extended_domain_bounds;

    //! If we have harmonic cells on the positive boundary.
    bool halo_up[DIMENSIONS];
    //! If we have harmonic cells on the negative boundary.
    bool halo_down[DIMENSIONS];
    //! If we have ghost cells on the positive boundary
    bool ghost_up[DIMENSIONS];
    //! If we have ghost cells on the negative boundary
    bool ghost_down[DIMENSIONS];

    //! Whether to use halo cells to reflect the periodic boundary conditions
    bool use_halo_cells;

  };

}
#endif // __DOMAIN_HPP__GFLOW__