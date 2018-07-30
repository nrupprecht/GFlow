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
    void initialize();

    // Pre-integrate calls sectorize
    virtual void pre_integrate();

    // Pre-forces is where resectorization may happen
    virtual void pre_forces();

  private:
    // --- Helper functions

    //! @brief Adds particles to the correct position(s) in the sectorization cells.
    //!
    //! Add a particle to the cell(s) it belongs in. A particle and its images can only 
    //! be in one central cell, but can also be in multiple harmonic cells. If it needs
    //! to be in a harmonic cell, we create a ghost particle in simdata, and put the
    //! ghost in the the harmonic cell. We also add an alias mapping from the ghost to the
    //! actual particle in the simdata.
    void addToCell(int);

    //! Get the (linear) index of the cell the position falls in
    int getCellIndex(RealType *);

    //! Calculate the data neccessary to run domains in parallel e.g. domain index, which domain this is, etc.
    void parallel_assignments();

    //! Check whether we need to remake the verletl lists
    void check_needs_remake();

    // Turns a linear cell index into a (DIMENSIONS)-dimensional index
    void linear_to_tuple(const int, int*);

    // Turns a (DIMENSIONS)-dimensional index into a linear cell index
    void tuple_to_linear(int&, const int*);

    void find_adjacent_cells(int[DIMENSIONS]);

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

    //! The bounds of the domain
    Bounds domain_bounds;
    //! The bounds of the domain, including ghost and harmonic cells
    Bounds extended_domain_bounds;
    //! The bounds of the entire simulation
    Bounds bounds;

    //! If we have harmonic cells on the positive boundary.
    bool harmonic_up[DIMENSIONS];
    //! If we have harmonic cells on the negative boundary.
    bool harmonic_down[DIMENSIONS];

    bool ghost_up[DIMENSIONS];
    bool ghost_down[DIMENSIONS];

    // --- For remakes

    //! The positions of the particles at the last verlet list creation
    RealType **xVL;
    int sizeXVL; // The size of the xVL array

  };

}
#endif // __DOMAIN_HPP__GFLOW__