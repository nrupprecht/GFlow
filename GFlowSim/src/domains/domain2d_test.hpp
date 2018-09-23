#ifndef __DOMAIN_2D_TEST_HPP__GFLOW__
#define __DOMAIN_2D_TEST_HPP__GFLOW__

#include "../base/domainbase.hpp"

namespace GFlowSimulation {

  /**
  *  @brief Test Cell class
  *
  */
  struct Cell2d {
    Cell2d() : x(nullptr), dx(nullptr), f(nullptr), sg(nullptr), array_length(0), loaded(false) {};

    // A list of all the particles in this cell
    vector<int> id_list;

    // Clear the id list
    void clear() { id_list.clear(); }

    // Add an element to the id list
    void add(int id) { id_list.push_back(id); }

    // The number of particles in the cell
    int size() { return id_list.size(); }

    // Helper arrays for calculating forces
    RealType **x, **f, *sg;
    RealType **dx; // Displacements between another particle these particles
    int array_length;
    bool loaded;
  };


  /**
  *  @brief Test 2d Domain class
  *
  */
  class Domain2d : public DomainBase {
  public:
    //! @brief Constructor
    Domain2d(GFlow*);

    //! Create cells, assign their neighbors, etc.
    virtual void initialize();

    // Pre-integrate calls sectorize
    virtual void pre_integrate();

    //! Exchange particles between processors
    virtual void exchange_particles() {};

    // --- Locator functions

    //! @brief Get all the particles within a radius of another particle
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.\n
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, RealType, vector<int>&) {};

    // --- Mutators

    //! @brief Set the cell size. 
    //!
    //! Really, this suggests a cell size. It must be larger than the minimum possible cell size, 
    //! and must evenly divide the size of the domain. Inherited from DomainBase.
    virtual void setCellSize(RealType) {};

    void calculateForces();

  private:

    //! @brief Remake the verlet lists for all the forces.
    //!
    //! Resectorizes the particles into cells and calculates verlet lists from the cell decomposition.
    virtual void construct();

    virtual bool check_needs_remake();

    inline void create_cells();

    inline void fill_cells();

    inline void load_cell(int, int);

    inline void load_cell(int);

    inline void load_cell(Cell2d&);

    inline void release_cell(int, int);

    inline void release_cell(int);

    inline int getCellIndex(RealType*);

    inline Cell2d& getCell(int, int);

    // A list of all the cells
    vector<Cell2d> cells;

  };

}
#endif // __DOMAIN_2D_TEST_HPP__GFLOW__