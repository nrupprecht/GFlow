
#ifndef __DOMAIN_TEST_HPP__GFLOW__
#define __DOMAIN_TEST_HPP__GFLOW__

#include "../base/domainbase.hpp"

namespace GFlowSimulation {

  inline RealType min(RealType *array, int length) {
    if (length==0) return 0;
    RealType m = array[0];
    for (int i=0; i<length; ++i)
      if (array[i]<m) m = array[i];
    return m;
  }


  struct CellTest {
    //! @brief Dimension setting constructor
    CellTest(int d);

    //! @brief Copy constructor
    CellTest(const CellTest& cell);

    //! @brief Operator equals.
    CellTest operator=(const CellTest& cell);

    int size() const;

    //! @brief The ids of the particles contained in the cell.
    vector<int> particle_ids;

    //! @brief Pointers to cells adjacent to the cell.
    vector<CellTest*> adjacent;

    //! @brief The number of dimensions of the simulation.
    int sim_dimensions;

    //! @brief True if particles in this cell should be create ghost cells.
    bool halo_adjacent;
  };



  class DomainTest : public DomainBase {
  public:
    DomainTest(GFlow*);

    virtual ~DomainTest() override;

    virtual void initialize() override;

    virtual void pre_integrate() override;

    //! Exchange particles between processors
    virtual void exchange_particles() override;

    //! @brief Get all the particles within a radius of another particle
    //!
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.\n
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, RealType, vector<int>&) override;

    virtual void removeOverlapping(RealType) override;

    virtual void construct() override;

    virtual void setCellSize(RealType) override;

  private:

    inline void create_cells();

    inline void update_cells();

    inline void clear_cells();

    inline void fill_cells();

    //! @brief Turns a linear cell index into a (DIMENSIONS)-dimensional index
    inline void linear_to_tuple(const int, int*);

    //! @brief Turns a (DIMENSIONS)-dimensional index into a linear cell index.
    inline void tuple_to_linear(int&, const int*);

    inline bool correct_index(int*);

    inline void get_cell_index_tuple(const RealType*, int*);

    inline int get_cell_index(const RealType*);

    //! @brief Array of products, used to compute linear indices from vectors or tuple indices.
    int *products = nullptr;

    //! @brief What type of borders there are in the "up" directions.
    //!
    //! 0 - None, 1 - Halo, 2 - Wrap.
    int *border_type_up   = nullptr;
    int *border_type_down = nullptr;

    //! @brief Whether the domain has been initialized or not.
    bool initialized = false;

    //! @brief The number of particles last time we built the cells.
    int number = 0;

    //! @brief A vector holding all the cells in the domain.
    vector<CellTest> cells;
  };

}
#endif // __DOMAIN_TEST_HPP__GFLOW__
