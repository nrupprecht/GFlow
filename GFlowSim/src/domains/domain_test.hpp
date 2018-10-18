
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
    vector<int> particle_ids;
    vector<CellTest*> adjacent;
  };



  class DomainTest : public DomainBase {
  public:
    DomainTest(GFlow*);

    ~DomainTest();

    virtual void initialize() override;

    virtual void pre_integrate() override;

    //! Exchange particles between processors
    virtual void exchange_particles() override;

    //! @brief Get all the particles within a radius of another particle
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.\n
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, RealType, vector<int>&) override;

    virtual void construct() override;

    virtual void setCellSize(RealType) override;

  private:

    inline void create_cells();

    inline void clear_cells();

    inline void fill_cells();

    //! @brief Turns a linear cell index into a (DIMENSIONS)-dimensional index
    inline void linear_to_tuple(const int, int*);

    //! @brief Turns a (DIMENSIONS)-dimensional index into a linear cell index.
    inline void tuple_to_linear(int&, const int*);

    inline bool correct_index(int*);

    inline int get_cell_index(const RealType*);

    int *border_type_up = nullptr;
    int *border_type_down = nullptr;

    //! @brief Whether the domain has been initialized or not.
    bool initialized = false;

    //! @brief A vector holding all the cells in the domain.
    vector<CellTest> cells;
  };

}
#endif // __DOMAIN_TEST_HPP__GFLOW__
