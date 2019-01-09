
#ifndef __DOMAIN_TEST_HPP__GFLOW__
#define __DOMAIN_TEST_HPP__GFLOW__

#include "../base/domainbase.hpp"
#include "cell.hpp"

namespace GFlowSimulation {

  class Domain : public DomainBase {
  public:
    Domain(GFlow*);

    virtual ~Domain() override;

    virtual void initialize() override;

    virtual void pre_integrate() override;

    //! \brief Exchange particles between processors
    virtual void exchange_particles() override;

    //! \brief Get all the particles within a radius of another particle
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

    //! \brief Turns a linear cell index into a (DIMENSIONS)-dimensional index
    inline void linear_to_tuple(const int, int*);

    //! \brief Turns a (DIMENSIONS)-dimensional index into a linear cell index.
    inline void tuple_to_linear(int&, const int*);

    inline bool correct_index(int*);

    inline void get_cell_index_tuple(const RealType*, int*);

    inline int get_cell_index(const RealType*);

    //! \brief Array of products, used to compute linear indices from vectors or tuple indices.
    int *products = nullptr;

    //! \brief What type of borders there are in the "up" directions.
    //!
    //! 0 - None, 1 - Halo, 2 - Wrap.
    int *border_type_up   = nullptr;
    int *border_type_down = nullptr;

    //! \brief Whether the domain has been initialized or not.
    bool initialized = false;

    //! \brief The number of particles last time we built the cells.
    int number = 0;

    //! \brief A vector holding all the cells in the domain.
    vector<Cell> cells;
  };

}
#endif // __DOMAIN_TEST_HPP__GFLOW__
