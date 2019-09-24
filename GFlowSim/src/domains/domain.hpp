
#ifndef __DOMAIN_TEST_HPP__GFLOW__
#define __DOMAIN_TEST_HPP__GFLOW__

#include "../base/domainbase.hpp"
#include "cell.hpp"

namespace GFlowSimulation {

  class Domain : public DomainBase {
  public:
    //! \brief Default constructor.
    Domain(GFlow*);

    //! \brief Initialize the domain.
    virtual void initialize() override;

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

    //! \brief Set the minimum cell size. 
    //!
    //! Causes a rebuild of all the cells.
    virtual void setCellSize(RealType) override;

    //! \brief Traverse the cells structure, calling a function on pairs of particles that are close to one another.
    //!
    //! The function that is passed in should expect to receive particles' id1, id2, wrapping type (0 - no wrapping
    //! required, 1 - wrapping required), radius of particle 1, radius of particle 2, and distance between particles.
    virtual void traversePairs(PairFunction) override;

  private:

    //! \brief Update the linked cells structure.
    virtual void structure_updates() override;

    //! \brief Calculates the domain cell dimensions, widths, and inverse widths given 
    //! that the cutoff has been calculated. 
    inline void calculate_domain_cell_dimensions() override;

    //! \brief Calculates the array of products used to convert between linear and tuple indices.
    inline void calculate_product_array();

    //! \brief Create the cells.
    inline void create_cells();

    //! \brief Updates the cells - this just clears the cells, fills them, and records how many particles there were.
    //!
    //! \see clear_cells
    //! \see fill_cells
    inline void update_cells();

    //! \brief Clear all the particle IDs from all the cells
    inline void clear_cells();

    //! \brief Put all particles in the domain into cells.
    inline void fill_cells();

    //! \brief Turns a linear cell index into a (DIMENSIONS)-dimensional index
    inline void linear_to_tuple(const int, vector<int>&);

    //! \brief Turns a (DIMENSIONS)-dimensional index into a linear cell index.
    inline void tuple_to_linear(int&, const vector<int>&);

    //! \brief Correct a linear index for wrapping. Returns true if the index is a valid cell.
    //!
    //! If the flag is set to false, we do not wrap positions.
    inline bool correct_index(vector<int>&, bool=true);

    //! \brief Get the tuple index of a cell that a position lies in.
    inline void get_cell_index_tuple(const RealType*, vector<int>&);

    //! \brief Get the linear index of the cell that a position lies within.
    inline int get_cell_index(const RealType*);

    //! \brief Add a particle to the cell it belongs in.
    inline void add_to_cell(const RealType*, int);

    //! \brief Helper array for doing indexing.
    vector<int> products;

    //! \brief The number of halo or ghost sectors added below.
    //!
    //! We don't actually need this number (as of now), we only need dim_shift_down, but we keep it for completeness.
    vector<int> dim_shift_up;

    //! \brief The number of halo or ghost sectors added below.
    vector<int> dim_shift_down;

    // //! \brief Array of products, used to compute linear indices from vectors or tuple indices.
    // int *products;

    // //! \brief What type of borders there are in the "up" directions.
    // //!
    // //! 0 - No ghost particles, 1 - Ghost particles, no wrapping, 2 - Ghost particles, wrapping.
    // int *border_type_up;

    // //! \brief What type of borders there are in the "down" directions.
    // //!
    // //! 0 - No ghost particles, 1 - Ghost particles, no wrapping, 2 - Ghost particles, wrapping.
    // int *border_type_down;

    // //! \brief The number of halo or ghost sectors added below.
    // //!
    // //! We don't actually need this number (as of now), we only need dim_shift_down, but we keep it for completeness.
    // int *dim_shift_up;
    // //! \brief The number of halo or ghost sectors added below.
    // int *dim_shift_down;

    //! \brief A vector holding all the cells in the domain.
    vector<Cell> cells;

    //! \brief Function used for setting an excluded particle when calling getAllWithin
    int _exclude = -1;
  };

}
#endif // __DOMAIN_TEST_HPP__GFLOW__
