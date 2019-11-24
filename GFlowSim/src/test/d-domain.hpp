#ifndef __D_DOMAIN_HPP__GFLOW__
#define __D_DOMAIN_HPP__GFLOW__

#include "container-layout.hpp"
#include "../utility/bounds.hpp"

namespace GFlowSimulation {

  template<int dims> using ivec = vec<dims, int>;

  template<int sim_dimensions, template<int> class Container=ParticleContainer> 
  class DomainD {
  public:

    DomainD(Container<sim_dimensions> *prts) : particles(prts), process_bounds(sim_dimensions), simulation_bounds(sim_dimensions) {};

    void initialize() {
      // Set target cell size. Cells must be at least this large.
      target_cell_size = 2*max_small_radius + skin_depth;

      // Set up bounds.
      int prod = 1;
      for (int d=0; d<sim_dimensions; ++d) {
        dimensions[d] = static_cast<int>(process_bounds.wd(0) / target_cell_size);
        widths[d] = process_bounds.wd(0) / dimensions[d];
        inv_widths[d] = 1./widths[d];
        products[d] = prod;
        prod *= dimensions[d];
      }



    }

    void structure_updates() {
      auto x = particles->X();
      int size = particles->size();

      if (first_occupant.size()<size) {
        first_occupant = vector<int>(size, -1);
        next_occupant = vector<int>(size, -1);
      }

      for (int i=0; i<size; ++i) {
        vec<sim_dimensions> xc = x[i];
        hadamard_equals(xc, inv_widths);
        
        // Convert to int.
        ivec<sim_dimensions> index;
        index.copy_vec_cast(xc);
        int linear = index*products;
      }
    }

  private:

    real skin_depth = 0.01;
    real max_small_radius = 0.05;
    real target_cell_size;

    //! \brief The k-th entry is the local id of the first occupant of the k-th cell, or -1 if the cell is empty.
    vector<int> first_occupant;
    //! \brief the k-th entry is the next particle in the same cell as the k-th particle, or -1 if there is no such particle.
    vector<int> next_occupant;

    //! \brief The container of particles.
    Container<sim_dimensions> *particles = nullptr;

    vec<sim_dimensions> widths;
    vec<sim_dimensions> inv_widths;

    //! \brief The dimensionality of the sector grid.
    ivec<sim_dimensions> dimensions;

    ivec<sim_dimensions> dim_shift_up[sim_dimensions];
    ivec<sim_dimensions> dims_shift_down[sim_dimensions];

    ivec<sim_dimensions> products;

    //! \brief The bounds of the domain
    Bounds process_bounds;
    
    //! \brief The bounds of the entire simulation
    Bounds simulation_bounds;

  };
}
#endif // __D_DOMAIN_HPP__GFLOW__