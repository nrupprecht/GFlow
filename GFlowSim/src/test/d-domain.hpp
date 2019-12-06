#ifndef __D_DOMAIN_HPP__GFLOW__
#define __D_DOMAIN_HPP__GFLOW__

#include "container-layout.hpp"
#include "../utility/bounds.hpp"

namespace GFlowSimulation {

  

  template<int dims, DataLayout layout> 
  class DomainD : public Base {
  public:

    typedef ParticleContainer<dims, layout> Container;

    DomainD(GFlow *gflow) : Base(gflow), process_bounds(dims), simulation_bounds(dims) {};

    void setContainer(Container *ptr) {
      particles = ptr;
    }

    void initialize() {
      // Set target cell size. Cells must be at least this large.
      target_cell_size = 2*max_small_radius + skin_depth;
      // Make sure simulation bounds is correct.
      simulation_bounds = gflow->getBounds();
      // Set the process bounds.
      process_bounds = simulation_bounds;

      // Set up bounds.
      int prod = 1;
      for (int d=0; d<dims; ++d) {
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

      // The low corner. Particle positions should be relative to this.
      vec<dims> corner(process_bounds.min);

      for (int i=0; i<size; ++i) {
        vec<dims> xc;
        xc = x(i);
        xc -= corner;
        hadamard_equals(xc, inv_widths);

        // Convert to int.
        ivec<dims> index;
        index.copy_vec_cast(xc);
        int linear = index*products;

        //cout << x(i) << " -> " << xc << " " << index << " " << linear << endl;
        cout << x(i) << endl;
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
    Container *particles = nullptr;

    vec<dims> widths;
    vec<dims> inv_widths;


    //! \brief The dimensionality of the sector grid.
    ivec<dims> dimensions;

    ivec<dims> dim_shift_up[dims];
    ivec<dims> dims_shift_down[dims];

    ivec<dims> products;

    //! \brief The bounds of the domain
    Bounds process_bounds;
    
    //! \brief The bounds of the entire simulation
    Bounds simulation_bounds;

  };
}
#endif // __D_DOMAIN_HPP__GFLOW__