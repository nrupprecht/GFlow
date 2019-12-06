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
      int number_of_cells = 1;
      for (int d=0; d<dims; ++d) {
        dimensions[d] = static_cast<int>(process_bounds.wd(0) / target_cell_size);
        widths[d] = process_bounds.wd(0) / dimensions[d];
        inv_widths[d] = 1./widths[d];
        products[d] = number_of_cells;
        number_of_cells *= dimensions[d];
      }

      // Set up arrays.
      first_occupant = vector<int>(number_of_cells, -1);
      number_of_neighbors = vector<int>(number_of_cells, 0);
      neighbor_cells = vector<int>(number_of_cells, 0);
    }

    void construct() {

      structure_updates();

      verlet_list.clear();

      // The interaction function
      auto interaction_function = [&] (int id1, int id2, int w_type, RealType, RealType, RealType) {
        verlet_list.push_back(std::make_pair(id1, id2));
      };

      traversePairs(interaction_function);

    }

    vector<pair<int, int> > verlet_list;

    void structure_updates() {
      auto x = particles->X();
      int size = particles->size();

      // Make sure array is large enough.
      if (next_occupant.size()<size) next_occupant = vector<int>(size, -1);

      // Reset
      for (int i=0; i<size; ++i) next_occupant[i] = -1;
      for (int i=0; i<first_occupant.size(); ++i) first_occupant[i] = -1;

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

        // Insert into cell
        next_occupant[i] = first_occupant[linear];
        first_occupant[linear] = i;
      }
    }

    void traversePairs(std::function<void(int, int, int, RealType, RealType, RealType)> body) {
      // Get particle data
      auto x = particles->X();
      auto r = particles->R();
      auto type = particles->Type();

      int id0, id1;

      // Walk through cells.
      ivec<dims> cell1(0,0);
      for (cell1[1]=0; cell1[1]<dimensions[1]; ++cell1[1])
        for (cell1[0]=0; cell1[0]<dimensions[0]; ++cell1[0]) {
          // First particle in the central cell.
          int linear = cell1*products;
          id0 = first_occupant[linear];
          if (id0==-1) continue;

          // Data for finding neighboring cells.
          int first_neighbor = neighbor_cells[linear];
          int n_neighbors = number_of_neighbors[linear];

          while (id0!=-1) {
            // Search through all pairs of particles in the same cell.
            id1 = next_occupant[id0]; // Next particle in the cell.
            while (id1!=-1) {
              // Distance between particles
              real rsqr = distanceSqr(x(id0), x(id1));
              if (rsqr < sqr(r(id0) + r(id1) + skin_depth))
                body(id0, id1, 0, r(id0), r(id1), rsqr);
              // Increment particle.
              id1 = next_occupant[id1];
            }

            // --- Search through neighboring cells.
            for (int i=0; i<n_neighbors; ++i) {
              // Linear address of a neighboring cell.
              int neighbor_linear = neighbor_list[first_neighbor + i];
              // First particle in the neighboring cell.
              id1 = first_occupant[neighbor_linear];
              if (id1<0) continue;
              // Compare all particles in that cell to particle id0.
              while (id1!=-1) {
                real rsqr = distanceSqr(x(id0), x(id1));
                // If the particles are close, we don't need to use minimum image convention.
                if (rsqr < sqr(r(id0) + r(id1) + skin_depth))
                  body(id0, id1, 0, r(id0), r(id1), rsqr);
                else
                  body(id0, id1, 1, r(id0), r(id1), rsqr);
                // Increment particle.
                id1 = next_occupant[id1];
              }
            }

            // Next particle in the central cell.
            id0 = next_occupant[id0];
          }
        }
    }

  private:
    //! \brief The skin depth.
    real skin_depth = 0.01;
    real max_small_radius = 0.05;
    real target_cell_size;

    int number_of_cells = 0;

    //! \brief The k-th entry is the local id of the first occupant of the k-th cell, or -1 if the cell is empty.
    vector<int> first_occupant;
    //! \brief the k-th entry is the next particle in the same cell as the k-th particle, or -1 if there is no such particle.
    vector<int> next_occupant;

    //! \brief The number of neighbor cells that a cell has.
    vector<int> number_of_neighbors;
    //! \brief The location of the first neighbor cell that a cell has.
    vector<int> neighbor_cells;
    //! \brief A list that contains all the neighbors of cells. The arrays neighbor_cells and number_of_neighbors can be used to find
    //! which cells in this list are your neighbors.
    vector<int> neighbor_list;

    //! \brief The container of particles.
    Container *particles = nullptr;

    //! \brief The widths of a single cell.
    vec<dims> widths;
    //! \brief The inverse widths of a single cell.
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