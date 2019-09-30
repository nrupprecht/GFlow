#include "domain-list.hpp"
// Other files
#include "../base/forcemaster.hpp"

namespace GFlowSimulation {

  DomainList::DomainList(GFlow*) : DomainBase(gflow) {};

  void DomainList::getAllWithin(int id, vector<int>& neighbors, RealType distance) {

  }

  void DomainList::getAllWithin(Vec, vector<int>& neighbors, RealType distance) {

  }

  void DomainList::constructFor(int id, bool insert) {

  }

  void DomainList::traversePairs(PairFunction body) {
    const int dimensions = 2;


    // Get the array of max cutoffs
    const vector<RealType> & max_cutoffs = forceMaster->getMaxCutoff();

    // Create some indexing arrays.
    int zero[dimensions];
    set1_vec<dimensions>(static_cast<int*>(zero), 0);
    int lower[dimensions], upper[dimensions];
    set1_vec<dimensions>(lower, -1);
    set1_vec<dimensions>(upper, 1);

    // Get some data.
    auto x = simData->X();
    auto sg = simData->Sg();
    auto type = simData->Type();

    // A function to call on each cell.
    auto neighbor_function = [&] (int *index) -> void {
      // The variable "index" is the index of the cell we are focusing on. Iterate through all particles in the cell.
      // Convert to linear index.
      int linear;
      tuple_to_linear(linear, index);
      // First particle in the cell.
      int id1 = cell_list[linear];
      // Cutoff data
      RealType *cutoffs_id1 = cutoff_grid[type[id1]];
      RealType sigma1 = sg[id1]*max_cutoffs[type[id1]];
      // If the cell is empty, return;
      if (id1<0) return;
      // Otherwise, do interactions.
      // Look at particles in the same cell.
      int id2 = particle_list[id1];
      while (id2!=-1) {
        // If the other particle is a large particle, it will take care of this interaction
        if (sg[id2]*max_cutoffs[type[id2]]<=max_small_sigma) {
          // Look for distance between particles
          RealType r2 = getDistanceSqrNoWrap(x[id1], x[id2], sim_dimensions);
          if (r2 < sqr((sg[id1] + sg[id2])*cutoffs_id1[type[id2]] + skin_depth))
            body(id1, id2, 0, sg[id1], sg[id2], r2);
        }
        // Iterate through linked list.
        id2 = particle_list[id2];
      }
      // Case: id1 is a Small particle
      if (sigma1<max_small_sigma) {

        // Function to run in the inner loop.
        auto inner_function = [&] (int *offset) -> void {
          // The variable "offset" is the cell offset.
          int other[dimensions];
          // Get the actual cell index.
          add_vec<dimensions>(offset, index, other);
          // Convert to linear index.
          int linear;
          tuple_to_linear(linear, other);
          // Check if this is a legitimate cell, or if the index needs wrapping.

          // FILL IN
        };

        // Iterate over neighboring cells.
        for_loop<dimensions>(inner_function, lower, upper);
      }
      // Case: id1 is a Large particle.
      else {

      }
    };

    // Loop through all interactions.
    for_loop<dimensions>(neighbor_function, zero, dims.data());
  }

  void DomainList::structure_updates() {
    // Clear cells
    for (auto& p : cell_list) p = -1;
    for (auto& p : particle_list) p = -1;

    // Fill cells
    for (int i=0; i<simData->size(); ++i) {
      // Get the cell the particle belongs in.
      int index = get_cell_index(simData->X(i));
      // Basically insert particle i into the front of the cell linked list.
      particle_list[i] = cell_list[index];
      cell_list[index] = i;
    }
  }

  void DomainList::calculate_domain_cell_dimensions() {
    for (int d=0; d<sim_dimensions; ++d) {
      dims[d] = static_cast<int>(process_bounds.wd(d)/target_cell_size);
      // Check that the bounds are good
      if (dims[d]<=0) throw BadBounds();
      // Min dims
      if (dims[d]<4) dims[d] = 4;
      // Set widths 
      widths[d] = process_bounds.wd(d)/dims[d];
      inverseW[d] = 1./widths[d];
      // Do border related work
      if (border_type_down[d]) {
        ++dims[d];
        dim_shift_down[d] = 1;
      }
      else dim_shift_down[d] = 0;
      if (border_type_up[d]) {
        ++dims[d];
        dim_shift_up[d] = 1;
      }
      else dim_shift_up[d] = 0;
    }
  }

  void DomainList::create_cells() {
    const int dimensions = 2;
    cell_list = vector<int>(product<dimensions>(dims), -1);
  }

}