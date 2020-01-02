#include "domain-list.hpp"
// Other files
#include "../base/forcemaster.hpp"

namespace GFlowSimulation {

  DomainList::DomainList(GFlow *gflow) : DomainBase(gflow) {};

  void DomainList::getAllWithin(int id, vector<int>& neighbors, RealType distance) {

  }

  void DomainList::getAllWithin(Vec, vector<int>& neighbors, RealType distance) {
    
  }

  void DomainList::constructFor(int id, bool insert) {
    
  }

  void DomainList::traversePairs(PairFunction body) {
    constexpr int dimensions = 2;

    if (sim_dimensions!=dimensions) throw BadDimension("Dimensions in domain-list is currently hard coded.");

    // Get the array of max cutoffs
    const vector<RealType> & max_cutoffs = forceMaster->getMaxCutoff();

    // Maximum reasonable distance.
    RealType max_reasonable = sqr(0.9*simulation_bounds.wd(0));
    for (int d=1; d<dimensions; ++d) {
      RealType mr = sqr(0.9*simulation_bounds.wd(d));
      if (mr<max_reasonable) max_reasonable = mr;
    }

    // Create some indexing arrays.
    int zero[dimensions];
    set1_vec<dimensions>(static_cast<int*>(zero), 0);
    int lower[dimensions], upper[dimensions];
    set1_vec<dimensions>(lower, -1);
    set1_vec<dimensions>(upper, 2); // For loop look for indices <2 i.e. <=1.

    // Get some data.
    auto x = simData->X();
    auto sg = simData->Sg();
    auto type = simData->Type();

    // Create stencil
    constexpr int stencil_size = power<dimensions>(3)/2; // A constexpr integer power function.
    int stencil[stencil_size][dimensions];
    // Lambda for creating a stencil
    auto create_stencil = [&] (int* index, int counts) -> void {
      copy_vec<dimensions>(index, stencil[counts]);
    }; 
    for_loop<dimensions>(create_stencil, lower, upper, stencil_size);

    // A function to call on each cell.
    auto neighbor_function = [&] (int *index) -> void {
      // The variable "index" is the index of the cell we are focusing on. Convert to linear index.
      int linear;
      tuple_to_linear(linear, index);
      
      // First particle in the cell.
      int id1 = cell_list[linear];

      // Iterate through all particles in the cell
      for (; id1!=-1; id1 = particle_list[id1]) {
        // Cutoff data
        RealType *cutoffs_id1 = cutoff_grid[type(id1)];
        RealType sigma1 = sg(id1)*max_cutoffs[type(id1)];

        // Iterate through all particles in the same cell.
        int id2 = particle_list[id1];
        for (; id2!=-1; id2 = particle_list[id2]) {
          // If the other particle is a large particle, it will take care of this interaction
          if (sg(id2)*max_cutoffs[type(id2)]<=max_small_sigma) {
            // Look for distance between particles
            RealType r2 = getDistanceSqrNoWrap(x(id1), x(id2), dimensions);
            if (r2 < sqr((sg(id1) + sg(id2))*cutoffs_id1[type(id2)] + skin_depth)) {
              body(id1, id2, 0, sg(id1), sg(id2), r2);
            }
          }
        }

        // Iterate through particles in neighboring cells.
        // Case: id1 is a Small particle
        if (sigma1<max_small_sigma) {

          for (int i=0; i<stencil_size; ++i) {
            // The variable "offset" is the cell offset.
            int other[dimensions];
            // Get the actual cell index.
            add_vec<dimensions>(stencil[i], index, other);

            // Check if this is a legitimate cell, or if the index needs wrapping.
            bool valid_cell = true;
            if (other[0]<0 || dims[0]<=other[0] || other[1]<0 || dims[1]<=other[1]) valid_cell = false;
            // ...

            if (valid_cell) {
              // Cell is now valid Convert to linear index.
              int other_linear;
              tuple_to_linear(other_linear, other);
            
              // Iterate through all particles in this other cell.
              int id2 = cell_list[other_linear];
              for (; id2!=-1; id2 = particle_list[id2]) {
                // If the other particle is a large particle, it will take care of this interaction
                if (sg(id2)*max_cutoffs[type(id2)]>max_small_sigma) continue;
                // Look for distance between particles
                RealType r2 = getDistanceSqrNoWrap(x(id1), x(id2), dimensions);
                if (r2 < sqr((sg(id1) + sg(id2))*cutoffs_id1[type(id2)] + skin_depth)) body(id1, id2, 0, sg(id1), sg(id2), r2);
                else if (r2>max_reasonable) body(id1, id2, 1, sg(id1), sg(id2), r2);
              }
            }

          }
        }
        // Case: id1 is a Large particle.
        else {
          
        }
      }
    };


    /*
    int index[2];
    for (index[1]=0; index[1]<dims[1]; ++index[1]) {
      for (index[0]=0; index[0]<dims[0]; ++index[0]) {
        neighbor_function(static_cast<int*>(index));
      }
    }
    */
    

    // Loop through all interactions.
    for_loop<dimensions>(neighbor_function, zero, dims.data());

  }

  void DomainList::structure_updates() {
    // Clear particle list
    if (particle_list.size()<simData->size_owned()) particle_list = vector<int>(simData->size_owned(), -1);
    else {
      for (auto& p : particle_list) p = -1;
    }

    // Clear cells
    for (auto& p : cell_list) p = -1;

    // Fill cells
    for (int i=0; i<simData->size_owned(); ++i) {
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