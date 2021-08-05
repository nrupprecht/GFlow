// Other files
#include "../gflow.hpp"
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../base/forcemaster.hpp"
#include "../base/interaction.hpp"
#include "../base/topology.hpp"
#include "../utility/generic-dimension.hpp"

template<int dimensions> 
DimsDomain::DimsDomain(GFlow *gflow) : DomainBase(gflow) {};

template<int dimensions> 
void DimsDomain::getAllWithin(int id1, vector<int>& neighbors, RealType distance) {
  
}

template<int dimensions> 
void DimsDomain::getAllWithin(Vec X, vector<int>& neighbors, RealType distance) {
  
}

template<int dimensions> 
void DimsDomain::constructFor(int id1, bool insert) {
  
}

template<int dimensions> 
void DimsDomain::traversePairs(std::function<void(int, int, int, RealType, RealType, RealType)> body) {
  // Get the array of max cutoffs
  const vector<RealType> & max_cutoffs = forceMaster->getMaxCutoff();

  // Checks
  if (cutoff_grid==nullptr) return;

  // Maximum reasonable distance.
  RealType max_reasonable = sqr(0.5*simulation_bounds.wd(0));
  for (int d=1; d<dimensions; ++d) {
    RealType mr = sqr(0.5*simulation_bounds.wd(d));
    if (mr<max_reasonable) max_reasonable = mr;
  }

  // Find potential neighbors
  auto x = simData->X();
  auto rd = simData->Sg();
  auto type = simData->Type();

  // Go through all the cells in the simulation.
  for (const auto &c : cells) {
    for (auto p=c.particle_ids.begin(); p!=c.particle_ids.end(); ++p) {
      // The id of the particle
      int id1 = *p;
      if (type(id1)<0) continue;
      RealType *cutoffs_id1 = cutoff_grid[type(id1)];
      RealType sigma1 = rd(id1)*max_cutoffs[type(id1)];

      // If sigma is <= than max_small_sigma, only look through cell stencil
      if (sigma1<=max_small_sigma) {
        // All other particles in the same sector
        for (auto q = p+1; q!=c.particle_ids.end(); ++q) {
          int id2 = *q;
          // If the other particle is a large particle, it will take care of this interaction
          if (type(id1)<0 || rd(id2)*max_cutoffs[type(id2)]>max_small_sigma) continue;
          // Look for distance between particles
          RealType r2 = getDistanceSqrNoWrap(x(id1), x(id2), dimensions);
          if (r2 < sqr((rd(id1) + rd(id2)*cutoffs_id1[type(id2)] + skin_depth)))
            body(id1, id2, 0, rd(id1), rd(id2), r2);
        }
        // Seach through list of adjacent cells
        for (const auto &d : c.adjacent)
          for (const auto id2 : d->particle_ids) {
            // If the other particle is a large particle, it will take care of this interaction
            if (type(id2)<0 || rd(id2)*max_cutoffs[type(id2)]>max_small_sigma) continue;
            // Look for distance between particles
            RealType r2 = getDistanceSqrNoWrap(x(id1), x(id2), dimensions);
            if (r2 < sqr((rd(id1) + rd(id2))*cutoffs_id1[type(id2)] + skin_depth))
              body(id1, id2, 0, rd(id1), rd(id2), r2);
            else if (r2>max_reasonable)
              body(id1, id2, 1, rd(id1), rd(id2), r2);                
          }
      }
    }
  }
}

template<int dimensions> 
void DimsDomain::traverseGhostPairs(PairFunction body) {
  // // Checks
  // if (cutoff_grid==nullptr) return;

  // // Maximum reasonable distance.
  // RealType max_reasonable = sqr(0.5*simulation_bounds.wd(0));
  // for (int d=1; d<dimensions; ++d) {
  //   RealType mr = sqr(0.5*simulation_bounds.wd(d));
  //   if (mr<max_reasonable) max_reasonable = mr;
  // }

  // // Get data.
  // auto x = simData->X(), x_g = simData->X<1>();
  // auto rd = simData->Sg(), rd_g = simData->Sg<1>();
  // auto type = simData->Type(), type_g = simData->Type<1>();

  // // Traverse all ghost particles.
  // for (int id1 = 0; id1<simData->size_ghosts(); ++id1) {
  //   // The index of the cell that ghost particle id1 is in.
  //   int index = get_halo_cell_index(x_g(id1));
  //   RealType *cutoffs_id1 = cutoff_grid[type_g(id1)];
  //   // Look at particles in the same cell and in neighboring cells.
  //   auto cell = cells[index];
  //   // Neighboring cells. The same cell is one of its "neighbors."
  //   for (auto c : cell.adjacent) {
  //     for (auto id2 : c->particle_ids) {
  //       RealType r2 = getDistanceSqrNoWrap(x_g(id1), x(id2), dimensions);
  //       if (r2 < sqr((rd_g(id1) + rd(id2))*cutoffs_id1[type(id2)] + skin_depth) || max_reasonable<r2)
  //         body(id1, id2, 2, rd_g(id1), rd(id2), r2);
  //     }
  //   }
  // }
}

template<int dimensions> 
inline void DimsDomain::structure_updates() {
  // We should have just done a particle removal, so we can use number, not size (since all arrays are compressed)
  int number = simData->number_owned(); // Ghost and halo particles are not binned.
  auto x = simData->X();

  // Make sure there is enough space in the links array.
  if (capacity<number) {
    if (links) delete [] links;
    links = new int[number];
  }

  // Clear cells, linked lists.
  for (int i=0; i<num_cells; ++i) heads[i] = -1;
  for (int i=0; i<number; ++i) links[i] = -1;

  real corner[dimensions];
  copy_vec<dimensions>(topology->getSimulationBounds().min, corner);
  
  int index[dimensions];
  // Bin all the particles
  for (int i=0; i<number; ++i) {
    // Get the index of the particle.
    subtract_vec<dimensions>(x(i), corner);
    // Add the particle to the linked list structure.

  }

}

template<int dimensions> 
void DimsDomain::calculate_domain_cell_dimensions() {
  for (int d=0; d<dimensions; ++d) {
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

template<int dimensions> 
inline void DimsDomain::create_cells() {
  // --- Create the cells
  // Get the total number of cells - The dims MUST be set first.
  num_cells = 1;
  for (int d=0; d<sim_dimensions; ++d) num_cells *= dims[d];

  if (heads) delete [] heads;
  heads = new int[num_cells];
}

template<int dimensions> 
inline bool DimsDomain::correct_index(vector<int>& index, bool wrap) {
  bool good_index = true;
  const BCFlag *bcs = gflow->getBCs();
  for (int d=0; d<dimensions; ++d) {
    // Going HIGHER than the border.
    if (index[d]>=dims[d]) {
      // Check if we should wrap
      if (bcs[d]==BCFlag::WRAP && border_type_up[d]==0 && wrap)
        index[d] -= dims[d];
      // If not, this is a bad index
      else good_index = false;
    }
    // Going LOWER than the border.
    if (index[d]<0) {
      // Check if we should wrap
      if (bcs[d]==BCFlag::WRAP && border_type_down[d]==0 && wrap) 
        index[d] += dims[d];
      // If not, this is a bad index
      else good_index = false;
    }
  }
  return good_index;
}

template<int dimensions> 
inline void DimsDomain::add_to_cell(const RealType *x, int id) {
  int linear = get_cell_index(x);
  // Stores the *local* id of the particle
  cells[linear].particle_ids.push_back(id);
}

template<int dimensions> 
inline bool DimsDomain::is_halo_cell(const vector<int>& index) {
  // Test of the cell is a halo cell.
  for (int d=0; d<dimensions; ++d) {
    if (index[d]<dim_shift_down[d] || dims[d]-dim_shift_up[d]<=index[d]) return true;
  }
  // If the cell failed all the tests, it is not a halo cell.
  return false;
}

