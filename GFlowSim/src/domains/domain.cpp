#include "domain.hpp"
// Other files
#include "../gflow.hpp"
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../base/forcemaster.hpp"
#include "../base/interaction.hpp"
#include "../base/topology.hpp"
#include "../utility/generic-dimension.hpp"

namespace GFlowSimulation {

  Domain::Domain(GFlow *gflow) : DomainBase(gflow) {};

  void Domain::getAllWithin(int id1, vector<int>& neighbors, RealType distance) {
    // Set exclude flag so the more general getAllWithin function doesn't count id1 as a valid particle.
    _exclude = id1;
    // Wrap the position vector with a vec object.
    Vec x(simData->X(id1), sim_dimensions);
    getAllWithin(x, neighbors, distance);
    // Important: upwrap the position vector so x destructing does not try to delete the position.
    x.unwrap();
    // Reset exclude flag.
    _exclude = -1;
  }

  void Domain::getAllWithin(Vec X, vector<int>& neighbors, RealType distance) {
    // Default distance
    if (distance<0) distance = 2*max_small_sigma;
    // Clear our the vector.
    neighbors.clear();

    // Set up
    vector<int> tuple1(sim_dimensions), tuple2(sim_dimensions), cell_index(sim_dimensions), center(sim_dimensions), search_dims(sim_dimensions);
    vector<RealType> dX(sim_dimensions);
    auto x = simData->X();
    auto type = simData->Type();

    // Calculate sweep "radius"
    int prod = 1;
    for (int d=0; d<sim_dimensions; ++d) {
      center[d] = 1 + static_cast<int>(ceil(distance*inverseW[d]));
      // Search dimensions can't be so large that any cells are searched more than once.
      search_dims[d] = min(2*center[d]+1, dims[d]);
      // Correct center based on actual search dimensions.
      center[d] = static_cast<int>(ceil(search_dims[d]/2.));
      prod *= search_dims[d];
    }

    // The tuple address of the cell the particle is in
    get_cell_index_tuple(X.data, cell_index);
    
    // Look in a hypercube.
    for (int j=0; j<prod; ++j) {
      // Turn j into a tuple
      getAddressCM(j, search_dims.data(), tuple1.data(), sim_dimensions);
      // Shift so it is a displacement
      subtractVec(tuple1.data(), center.data(), tuple2.data(), sim_dimensions);
      // Get the cell at the displacement from the particle's cell
      addVec(tuple2.data(), cell_index.data(), tuple1.data(), sim_dimensions);
      // If the cell is valid, look for particles
      if(correct_index(tuple1)) {
        // Get the linear address of the other cell
        int linear;
        tuple_to_linear(linear, tuple1);
        for (auto &id2 : cells[linear].particle_ids) {
          // Don't count the excluded particle.
          if (id2==_exclude) continue;
          // If the other particle is a larger particle, it will take care of this interaction
          Base::gflow->getDisplacement(X.data, Base::simData->X(id2), dX.data());
          RealType r = magnitudeVec(dX.data(), sim_dimensions);
          if (r<distance) neighbors.push_back(id2);
        }
      }
    }
  }

  void Domain::constructFor(int id1, bool insert) {
    // Find all neighor particles.
    auto x = simData->X();
    auto sg = simData->Sg();
    vector<int> neighbors;
    // Search for nearby particles.
    Vec X(sim_dimensions);
    X.wrap(x[id1], sim_dimensions);
    getAllWithin(X, neighbors, 2*sg[id1] + skin_depth);
    X.unwrap();
    // Build lists for the particle. List type 2.
    for (auto id2 : neighbors) {
      RealType r = gflow->getDistance(x[id1], x[id2]);
      if (r<sg[id1]+sg[id2]+skin_depth)
        pair_interaction(id1, id2, 2);
    }
    // If insert flag is set to true, insert into cells
    if (insert) add_to_cell(x[id1], id1);
  }

  void Domain::traversePairs(std::function<void(int, int, int, RealType, RealType, RealType)> body) {
    // Get the array of max cutoffs
    const vector<RealType> & max_cutoffs = forceMaster->getMaxCutoff();

    // Checks
    if (cutoff_grid==nullptr) return;

    // Maximum reasonable distance.
    RealType max_reasonable = sqr(0.5*simulation_bounds.wd(0));
    for (int d=1; d<sim_dimensions; ++d) {
      RealType mr = sqr(0.5*simulation_bounds.wd(d));
      if (mr<max_reasonable) max_reasonable = mr;
    }

    // Tuples
    vector<int> tuple1(sim_dimensions), tuple2(sim_dimensions), cell_index(sim_dimensions), center(sim_dimensions), search_dims(sim_dimensions);

    // Find potential neighbors
    auto sg = simData->Sg();
    auto x = simData->X();
    auto type = simData->Type();

    // Go through all the cells in the simulation.
    for (const auto &c : cells) {
      for (auto p=c.particle_ids.begin(); p!=c.particle_ids.end(); ++p) {
        // The id of the particle
        int id1 = *p;
        if (type[id1]<0) continue;
        RealType *cutoffs_id1 = cutoff_grid[type[id1]];
        RealType sigma1 = sg[id1]*max_cutoffs[type[id1]];

        // If sigma is <= than max_small_sigma, only look through cell stencil
        if (sigma1<=max_small_sigma) {
          // All other particles in the same sector
          for (auto q = p+1; q!=c.particle_ids.end(); ++q) {
            int id2 = *q;
            // If the other particle is a large particle, it will take care of this interaction
            if (type[id1]<0 || sg[id2]*max_cutoffs[type[id2]]>max_small_sigma) continue;
            // Look for distance between particles
            RealType r2 = getDistanceSqrNoWrap(x[id1], x[id2], sim_dimensions);
            if (r2 < sqr((sg[id1] + sg[id2])*cutoffs_id1[type[id2]] + skin_depth))
              body(id1, id2, 0, sg[id1], sg[id2], r2);
          }
          // Seach through list of adjacent cells
          for (const auto &d : c.adjacent)
            for (const auto id2 : d->particle_ids) {
              // If the other particle is a large particle, it will take care of this interaction
              if (type[id2]<0 || sg[id2]*max_cutoffs[type[id2]]>max_small_sigma) continue;
              // Look for distance between particles
              RealType r2 = getDistanceSqrNoWrap(x[id1], x[id2], sim_dimensions);
              if (r2 < sqr((sg[id1] + sg[id2])*cutoffs_id1[type[id2]] + skin_depth))
                body(id1, id2, 0, sg[id1], sg[id2], r2);
              else if (r2>max_reasonable)
                body(id1, id2, 1, sg[id1], sg[id2], r2);                
            }
        }
        
        // If sigma is > max_small_sigma, we have to look through more cells
        else {
          // Calculate sweep "radius"
          RealType search_width = 2*sigma1+skin_depth;
          int prod = 1;
          for (int d=0; d<sim_dimensions; ++d) {
            center[d] = static_cast<int>(ceil(search_width/widths[d]));
            // Search dimensions can't be so large that any cells are searched more than once.
            // Note: this needs to be 2*(center[d]+1) not 2*center[d] + 1. This caused an error before.
            search_dims[d] = min(2*(center[d]+1), dims[d]);
            // Correct center based on actual search dimensions.
            center[d] = static_cast<int>(ceil(search_dims[d]/2.));
            prod *= search_dims[d];
          }

          // The tuple address of the cell the particle is in
          get_cell_index_tuple(x[id1], cell_index);
          // Look in a hypercube.
          for (int j=0; j<prod; ++j) {
            // Turn j into a tuple
            getAddressCM(j, search_dims.data(), tuple1.data(), sim_dimensions);
            // Shift so it is a displacement
            subtractVec(tuple1.data(), center.data(), tuple2.data(), sim_dimensions);
            // Get the cell at the displacement from the particle's cell
            addVec(tuple2.data(), cell_index.data(), tuple1.data(), sim_dimensions);
            // If the cell is valid, look for particles
            if(correct_index(tuple1)) {
              // Get the linear address of the other cell
              int linear;
              tuple_to_linear(linear, tuple1);
              for (auto &id2 : cells[linear].particle_ids) {
                // If the other particle is a larger particle, it will take care of this interaction
                if (id1==id2 || sg[id2]>sg[id1]) continue;
                RealType r2 = getDistanceSqrNoWrap(x[id1], x[id2], sim_dimensions);
                if (r2 < sqr((sigma1 + sg[id2])*cutoffs_id1[type[id2]] + skin_depth))
                  body(id1, id2, 0, sg[id1], sg[id2], r2);
                else if (max_reasonable<r2)
                  body(id1, id2, 1, sg[id1], sg[id2], r2);
              }
            }
          } 
        }
      }
    }
  }

  inline void Domain::structure_updates() {
    // Clear cells
    for (auto &c : cells) c.particle_ids.clear();
    // We should have just done a particle removal, so we can use number, not size (since all arrays are compressed)
    int number = Base::simData->number();
    auto x = simData->X();
    auto type = Base::simData->Type();
    // Bin all the particles
    for (int i=0; i<number; ++i)
      if (0<=type[i] && forceMaster->typeInteracts(type[i])) 
        add_to_cell(x[i], i);
  }

  void Domain::calculate_domain_cell_dimensions() {
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

  inline void Domain::create_cells() {
    // --- Create the cells
    // Get the total number of cells - The dims MUST be set first.
    const int size = getNumCells();
    cells = vector<Cell>(size);

    // Holder for tuple index
    vector<int> tuple1(sim_dimensions), tuple2(sim_dimensions);
    
    // --- Create a neighborhood stencil to help us find adjacent cells.
    RealType min_small_cutoff = 2*max_small_sigma + skin_depth;
    vector<vector<int> > neighbor_indices;
    vector<int> little_dims(sim_dimensions), center(sim_dimensions);
    int volume = 1;
    for (int d=0; d<sim_dimensions; ++d) {
      int sweep = ceil(min_small_cutoff/widths[d]);
      little_dims[d] = 2*sweep+1;
      center[d] = sweep;
      volume *= little_dims[d];
    }

    // Create stencil
    for (int i=0; i<floor(0.5*volume); ++i) {
      getAddressCM(i, little_dims.data(), tuple1.data(), sim_dimensions);
      subtractVec(tuple1.data(), center.data(), tuple2.data(), sim_dimensions);
      int dr2 = dotVec<RealType, int, RealType>(tuple2.data(), widths.data(), sim_dimensions);
      if (dr2<sqr(min_small_cutoff)) neighbor_indices.push_back(tuple2);
    }

    // At this point, the vector neighor_indices now contains all the relative indices of neighbors cells.
    // Some "neighbors" might be out of bounds, but correct_index will take care of this.

    // --- Assign cell types, adjacent cells, set cell bounds
    for (int c=0; c<size; ++c) {
      linear_to_tuple(c, tuple1);

      // Set cell neighbors
      for (auto neigh : neighbor_indices) {
        // Get the potential neighbor cell's index using the stencil.
        addVec(tuple1.data(), neigh.data(), tuple2.data(), sim_dimensions);
        // The index may fall out of bounds.
        if (correct_index(tuple2, true)) {
          int linear;
          tuple_to_linear(linear, tuple2);
          cells[c].adjacent.push_back(&cells[linear]);
        }
      }
    }
  }

  inline bool Domain::correct_index(vector<int>& index, bool wrap) {
    bool good_index = true;
    const BCFlag *bcs = gflow->getBCs();
    for (int d=0; d<sim_dimensions; ++d) {
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

  inline void Domain::add_to_cell(const RealType *x, int id) {
    int linear = get_cell_index(Base::simData->X(id));
    // Stores the *local* id of the particle
    cells[linear].particle_ids.push_back(id);
  }

}
