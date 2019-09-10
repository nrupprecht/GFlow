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

  Domain::Domain(GFlow *gflow) : DomainBase(gflow) {
    // Allocate arrays
    border_type_up = new int[sim_dimensions];
    border_type_down = new int[sim_dimensions];
    dim_shift_up = new int[sim_dimensions];
    dim_shift_down = new int[sim_dimensions];
    products = new int[sim_dimensions+1];
    // Initialize to zero
    zeroVec(border_type_up, sim_dimensions);
    zeroVec(border_type_down, sim_dimensions);
    zeroVec(dim_shift_up, sim_dimensions);
    zeroVec(dim_shift_down, sim_dimensions);
    zeroVec(products, sim_dimensions);
  };

  Domain::~Domain() {
    // Delete this object's data.
    if (border_type_up)   delete [] border_type_up;
    if (border_type_down) delete [] border_type_down;
    if (dim_shift_up)     delete [] dim_shift_up;
    if (dim_shift_down)   delete [] dim_shift_down;
    if (products)         delete [] products;
    // Null pointers, just in case.
    border_type_up = nullptr;
    border_type_down = nullptr;
    dim_shift_up = nullptr;
    dim_shift_down = nullptr;
    products = nullptr;
  }

  void Domain::initialize() {
    //! \todo Domains are generally initialized several times, though they doesn't need to be. Find a way to check whether we
    //! actually need to initialize the domain.

    // Common initialization tasks. Calculates max small sigma.
    DomainBase::initialize();

    // If bounds are unset, then don't make sectors. We cannot initialize if simdata is null
    if (process_bounds.vol()<=0 || isnan(process_bounds.vol()) || simData==nullptr || simData->size()==0) return;

    // Assign border types - do this before creating cells.
    assign_border_types();

    // Calculate skin depth
    calculate_skin_depth();

    // Use max_small_sigma. It is important that all processors share a consistent min_small_cutoff.
    // This can be achieved by sharing a max_small_sigma, and skin_depth.
    target_cell_size = min_small_cutoff = 2*max_small_sigma+skin_depth;

    // Calculate cell grid data
    calculate_domain_cell_dimensions();

    // Create the cells
    create_cells();

    // Construct the interaction handlers for the forces
    construct();

    // The domain has been initialized
    initialized = true;
  }

  void Domain::getAllWithin(int id1, vector<int>& neighbors, RealType distance) {
    // Default distance
    if (distance<0) distance = 2*max_small_sigma;

    // Set up
    int *tuple1 = new int[sim_dimensions], *tuple2 = new int[sim_dimensions];
    int *cell_index = new int[sim_dimensions], *center = new int[sim_dimensions];
    
    RealType **x = Base::simData->X();
    int *type = Base::simData->Type();

    // Get the boundary conditions
    RealType *dX = new RealType[sim_dimensions];
    int *search_dims = new int[sim_dimensions];

    // Calculate sweep "radius"
    int prod = 1;
    for (int d=0; d<sim_dimensions; ++d) {
      center[d] = static_cast<int>(ceil(distance/widths[d]));
      search_dims[d] = 2*center[d]+1;
      prod *= search_dims[d];
    }

    // The tuple address of the cell the particle is in
    get_cell_index_tuple(x[id1], cell_index);
    
    // Look in a hypercube.
    for (int j=0; j<prod; ++j) {
      // Turn j into a tuple
      getAddressCM(j, search_dims, tuple1, sim_dimensions);
      // Shift so it is a displacement
      subtractVec(tuple1, center, tuple2, sim_dimensions);
      // Get the cell at the displacement from the particle's cell
      addVec(tuple2, cell_index, tuple1, sim_dimensions);
      // If the cell is valid, look for particles
      if(correct_index(tuple1)) {
        // Get the linear address of the other cell
        int linear;
        tuple_to_linear(linear, tuple1);
        for (auto &id2 : cells[linear].particle_ids) {
          // Don't count yourself.
          if (id1==id2) continue;
          // If the other particle is a larger particle, it will take care of this interaction
          Base::gflow->getDisplacement(Base::simData->X(id1), Base::simData->X(id2), dX);
          RealType r = magnitudeVec(dX, sim_dimensions);
          if (r<distance) neighbors.push_back(id2);
        }
      }
    }

    // Clean up
    delete [] tuple1;
    delete [] tuple2;
    delete [] cell_index;
    delete [] center;
    delete [] dX;
    delete [] search_dims;
  }

  void Domain::getAllWithin(Vec X, vector<int>& neighbors, RealType distance) {
    // Default distance
    if (distance<0) distance = 2*max_small_sigma;

    // Set up
    int *tuple1 = new int[sim_dimensions], *tuple2 = new int[sim_dimensions];
    int *cell_index = new int[sim_dimensions], *center = new int[sim_dimensions];
    
    RealType **x = Base::simData->X();
    int *type = Base::simData->Type();

    // Get the boundary conditions
    RealType *dX = new RealType[sim_dimensions];
    int *search_dims = new int[sim_dimensions];


    // Calculate sweep "radius"
    int prod = 1;
    for (int d=0; d<sim_dimensions; ++d) {
      center[d] = static_cast<int>(ceil(distance/widths[d]));
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
      getAddressCM(j, search_dims, tuple1, sim_dimensions);
      // Shift so it is a displacement
      subtractVec(tuple1, center, tuple2, sim_dimensions);
      // Get the cell at the displacement from the particle's cell
      addVec(tuple2, cell_index, tuple1, sim_dimensions);
      // If the cell is valid, look for particles
      if(correct_index(tuple1)) {
        // Get the linear address of the other cell
        int linear;
        tuple_to_linear(linear, tuple1);
        for (auto &id2 : cells[linear].particle_ids) {
          // If the other particle is a larger particle, it will take care of this interaction
          Base::gflow->getDisplacement(X.data, Base::simData->X(id2), dX);
          RealType r = magnitudeVec(dX, sim_dimensions);
          if (r<distance) neighbors.push_back(id2);
        }
      }
    }

    // Clean up
    delete [] tuple1;
    delete [] tuple2;
    delete [] cell_index;
    delete [] center;
    delete [] dX;
    delete [] search_dims;
  }

  void Domain::removeOverlapping(RealType factor) {
    // Update particles in the cells
    update_cells();

    // Tuples
    int *tuple1 = new int[sim_dimensions], *tuple2 = new int[sim_dimensions];
    int *cell_index = new int[sim_dimensions], *center = new int[sim_dimensions];
    
    // Find potential neighbors
    RealType *sg = Base::simData->Sg();
    RealType **x = Base::simData->X();
    int *type = Base::simData->Type();

    // Get the boundary conditions
    const BCFlag *bcs = gflow->getBCs();
    RealType *dX = new RealType[sim_dimensions];
    int *search_dims = new int[sim_dimensions];

    // Go through all the cells
    for (const auto &c : cells) {
      for (auto p=c.particle_ids.begin(); p!=c.particle_ids.end(); ++p) {
        int id1 = *p;
        if (type[id1]<0) continue;
        // If sigma is <= than min_small_sigma, only look through cell stencil
        if (sg[id1]<=max_small_sigma) {
          // All other particles in the same sector
          auto q = p;
          ++q;
          for (; q!=c.particle_ids.end(); ++q) {
            int id2 = *q;
            // If the other particle is a large particle, it will take care of this interaction
            if (sg[id2]>max_small_sigma || type[id2]<0) continue;
            // Get distance between particles
            subtractVec(x[id1], x[id2], dX, sim_dimensions);
            RealType r = magnitudeVec(dX, sim_dimensions);
            RealType overlap = sg[id1] + sg[id2] - r;
            if (overlap/min(sg[id1], sg[id2]) > factor) {
              Base::simData->markForRemoval(sg[id1]>sg[id2] ? id2 : id1);
            }
          }

          // Seach through list of adjacent cells
          for (const auto &d : c.adjacent)
            for (const auto id2 : d->particle_ids) {
              // If the other particle is a large particle, it will take care of this interaction
              if (sg[id2]>max_small_sigma) continue;
              // Look for distance between particles
              Base::gflow->getDisplacement(Base::simData->X(id1), Base::simData->X(id2), dX);
              RealType r = magnitudeVec(dX, sim_dimensions);

              RealType overlap = sg[id1] + sg[id2] - r;
              if (overlap/min(sg[id1], sg[id2]) > factor) 
                Base::simData->markForRemoval(sg[id1]>sg[id2] ? id2 : id1);
            }
        }
        // If sigma is > min_small_sigma, we have to look through more cells
        else {
          // Calculate sweep "radius"
          RealType search_width = 2*sg[id1]+skin_depth;
          int prod = 1;
          for (int d=0; d<sim_dimensions; ++d) {
            center[d] = static_cast<int>(ceil(search_width/widths[d]));
            // Search dimensions can't be so large that any cells are searched more than once.
            search_dims[d] = min(2*center[d]+1, dims[d]);
            // Correct center based on actual search dimensions.
            center[d] = static_cast<int>(ceil(search_dims[d]/2.));
            prod *= search_dims[d];
          }

          // The tuple address of the cell the particle is in
          get_cell_index_tuple(x[id1], cell_index);
          
          // Look in a hypercube.
          for (int j=0; j<prod; ++j) {
            // Turn j into a tuple
            getAddressCM(j, search_dims, tuple1, sim_dimensions);
            // Shift so it is a displacement
            subtractVec(tuple1, center, tuple2, sim_dimensions);
            // Get the cell at the displacement from the particle's cell
            addVec(tuple2, cell_index, tuple1, sim_dimensions);

            // If the cell is valid, look for particles
            if(correct_index(tuple1)) {
              // Get the linear address of the other cell
              int linear;
              tuple_to_linear(linear, tuple1);
              for (auto &id2 : cells[linear].particle_ids) {
                // If the other particle is a larger particle, it will take care of this interaction
                if (id1==id2 || sg[id2]>sg[id1] || (sg[id1]==sg[id2] && id1<id2)) continue; // IF TWO PARTICLES ARE THE SAME SIZE, ERROR
                gflow->getDisplacement(x[id1], x[id2], dX);
                RealType r = magnitudeVec(dX, sim_dimensions);

                RealType overlap = sg[id1] + sg[id2] - r;
                if (overlap/min(sg[id1], sg[id2]) > factor) {
                  simData->markForRemoval(id2);
                }
              }
            }
          }
          
        }
      }
    }

    // Remove particles
    Base::simData->doParticleRemoval();

    // Clean up
    delete [] tuple1;
    delete [] tuple2;
    delete [] cell_index;
    delete [] center;
    delete [] dX;
    delete [] search_dims;
  }

  void Domain::construct() {
    // Domain base common tasks
    DomainBase::construct();
    
    // Start timer.
    start_timer();

    // Update particles in the cells
    update_cells();

    // Maximum reasonable distance.
    RealType max_reasonable = sqr(0.9*simulation_bounds.wd(0));
    for (int d=1; d<sim_dimensions; ++d) {
      RealType mr = sqr(0.9*simulation_bounds.wd(d));
      if (mr<max_reasonable) max_reasonable = mr;
    }

    // If there are no interaction, we don't need to make any verlet lists
    if (gflow->getInteractions().empty()) return;
    // If there is no force master, return
    if (forceMaster==nullptr) return;
    // Make sure force master has interaction array set up
    forceMaster->initialize_does_interact();
    // Get the array of max cutoffs
    const vector<RealType> & max_cutoffs = forceMaster->getMaxCutoff();

    // Tuples
    vector<int> tuple1(sim_dimensions), tuple2(sim_dimensions);
    //int *tuple1 = new int[sim_dimensions], *tuple2 = new int[sim_dimensions];
    int *cell_index = new int[sim_dimensions], *center = new int[sim_dimensions];
    Vec dx(sim_dimensions);

    // Find potential neighbors
    RealType *sg = simData->Sg();
    RealType **x = simData->X();
    int    *type = simData->Type();
    // Get the boundary conditions
    const BCFlag *bcs = gflow->getBCs();
    vector<int> search_dims(sim_dimensions, 0);

    for (const auto &c : cells) {
      for (auto p=c.particle_ids.begin(); p!=c.particle_ids.end(); ++p) {
        // The id of the particle
        int id1 = *p;
        RealType *cutoffs_id1 = cutoff_grid[type[id1]];
        RealType sigma1 = sg[id1]*max_cutoffs[type[id1]];

        // If sigma is <= than max_small_sigma, only look through cell stencil
        if (sigma1<=max_small_sigma) {
          // All other particles in the same sector
          auto q = p;
          ++q;
          for (; q!=c.particle_ids.end(); ++q) {
            int id2 = *q;
            // If the other particle is a large particle, it will take care of this interaction
            if (sg[id2]*max_cutoffs[type[id2]]>max_small_sigma) continue;
            // Look for distance between particles
            RealType r2 = getDistanceSqrNoWrap(x[id1], x[id2], sim_dimensions);
            if (r2 < sqr((sg[id1] + sg[id2])*cutoffs_id1[type[id2]] + skin_depth)) {
              pair_interaction(id1, id2);
            }
          }
          // Seach through list of adjacent cells
          for (const auto &d : c.adjacent)
            for (const auto id2 : d->particle_ids) {
              // If the other particle is a large particle, it will take care of this interaction
              if (sg[id2]*max_cutoffs[type[id2]]>max_small_sigma) continue;
              // Look for distance between particles
              RealType r2 = getDistanceSqrNoWrap(x[id1], x[id2], sim_dimensions);
              if (r2 < sqr((sg[id1] + sg[id2])*cutoffs_id1[type[id2]] + skin_depth))
                pair_interaction(id1, id2);
              else if (r2>max_reasonable) pair_interaction(id1, id2, 1);
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
            subtractVec(tuple1.data(), center, tuple2.data(), sim_dimensions);
            // Get the cell at the displacement from the particle's cell
            addVec(tuple2.data(), cell_index, tuple1.data(), sim_dimensions);
            // If the cell is valid, look for particles
            if(correct_index(tuple1.data())) {
              // Get the linear address of the other cell
              int linear;
              tuple_to_linear(linear, tuple1.data());
              for (auto &id2 : cells[linear].particle_ids) {
                // If the other particle is a larger particle, it will take care of this interaction
                if (id1==id2 || sg[id2]>sg[id1]) continue;
                RealType r2 = getDistanceSqrNoWrap(x[id1], x[id2], sim_dimensions);
                if (r2 < sqr((sigma1 + sg[id2])*cutoffs_id1[type[id2]] + skin_depth))
                  pair_interaction(id1, id2);
                else if (max_reasonable<r2) pair_interaction(id1, id2, 1);
              }
            }
          }
        }
      }
    }
    
    // Close all
    forceMaster->close();

    // Clean up
    delete [] cell_index;
    delete [] center;

    // Start timer.
    stop_timer();
  }

  void Domain::constructFor(int id, bool insert) {

  }

  void Domain::setCellSize(RealType) {
    // @todo Implement this.
  }

  inline void Domain::update_cells() {
    clear_cells();
    fill_cells();
  }


  void Domain::migrate_particles() {
    
  }

  void Domain::construct_halo_particles() {

  }

  void Domain::construct_ghost_particles() {

  }

  inline void Domain::assign_border_types() {
    // Get the boundary condition flags
    const BCFlag *bcs = Base::gflow->getBCs(); 
    //! \todo Use topology object to determine border types.
    //! For now, just use halo cells whenever possible.
    for (int d=0; d<sim_dimensions; ++d) {
      if (bcs[d]==BCFlag::WRAP) {
        // This processor takes up the whole bounds in this dimension, there are no processors above or below this one.
        if (process_bounds.wd(d)==simulation_bounds.wd(d)) {
          border_type_up[d]   = 0;
          border_type_down[d] = 0;
        }
        // There are one or more processors above and below this one.
        else {
          // If this processor expects ghost particles that are wrapped.
          if (process_bounds.max[d]==simulation_bounds.max[d]) border_type_up[d] = 2;
          // Ghost particles, but not wrapped ones.
          else border_type_up[d] = 1;

          // If this processor expects ghost particles that are wrapped.
          if (process_bounds.min[d]==simulation_bounds.min[d]) border_type_down[d] = 2;
          // Ghost particles, but not wrapped ones.
          else border_type_down[d] = 1;
        }
      }
      else {
        // Are there ghost particles?
        if (process_bounds.max[d]==simulation_bounds.max[d]) border_type_up[d] = 0;
        else border_type_up[d] = 1;
        // Are there ghost particles?
        if (process_bounds.min[d]==simulation_bounds.min[d]) border_type_down[d] = 0;
        else border_type_down[d] = 1;
      }
    }
  }

  inline void Domain::calculate_skin_depth() {
    RealType rho = simData->size() / process_bounds.vol();
    RealType candidate = inv_sphere_volume((target_list_size)/rho + 0.5*sphere_volume(max_small_sigma, sim_dimensions), sim_dimensions) - 2*max_small_sigma;
    skin_depth = max(static_cast<RealType>(0.5 * max_small_sigma), candidate);
    // Use the same skin depth on all processors - take the average.
    if (topology) {
      MPIObject::mpi_sum(skin_depth);
      skin_depth /= static_cast<RealType>(topology->getNumProc());
    }
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

    // Initialize products array
    calculate_product_array();
  }

  inline void Domain::calculate_product_array() {
    products[sim_dimensions] = 1;
    for (int d=sim_dimensions-1; d>=0; --d)
      products[d] = dims[d]*products[d+1];
  }

  inline void Domain::create_cells() {
    // --- Create the cells
    // Get the total number of cells - The dims MUST be set first.
    const int size = getNumCells();
    cells = vector<Cell>(size, Cell(sim_dimensions));

    // Holder for tuple index
    int *tuple1 = new int[sim_dimensions], *tuple2 = new int[sim_dimensions];
    
    // --- Create a neighborhood stencil to help us find adjacent cells.
    int sweep = ceil(min_small_cutoff/min(widths, sim_dimensions));
    vector<int*> neighbor_indices;
    int *little_dims = new int[sim_dimensions], *center = new int[sim_dimensions];
    for (int d=0; d<sim_dimensions; ++d) little_dims[d] = 2*sweep+1;
    for (int d=0; d<sim_dimensions; ++d) center[d]      = sweep;

    // Create stencil
    for (int i=0; i<floor(0.5*pow(2*sweep+1, sim_dimensions)); ++i) {
      getAddressCM(i, little_dims, tuple1, sim_dimensions);
      subtractVec(tuple1, center, tuple2, sim_dimensions);
      int dr2 = sqr(tuple2, sim_dimensions);
      if (dr2<sqr(sweep+1)) {
        int *new_index = new int[sim_dimensions];
        copyVec(tuple2, new_index, sim_dimensions);
        neighbor_indices.push_back(new_index);
      }
    }

    // At this point, the vector neighor_indices now contains all the relative indices of neighbors cells.
    // Some "neighbors" might be out of bounds, but correct_index will take care of this.

    // --- Assign cell types, adjacent cells, set cell bounds
    for (int c=0; c<size; ++c) {
      linear_to_tuple(c, tuple1);

      // Set cell neighbors
      for (auto n : neighbor_indices) {
        addVec(tuple1, n, tuple2, sim_dimensions);
        
        if (correct_index(tuple2, true)) {
          int linear;
          tuple_to_linear(linear, tuple2);
          cells[c].adjacent.push_back(&cells[linear]);
        }
      }
    }


    // Clean up 
    delete [] tuple1;
    delete [] tuple2;
    delete [] little_dims;
    delete [] center;
    for (auto &n : neighbor_indices) delete [] n;
  }

  inline void Domain::clear_cells() {
    for (auto &c : cells) 
      c.particle_ids.clear();
  }

  inline void Domain::fill_cells() {
    // We should have just done a particle removal, so we can use number, not size (since all arrays are compressed)
    RealType **x = simData->X();
    int number = Base::simData->number();
    int *type = Base::simData->Type();
    int *tuple = new int[sim_dimensions], linear;
    // Bin all the particles
    for (int i=0; i<number; ++i) {
      if (0<=type[i] && forceMaster->typeInteracts(type[i])) 
        add_to_cell(x[i], i);
    }

    // Clean up
    delete [] tuple;
  }

  // Turns a linear cell index into a (DIMENSIONS)-dimensional index
  inline void Domain::linear_to_tuple(int linear, int *tuple) {
    // Same as the get address function in
    getAddressCM(linear, dims, tuple, sim_dimensions); // We need to use the column major form
  }

  // Turns a (DIMENSIONS)-dimensional index into a linear cell index
  // This is column major form.
  inline void Domain::tuple_to_linear(int &linear, const int *tuple) {
    // Product lambda
    linear = 0;
    for (int d=0; d<sim_dimensions; ++d)
      linear += tuple[d]*products[d+1];
  }

  inline bool Domain::correct_index(int *index, bool wrap) {
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

  inline void Domain::get_cell_index_tuple(const RealType *x, int *index) {
    for (int d=0; d<sim_dimensions; ++d)
      index[d] = static_cast<int>((x[d] - process_bounds.min[d])*inverseW[d]) + dim_shift_down[d];
  }

  int Domain::get_cell_index(const RealType *x) {
    int linear = 0;
    for (int d=0; d<sim_dimensions; ++d) {
      RealType index = static_cast<int>((x[d] - process_bounds.min[d])*inverseW[d]) + dim_shift_down[d];
      if (index>=dims[d]) index = dims[d]-1;
      else if (index<0)   index = 0;
      linear += index*products[d+1];
    }
    // Return the index
    return linear;

    //return get_index<2>(x, process_bounds.min, inverseW, dim_shift_down, dims, products);
  }

  inline void Domain::add_to_cell(const RealType *x, int id) {
    int linear = get_cell_index(Base::simData->X(id));
    // Stores the *local* id of the particle
    cells[linear].particle_ids.push_back(id);
  }

}
