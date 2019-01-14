#include "domain.hpp"
// Other files
#include "../gflow.hpp"
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../interactionhandlers/verletlist.hpp"
#include "../base/forcemaster.hpp"

#include "../base/interaction.hpp"
#include "../base/interactionhandler.hpp"

namespace GFlowSimulation {

  Domain::Domain(GFlow *gflow) : DomainBase(gflow) {
    border_type_up = new int[sim_dimensions];
    border_type_down = new int[sim_dimensions];
    products = new int[sim_dimensions+1];
    // Initialize values
    for (int d=0; d<sim_dimensions; ++d) {
      border_type_up[d]   = 0;
      border_type_down[d] = 0;
    }
  };

  Domain::~Domain() {
    // Delete this object's data
    if (border_type_up)   delete [] border_type_up;
    if (border_type_down) delete [] border_type_down;
    if (products)         delete [] products;
    border_type_up = nullptr;
    border_type_down = nullptr;
    products = nullptr;
  }

  void Domain::initialize() {
    //! @todo Domains are generally initialized several times, though they doesn't need to be. Find a way to check whether we
    //! actually need to initialize the domain.

    // Common initialization tasks
    Base::initialize();

    // Get the simulation bounds from gflow
    bounds = Base::gflow->getBounds();
    // Get the bounds from the gflow object - for now assumes this is the only domain, so bounds==domain_bounds
    domain_bounds = gflow->getBounds();

    // If bounds are unset, then don't make sectors. We cannot initialize if simdata is null
    if (domain_bounds.vol()<=0 || simData==nullptr) return; 

    // Assign border types - do this before creating cells.
    assign_border_types();

    // Calculate the maxiumu "small sigma"
    calculate_max_small_sigma();

    // Use max_small_sigma
    cutoff = minCutoff = 2*max_small_sigma+skin_depth;

    calculate_domain_cell_dimensions();
    
    // Initialize products array
    products[sim_dimensions] = 1;
    for (int d=sim_dimensions-1; d>=0; --d)
      products[d] = dims[d]*products[d+1];

    //! @todo Processors might need to communicate with one another about what they chose at this point

    // Create the cells
    create_cells();

    // Construct the interaction handlers for the forces
    number = 0; // To make sure we do a full clear, fill when we construct.
    construct();

    // The domain has been initialized
    initialized = true;
  }

  void Domain::pre_integrate() {
    // Reset time points
    lastCheck  = -1.;
    lastUpdate = -1.;
    updateDelay = 1.0e-4;
  }

  void Domain::exchange_particles() {
    //! @todo Implement this.
  }

  void Domain::getAllWithin(int, RealType, vector<int>&) {
    //! @todo Implement this.
  }

  void Domain::removeOverlapping(RealType factor) {
    // Domain base common tasks
    DomainBase::construct();
    
    // Update particles in the cells
    update_cells();

    RealType max_reasonable = sqr(0.9*bounds.wd(0));

    // Tuples
    int *tuple1 = new int[sim_dimensions], *tuple2 = new int[sim_dimensions];
    int *cell_index = new int[sim_dimensions], *center = new int[sim_dimensions];

    // Find potential neighbors
    RealType *sg = Base::simData->Sg();
    RealType **x = Base::simData->X();

    // Get the boundary conditions
    const BCFlag *bcs = gflow->getBCs();
    RealType *dX = new RealType[sim_dimensions];
    int *search_dims = new int[sim_dimensions];

    // Go through all the cells
    for (const auto &c : cells) {
      for (auto p=c.particle_ids.begin(); p!=c.particle_ids.end(); ++p) {
        int id1 = *p;
        // If sigma is <= than min_small_sigma, only look through cell stencil
        if (sg[id1]<=max_small_sigma) {
          // All other particles in the same sector
          auto q = p;
          ++q;
          for (; q!=c.particle_ids.end(); ++q) {
            int id2 = *q;
            // If the other particle is a large particle, it will take care of this interaction
            if (sg[id2]>max_small_sigma) continue;
            // Get distance between particles
            subtractVec(x[id1], x[id2], dX, sim_dimensions);
            RealType r = magnitudeVec(dX, sim_dimensions);
            RealType overlap = sg[id1] + sg[id2] - r;
            if (overlap/min(sg[id1], sg[id2]) > factor)
              Base::simData->markForRemoval(sg[id1]>sg[id2] ? id2 : id1);
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
                // If the other particle is a larger particle, it will take care of this interaction
                if (id1==id2 || sg[id2]>sg[id1] || (sg[id1]==sg[id2] && id1<id2)) continue; // IF TWO PARTICLES ARE THE SAME SIZE, ERROR
                Base::gflow->getDisplacement(Base::simData->X(id1), Base::simData->X(id2), dX);
                RealType r = magnitudeVec(dX, sim_dimensions);
                RealType overlap = sg[id1] + sg[id2] - r;
                if (overlap/sg[id2] > factor)
                  simData->markForRemoval(id2);
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

    // Update particles in the cells
    update_cells();

    // Set a "maximum reasonable distance"
    RealType max_reasonable = sqr(0.9*bounds.wd(0));

    if (gflow->getInteractions().empty()) return;

    // A tuple
    int *tuple1 = new int[sim_dimensions], *tuple2 = new int[sim_dimensions];
    int *cell_index = new int[sim_dimensions], *center = new int[sim_dimensions];

    // Find potential neighbors
    RealType *sg = Base::simData->Sg();
    RealType **x = Base::simData->X();
    // Get the boundary conditions
    const BCFlag *bcs = gflow->getBCs();
    int *search_dims = new int[sim_dimensions];
    for (const auto &c : cells) {
      for (auto p=c.particle_ids.begin(); p!=c.particle_ids.end(); ++p) {
        // The id of the particle
        int id1 = *p;
        // If sigma is <= than min_small_sigma, only look through cell stencil
        if (sg[id1]<=max_small_sigma) {
          // All other particles in the same sector
          auto q = p;
          ++q;
          for (; q!=c.particle_ids.end(); ++q) {
            int id2 = *q;
            RealType r2 = getDistanceSqrNoWrap(x[id1], x[id2], sim_dimensions);
            if (r2 < sqr(sg[id1] + sg[id2] + skin_depth))
              pair_interaction(id1, id2);
          }
          // Seach through list of adjacent cells
          for (const auto &d : c.adjacent)
            for (const auto id2 : d->particle_ids) {
              // If the other particle is a large particle, it will take care of this interaction
              if (sg[id2]>max_small_sigma) continue;
              // Look for distance between particles
              RealType r2 = getDistanceSqrNoWrap(x[id1], x[id2], sim_dimensions);
              if (r2 < sqr(sg[id1] + sg[id2] + skin_depth) || max_reasonable<r2)
                pair_interaction(id1, id2);
            }
        }
        
        // If sigma is > min_small_sigma, we have to look through more cells
        else {
          // Calculate sweep "radius"
          RealType search_width = 2*sg[id1]+skin_depth;
          int prod = 1;
          for (int d=0; d<sim_dimensions; ++d) {
            center[d] = static_cast<int>(ceil(search_width/widths[d]));
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
                // If the other particle is a larger particle, it will take care of this interaction
                if (id1==id2 || sg[id2]>sg[id1]) continue;
                RealType r2 = getDistanceSqrNoWrap(x[id1], x[id2], sim_dimensions);
                if (r2 < sqr(sg[id1] + sg[id2] + skin_depth) || max_reasonable<r2)
                  pair_interaction(id1, id2);
              }
            }
          }
        }
      }
    }
    // Close all
    Base::forceMaster->close();
    // Clean up
    delete [] tuple1;
    delete [] tuple2;
    delete [] cell_index;
    delete [] center;
    delete [] search_dims;
  }

  void Domain::setCellSize(RealType) {
    // @todo Implement this.
  }

  inline void Domain::update_cells() {
    clear_cells();
    fill_cells();
    // Record how many particles there were
    number = simData->number();
  }

  inline void Domain::assign_border_types() {
    const BCFlag *bcs = Base::gflow->getBCs(); // Get the boundary condition flags
    //! \todo Use topology object to determine border types.
    //! For now, just use halo cells whenever possible.
    for (int d=0; d<sim_dimensions; ++d) {
      if (bcs[d]==BCFlag::WRAP) {
        border_type_up[d] = 1;
        border_type_down[d] = 1;
      }
      else {
        border_type_up[d] = 0;
        border_type_down[d] = 0;
      }
    }
  }

  inline void Domain::calculate_domain_cell_dimensions() {
    for (int d=0; d<sim_dimensions; ++d) {
      dims[d] = static_cast<int>(domain_bounds.wd(d)/cutoff);
      // Check that the bounds are good
      if (dims[d]<=0) throw BadBounds();
      widths[d] = domain_bounds.wd(d)/dims[d];
      inverseW[d] = 1./widths[d];
    }
  }

  inline void Domain::create_cells() {
    // --- Create the cells
    // Get the total number of cells - The dims MUST be set first.
    const int size = getNumCells();
    cells = vector<Cell>(size, Cell(sim_dimensions));

    // Holder for tuple index
    int *tuple1 = new int[sim_dimensions], *tuple2 = new int[sim_dimensions];
    // --- Create a neighborhood stencil to help us find adjacent cells.
    int sweep = ceil(minCutoff/min(widths, sim_dimensions));
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
    for (auto &c : cells) c.particle_ids.clear();
  }

  inline void Domain::fill_cells() {
    // We should have just done a particle removal, so we can use number, not size (since all arrays are compressed)
    RealType **x = simData->X();
    int number = Base::simData->number();
    int *tuple = new int[sim_dimensions], linear;

    // Bin all the particles
    for (int i=0; i<number; ++i) {
      linear = get_cell_index(x[i]);
      // Stores the *local* id of the particle
      cells[linear].particle_ids.push_back(i);
    }
    
    // do_halo_assignment();

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

  inline bool Domain::correct_index(int *tuple, bool wrap) {
    bool good_index = true;
    const BCFlag *bcs = gflow->getBCs();
    for (int d=0; d<sim_dimensions; ++d) {
      if (tuple[d]>=dims[d]) {
        if (bcs[d]==BCFlag::WRAP && wrap)
          tuple[d] -= dims[d];
        else good_index = false;
      }
      if (tuple[d]<0) {
        if (bcs[d]==BCFlag::WRAP && wrap) tuple[d] += dims[d];
        else good_index = false;
      }
    }
    return good_index;
  }

  inline void Domain::get_cell_index_tuple(const RealType *x, int *index) {
    for (int d=0; d<sim_dimensions; ++d) {
      index[d] = static_cast<int>((x[d] - domain_bounds.min[d])*inverseW[d]);

      // Even when wrapping, rounding errors or drifting can cause index to be too large or small.
      if (index[d]>=dims[d]) index[d] = dims[d]-1;
      else if (index[d]<0)   index[d] = 0;
    }
  }

  int Domain::get_cell_index(const RealType *x) {
    int linear = 0;
    for (int d=0; d<sim_dimensions; ++d) {
      RealType dth = static_cast<int>((x[d] - domain_bounds.min[d])*inverseW[d]);
      if (dth>=dims[d]) dth = dims[d]-1;
      else if (dth<0)   dth = 0;
      linear += dth*products[d+1];
    }
    // Return the index
    return linear;
  }

  inline void Domain::add_to_cell(const RealType *x, int id) {
    int linear = get_cell_index(Base::simData->X(id));
    // Stores the *local* id of the particle
    cells[linear].particle_ids.push_back(id);
  }

  inline void Domain::do_halo_assignment() {
    // TWO DIMENSIONS
    int count = 0, linear = 0, *tuple = new int[sim_dimensions];

    vector<int> plus_list, minus_list;
    RealType *dR = new RealType[sim_dimensions];
    
    if (border_type_up[0]) { // Wrap in X direction
      // Bottom
      for (int i=0; i<dims[0]; ++i) {
        tuple[0] = 0; tuple[1] = i;
        tuple_to_linear(linear, tuple);
        // Create halo particle
        for (auto id : cells[linear].particle_ids) {
          ++count;
          plus_list.push_back(id);
        }
      }
      // Top
      for (int i=0; i<dims[0]; ++i) {
        tuple[0] = dims[1]-1; tuple[1] = i;
        tuple_to_linear(linear, tuple);
        for (auto id : cells[linear].particle_ids) {
          ++count;
          minus_list.push_back(id);
        }
      }
      // Add particles to simdata and cells
      dR[0] = domain_bounds.wd(0); dR[1] = 0;
      halo_list_add(plus_list,  dR);
      dR[0] = -domain_bounds.wd(0); dR[1] = 0;
      halo_list_add(minus_list, dR);
      plus_list.clear(); minus_list.clear();
    }

    /*
    if (border_type_up[1]) {
      // Bottom
      for (int i=0; i<dims[1]; ++i) {
        tuple[0] = i; tuple[1] = 0;
        tuple_to_linear(linear, tuple);
        for (auto id : cells[linear].particle_ids) {
          ++count;
          plus_list.push_back(id);
        }
      }
      // Top
      for (int i=0; i<dims[1]; ++i) {
        tuple[0] = i; tuple[1] = dims[0]-1;
        tuple_to_linear(linear, tuple);
        for (auto id : cells[linear].particle_ids) {
          ++count;
          minus_list.push_back(id);
        }
      }
      // Add particles to simdata and cells
      dR[0] = 0; dR[1] = domain_bounds.wd(1);
      halo_list_add(plus_list,  dR);
      dR[0] = 0; dR[1] = -domain_bounds.wd(1);
      halo_list_add(minus_list, dR);
      plus_list.clear(); minus_list.clear();
    }
    */

    //! \todo Check edge cells for particles that should produce ghost particles
    // ---

    // Clean up
    delete [] tuple;
  }

  inline void Domain::halo_list_add(const vector<int>& id_list, RealType *displacement) {
    for (auto id : id_list) {
      Base::simData->create_halo_of(id, displacement);
      int id2 = Base::simData->size()-1; // The id of the halo particle
      add_to_cell(Base::simData->X(id2), id2); // Add the halo particle to the cells
    }
  }

  void Domain::calculate_max_small_sigma() {
    // Find average sigma
    RealType sigma = 0, max_sigma = 0;
    for (int n=0; n<Base::simData->size(); ++n) {
      if (Base::simData->Type(n)<0) continue;
      RealType s = Base::simData->Sg(n);
      sigma += s;
      if (s>max_sigma) max_sigma = s;
    }
    sigma /= Base::simData->number();
    // Threshhold sigma is between average and maximum sigma
    RealType threshold = 0.5*(sigma + max_sigma), max_under = sigma;
    if (threshold!=sigma) {
      for (int n=0; n<Base::simData->size(); ++n) {
        if (Base::simData->Type(n)<0) continue;
        RealType s = Base::simData->Sg(n);
        if (s<threshold && max_under<s) max_under = s;
      }
    }
    max_small_sigma = 1.025*max_under;
  }

}
