#include "domain_test.hpp"
// Other files
#include "../gflow.hpp"
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../interactionhandlers/verletlist.hpp"
#include "../base/forcemaster.hpp"

namespace GFlowSimulation {

  DomainTest::DomainTest(GFlow *gflow) : DomainBase(gflow) {};

  DomainTest::~DomainTest() {
    if (border_type_up)   delete [] border_type_up;
    if (border_type_down) delete [] border_type_down;
  }

  void DomainTest::initialize() {
    //! @todo Domains are generally initialized several times, though they doesn't need to be. Find a way to check whether we
    //! actually need to initialize the domain.

    // Common initialization tasks
    Base::initialize();

    // Get the simulation bounds from gflow
    bounds = Base::gflow->getBounds();
    // Get the bounds from the gflow object - for now assumes this is the only domain, so bounds==domain_bounds
    domain_bounds = gflow->getBounds();

    RealType target_width = (2*0.05+skin_depth);
    max_small_sigma = 0.05;

    // Use max_small_sigma
    for (int d=0; d<sim_dimensions; ++d) {
      dims[d] = static_cast<int>(domain_bounds.wd(d)/target_width);
      widths[d] = domain_bounds.wd(d)/dims[d];
      inverseW[d] = 1./widths[d];
    }
    minCutoff = 2*max_small_sigma + skin_depth;

    // Create the cells
    create_cells();

    // Construct the interaction handlers for the forces
    construct();

    // The domain has been initialized
    initialized = true;
  }

  void DomainTest::pre_integrate() {
    // Reset time points
    lastCheck  = -1.;
    lastUpdate = -1.;
    updateDelay = 1.0e-4;
  }

  void DomainTest::exchange_particles() {
    //! @todo Implement this.
  }

  void DomainTest::getAllWithin(int, RealType, vector<int>&) {
    //! @todo Implement this.
  }

  void DomainTest::construct() {
    // Domain base common tasks
    DomainBase::construct();
    // Clear out the cells
    clear_cells();
    // Fill the cells with particles
    fill_cells();

    RealType max_reasonable = sqr(0.9*bounds.wd(0));

    // Find potential neighbors
    RealType *sg = Base::simData->Sg();
    RealType **x = Base::simData->X();
    // Get the boundary conditions
    const BCFlag *bcs = gflow->getBCs();
    RealType dX[DIMENSIONS];
    for (const auto &c : cells) {
      for (auto p=c.particle_ids.begin(); p!=c.particle_ids.end(); ++p) {
        int id1 = *p;
        // All other particles in the same sector
        auto q = p;
        ++q;
        for (; q!=c.particle_ids.end(); ++q) {
          int id2 = *q;
          subtractVec(x[id1], x[id2], dX);
          RealType r2 = sqr(dX);
          if (r2 < sqr(sg[id1] + sg[id2] + skin_depth))
            pair_interaction(id1, id2);
        }

        // If sigma is <= than min_small_sigma, only look through cell stencil
        if (sg[id1]<=max_small_sigma) {
          // Seach through list of adjacent cells
          for (const auto &d : c.adjacent)
            for (const auto id2 : d->particle_ids) {
              // If the other particle is a large particle, it will take care of this interaction
              if (sg[id2]>max_small_sigma) continue;
              
              // Look for distance between particles
              subtractVec(x[id1], x[id2], dX);

              //getDisplacement(x[id1], x[id2], dX, bounds, bcs);
              RealType r2 = sqr(dX);

              if (r2 < sqr(sg[id1] + sg[id2] + skin_depth) || max_reasonable<r2)
                pair_interaction(id1, id2);

              /*
              if (r2 < sqr(sg[id1] + sg[id2] + skin_depth)) {
                pair_interaction(id1, id2);
              }
              */
            }
        }
        // If sigma is > min_small_sigma, we have to look through more cells
        else {
          throw false; // UNIMPLEMENTED
        }
      }
    }

    // Close all
    Base::forceMaster->close();
  }

  void DomainTest::setCellSize(RealType) {
    // @todo Implement this.
  }

  inline void DomainTest::create_cells() {
    // Initialize border record
    border_type_up = new int[sim_dimensions];
    border_type_down = new int[sim_dimensions];

    // --- Determine border type
    const BCFlag *bcs = Base::gflow->getBCs(); // Get the boundary condition flags
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

    // --- Create the cells
    // Get the total number of cells - The dims MUST be set first.
    const int size = getNumCells();
    cells = vector<CellTest>(size, CellTest());

    // Holder for tuple index
    int *tuple1 = new int[sim_dimensions], *tuple2 = new int[sim_dimensions];

    // --- Create a neighborhood stencil to help us find adjacent cells
    int sweep = ceil(minCutoff/min(widths, sim_dimensions));
    vector<int*> neighbor_indices;
    int *little_dims = new int[sim_dimensions], *center = new int[sim_dimensions];
    for (int d=0; d<sim_dimensions; ++d) little_dims[d] = 2*sweep+1;
    for (int d=0; d<sim_dimensions; ++d) center[d]      = sweep;
    // Create stencil
    for (int i=0; i<floor(0.5*pow(2*sweep+1, sim_dimensions)); ++i) {
      getAddressCM(i, little_dims, tuple1);
      subtractVec(tuple1, center, tuple2);
      int dr2 = sqr(tuple2);
      if (dr2<sqr(sweep+1)) {
        int *new_index = new int[sim_dimensions];
        copyVec(tuple2, new_index);
        neighbor_indices.push_back(new_index);
      }
    }

    // --- Assign cell types, adjacent cells
    for (int c=0; c<size; ++c) {
      linear_to_tuple(c, tuple1);
      for (auto n : neighbor_indices) {
        addVec(tuple1, n, tuple2);
        if (correct_index(tuple2)) {
          int linear;
          tuple_to_linear(linear, tuple2);
          cells[c].adjacent.push_back(&cells[linear]);
        }
      }
    }

    // Clean up 
    delete [] tuple1;
    delete [] tuple2;
    for (auto &n : neighbor_indices) delete [] n;
    delete [] little_dims;
    delete [] center;
  }

  inline void DomainTest::clear_cells() {
    for (auto &c : cells) c.particle_ids.clear();
  }

  inline void DomainTest::fill_cells() {
    RealType **X = simData->X();
    int number = Base::simData->number;
    // Bin all the particles
    for (int i=0; i<number; ++i) {
      int linear = get_cell_index(X[i]);
      // Stores the *local* id of the particle
      cells[linear].particle_ids.push_back(i);
    }
  }

  // Turns a linear cell index into a (DIMENSIONS)-dimensional index
  inline void DomainTest::linear_to_tuple(int linear, int *tuple) {
    // Same as the get address function in
    getAddressCM(linear, dims, tuple); // We need to use the column major form
  }

  // Turns a (DIMENSIONS)-dimensional index into a linear cell index
  // This is column major form.
  inline void DomainTest::tuple_to_linear(int &linear, const int *tuple) {
    // Product lambda
    auto product = [&] (int n) -> int {
      int total = 1;
      for (int i=n; i<sim_dimensions; ++i) total *= dims[i];
      return total;
    };

    linear = 0;
    for (int d=0; d<sim_dimensions; ++d)
      linear += tuple[d]*product(d+1);
  }

  inline bool DomainTest::correct_index(int *tuple) {
    bool good_index = true;
    const BCFlag *bcs = gflow->getBCs();
    for (int d=0; d<sim_dimensions; ++d) {
      if (tuple[d]>=dims[d]) {
        if (bcs[d]==BCFlag::WRAP) tuple[d] -= dims[d];
        else good_index = false;
      }
      if (tuple[d]<0) {
        if (bcs[d]==BCFlag::WRAP) tuple[d] += dims[d];
        else good_index = false;
      }
    }
    return good_index;
  }

  int DomainTest::get_cell_index(const RealType *x) {
    int index[DIMENSIONS], linear;
    for (int d=0; d<DIMENSIONS; ++d) {
      index[d] = static_cast<int>((x[d] - domain_bounds.min[d])*inverseW[d]);

      // Increment index if there are halo or ghost cells
      // if (border_type_down[d]) ++index[d];

      // Even when wrapping, rounding errors (I assume) can cause index to be too large.
      // When not wrapping, particles could be outside the sectorization grid
      // This also may be useful later for parallel things, if a particle should interact with
      // particles in this domain, but is so big that it is outside of the ghost cells.
      if (index[d]>=dims[d]) index[d] = dims[d]-1; 
      else if (index[d]<0)   index[d] = 0;
    }

    // Get the linear index
    tuple_to_linear(linear, index);

    // Return the index
    return linear;
  }

}