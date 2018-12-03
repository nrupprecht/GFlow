#include "domain_test.hpp"
// Other files
#include "../gflow.hpp"
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../interactionhandlers/verletlist.hpp"
#include "../base/forcemaster.hpp"

#include "../base/interaction.hpp"
#include "../base/interactionhandler.hpp"

namespace GFlowSimulation {

  // Dimension setting constructor
  CellTest::CellTest(int d) : sim_dimensions(d) {}

  // Copy constructor
  CellTest::CellTest(const CellTest& cell) : sim_dimensions(cell.sim_dimensions) {};

  CellTest CellTest::operator=(const CellTest& cell) {
    // Set data
    particle_ids = cell.particle_ids;
    adjacent = cell.adjacent;
    sim_dimensions = cell.sim_dimensions;
    // Return
    return *this;
  }

  int CellTest::size() const {
    return particle_ids.size();
  }

  // --------------

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

    // Find average sigma
    RealType sigma = 0, max_sigma = 0;
    for (int n=0; n<Base::simData->number; ++n) {
      RealType s = Base::simData->Sg(n);
      sigma += s;
      if (s>max_sigma) max_sigma = s;
    }
    sigma /= Base::simData->number;
    // Threshhold sigma is between average and maximum sigma
    RealType threshold = 0.5*(sigma + max_sigma), max_under = sigma;
    if (threshold!=sigma) {
      for (int n=0; n<Base::simData->number; ++n) {
        RealType s = Base::simData->Sg(n);
        if (s<threshold && max_under<s) max_under = s;
      }
    }
    max_small_sigma = 1.025*max_under;
    RealType target_width = (2*max_small_sigma+skin_depth);

    // Use max_small_sigma
    for (int d=0; d<sim_dimensions; ++d) {
      dims[d] = static_cast<int>(domain_bounds.wd(d)/target_width);
      if (dims[d]<=0) throw BadBounds(); // Make sure bin numbers are positive
      widths[d] = domain_bounds.wd(d)/dims[d];
      inverseW[d] = 1./widths[d];
    }
    minCutoff = 2*max_small_sigma + skin_depth;

    // Create the cells
    create_cells();

    // Construct the interaction handlers for the forces
    number = 0; // To make sure we do a full clear, fill when we construct.
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

  void DomainTest::removeOverlapping(RealType factor) {
    // Domain base common tasks
    DomainBase::construct();
    // Clear out the cells
    clear_cells();
    // Fill the cells with particles
    fill_cells();

    RealType max_reasonable = sqr(0.9*bounds.wd(0));

    // A tuple
    int *tuple1 = new int[sim_dimensions], *tuple2 = new int[sim_dimensions];
    int *cell_index = new int[sim_dimensions], *center = new int[sim_dimensions];

    // A list of ids of particles to remove
    std::set<int> removeList;

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
            subtractVec(x[id1], x[id2], dX, sim_dimensions);
            RealType r = magnitudeVec(dX, sim_dimensions);
            RealType overlap = sg[id1] + sg[id2] - r;
            if (overlap/min(sg[id1], sg[id2]) > factor)
              removeList.insert(sg[id1]>sg[id2] ? id2 : id1);
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
                removeList.insert(sg[id1]>sg[id2] ? id2 : id1);
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
                if (overlap/min(sg[id1], sg[id2]) > factor)
                  removeList.insert(sg[id1]>sg[id2] ? id2 : id1);
              }
            }
          }
        }
      }
    }
    // Remove all the particles that we need to 
    for (auto id : removeList)
      Base::simData->removeParticle(id);
    // Clean up
    delete [] tuple1;
    delete [] tuple2;
    delete [] cell_index;
    delete [] center;
    delete [] dX;
    delete [] search_dims;
  }

  void DomainTest::construct() {
    // Domain base common tasks
    DomainBase::construct();

    // Update particles in the cells
    update_cells();

    // Set a "maximum reasonable distance"
    RealType max_reasonable = sqr(0.9*bounds.wd(0));

    // A tuple
    int *tuple1 = new int[sim_dimensions], *tuple2 = new int[sim_dimensions];
    int *cell_index = new int[sim_dimensions], *center = new int[sim_dimensions];

    InteractionHandler *hs = nullptr;
    if (gflow->getInteractions().empty()) return;
    else hs = gflow->getInteractions()[0]->getInteractionHandler();

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
            if (r2 < sqr(sg[id1] + sg[id2] + skin_depth)) {
              //hs->addPair(id1, id2);
              pair_interaction(id1, id2);
            }
          }
          // Seach through list of adjacent cells
          for (const auto &d : c.adjacent)
            for (const auto id2 : d->particle_ids) {
              // If the other particle is a large particle, it will take care of this interaction
              if (sg[id2]>max_small_sigma) continue;
              // Look for distance between particles
              RealType r2 = getDistanceSqrNoWrap(x[id1], x[id2], sim_dimensions);
              if (r2 < sqr(sg[id1] + sg[id2] + skin_depth) || max_reasonable<r2) {
                //hs->addPair(id1, id2);
                pair_interaction(id1, id2);
              }
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

  void DomainTest::setCellSize(RealType) {
    // @todo Implement this.
  }

  inline void DomainTest::update_cells() {
    clear_cells();
    fill_cells();
    number = simData->number;
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
    cells = vector<CellTest>(size, CellTest(sim_dimensions));

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
    RealType **x = simData->X();
    int number = Base::simData->number;
    // Bin all the particles
    for (int i=0; i<number; ++i) {
      int linear = get_cell_index(x[i]);
      // Stores the *local* id of the particle
      cells[linear].particle_ids.push_back(i);
    }
  }

  // Turns a linear cell index into a (DIMENSIONS)-dimensional index
  inline void DomainTest::linear_to_tuple(int linear, int *tuple) {
    // Same as the get address function in
    getAddressCM(linear, dims, tuple, sim_dimensions); // We need to use the column major form
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
        if (bcs[d]==BCFlag::WRAP)
          tuple[d] -= dims[d];
        else good_index = false;
      }
      if (tuple[d]<0) {
        if (bcs[d]==BCFlag::WRAP) tuple[d] += dims[d];
        else good_index = false;
      }
    }
    return good_index;
  }

  inline void DomainTest::get_cell_index_tuple(const RealType *x, int *index) {
    for (int d=0; d<sim_dimensions; ++d) {
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
  }

  int DomainTest::get_cell_index(const RealType *x) {

    int linear = 0;

    // Product lambda
    auto product = [&] (int n) -> int {
      int total = 1;
      for (int i=n; i<sim_dimensions; ++i) total *= dims[i];
      return total;
    };

    for (int d=0; d<sim_dimensions; ++d) {
      RealType dth = static_cast<int>((x[d] - domain_bounds.min[d])*inverseW[d]);
      if (dth>=dims[d]) dth = dims[d]-1; 
      else if (dth<0)   dth = 0;
      linear += dth*product(d+1);
    }


    /*
    int index[DIMENSIONS], linear;

    // Get the tuple index of the cell
    get_cell_index_tuple(x, index);

    // Get the linear index
    tuple_to_linear(linear, index);
  */

    // Return the index
    return linear;
  }

}