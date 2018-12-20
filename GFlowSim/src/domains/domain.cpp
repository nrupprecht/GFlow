#include "domain.hpp"
// Other files
#include "../gflow.hpp"
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../interactionhandlers/verletlist.hpp"
#include "../base/forcemaster.hpp"

namespace GFlowSimulation {

  Domain::Domain(GFlow *gflow) : DomainBase(gflow), use_halo_cells(false), extended_domain_bounds(Bounds(2)) {
    // --- Find what the domain bounds are
    // For now, we are only running on one core, this is the only domain
    domain_bounds = gflow->getBounds();
    extended_domain_bounds = Bounds(sim_dimensions);

    domains_up = new int[sim_dimensions];
    domains_down = new int[sim_dimensions];
    domain_dims = new int[sim_dimensions];
    domain_index = new int[sim_dimensions];

    halo_up = new bool[sim_dimensions];
    halo_down = new bool[sim_dimensions];
    ghost_up = new bool[sim_dimensions];
    ghost_down = new bool[sim_dimensions];

    // Set some default values
    setVec(domains_up, -1, sim_dimensions);
    setVec(domains_down, -1, sim_dimensions);
    setVec(domain_dims, 0, sim_dimensions);
    setVec(domain_index, 0, sim_dimensions);
    setVec(halo_up, false, sim_dimensions);
    setVec(halo_down, false, sim_dimensions);
    setVec(ghost_up, false, sim_dimensions);
    setVec(ghost_down, false, sim_dimensions);

    #if USE_MPI==1
    #if _CLANG_ == 1
    // MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    #else
    // MPI::Init(argc, argv);
    rank = MPI::COMM_WORLD.Get_rank();
    numProc = MPI::COMM_WORLD.Get_size();
    #endif
    // Initialize flag
    parallel_init = false;
    #endif
  };

  Domain::~Domain() {
    delete [] domains_up;
    delete [] domains_down;
    delete [] domain_dims;
    delete [] domain_index;
    delete [] halo_up;
    delete [] halo_down;
    delete [] ghost_up;
    delete [] ghost_down;
  }

  void Domain::initialize() {
    // First call Base's initialize function
    Base::initialize();

    // Get the simulation bounds from gflow
    bounds = Base::gflow->getBounds();
    domain_bounds = bounds;

    // If bounds are unset, then don't make sectors
    if (domain_bounds.vol()<=0) return;

    // We cannot initialize if simdata is null
    if (simData==nullptr) return;

    // --- Calculate cutoff
    RealType *sg = simData->Sg();
    int *type = simData->Type();
    int size = simData->size();

    // We may have already set the cutoff, and don't want it to be overwritten
    if (cutoff <= minCutoff) {
      // Largest radius
      RealType maxSigma(0);
      // Look for largest and second largest radii
      for (int i=0; i<size; ++i)
        if (type[i]>-1 && sg[i]>maxSigma) maxSigma = sg[i];
      max_small_sigma = maxSigma;
      minCutoff = 2*max_small_sigma + skin_depth; // Cutoff radius

      // The actual cutoff distance is some multiple of the minimum cutoff
      cutoff = minCutoff*cutoffFactor;
      // Create the cells
      create_cells();
    }
    // Make verlet lists
    construct();
  }

  void Domain::pre_integrate() {
    // Reset time points
    lastCheck  = -1.;
    lastUpdate = -1.;
    updateDelay = 1.0e-4;
  }

  void Domain::exchange_particles() {
    //! @todo Implement this
  }

  void Domain::getAllWithin(int id, RealType radius, vector<int>& neighbors) {
    // Clear the vector
    neighbors.clear();

    // Find where particle [id] is
    RealType *X = Base::simData->X(id);
    RealType *displacement = new RealType[sim_dimensions];
    int linear = getCellIndex(X), lin;
    int *index = new int[sim_dimensions], *d_index = new int[sim_dimensions], *other_index = new int[sim_dimensions];
    linear_to_tuple(linear, index);
    const BCFlag *bcs = Base::gflow->getBCs();

    bool is_special = ( cells[linear].cellType==CellType::Halo || cells[linear].cellType==CellType::Ghost );
    // --- DO THIS FOR NOW
    // Go through all adjacent cells
    for (int c=0; c<pow(3, sim_dimensions); ++c) {
      // Convert to a base 3 number (-1, 0, 1)
      int c0=c;
      for (uint d=0; d<sim_dimensions; ++d) {
        d_index[d] = (c0%3) - 1;
        c0 /= 3;
      }
      // Add the index displacement to the index to get a (potentially) adjacent index
      addVec(index, d_index, other_index, sim_dimensions);
      // Check if this cell really should be adjacent and if its a real cell
      if (correct_index(other_index, is_special)) {
        tuple_to_linear(lin, other_index);
        for (const auto n_id : cells[lin].id_list) {
          if (id==n_id) continue;
          getDisplacement(Base::simData->X(id), Base::simData->X(n_id), displacement, bounds, bcs, sim_dimensions);
          // Check that the particle is close enough
          if (sqr(displacement, sim_dimensions)<sqr(radius))
            neighbors.push_back(n_id);
        }
      }
    }
    // Clean up
    delete [] displacement;
    delete [] index;
    delete [] d_index;
    delete [] other_index;
  }

  void Domain::setSkinDepth(RealType s) {
    DomainBase::setSkinDepth(s);
    // Have to recalculate everything since we messed with the skin depth
    initialize();
  }

  void Domain::setCellSize(RealType cell_size) {
    // We don't set a cutoff less than the minimum
    if (cell_size<minCutoff) return;
    // Set the target cutoff
    cutoff = cell_size;
    // Recreate cells
    create_cells();
  }

  void Domain::addToCell(int id) {
    int index = getCellIndex(simData->X(id));
    Cell &cell = cells[index];
    cell.add(id);
  }

  int Domain::getCellIndex(RealType *x) {

    int linear = 0;
    // Product lambda
    auto product = [&] (int n) -> int {
      int total = 1;
      for (int i=n; i<sim_dimensions; ++i) total *= dims[i];
      return total;
    };
    // Calculate
    for (int d=0; d<sim_dimensions; ++d) {
      RealType dth = static_cast<int>((x[d] - domain_bounds.min[d])*inverseW[d]);
      if (dth>=dims[d]) dth = dims[d]-1; 
      else if (dth<0)   dth = 0;
      linear += dth*product(d+1);
    }
    return linear;
  }

  void Domain::construct() {
    // Setup and common tasks
    DomainBase::construct();

    // Fill the linked cells - this is where particles outside the domain, or in the 
    // halo region are detected.
    fill_cells();

    // Get the boundary conditions
    const BCFlag *bcs = gflow->getBCs();

    // A simple displacement function that we might get to use
    auto Displacement_Halo = [&] (const RealType *x1, const RealType *x2, 
      RealType *d, const Bounds &bounds, const BCFlag *bcs) 
    {
      for (int i=0; i<sim_dimensions; ++i)
        d[i] = x1[i] - x2[i];
    };

    // Decide which if we can use a simple displacement function
    bool no_wrap = true;
    for (int d=0; d<sim_dimensions; ++d) 
      if (bcs[d]==BCFlag::WRAP) no_wrap = false;
    // Choose displacement function
    //void (*Displacement) (const RealType*, const RealType*, RealType*, const Bounds&, const BCFlag*); 
    //Displacement = (use_halo_cells || no_wrap) ? Displacement_Halo : getDisplacement;

    RealType *sg = Base::simData->Sg();

    // --- Calculate verlet lists
    RealType *displacement = new RealType[sim_dimensions];
    for (const auto &cell1 : cells) {
      // Loop through every particle in the cell
      for (auto i1 = cell1.id_list.cbegin(); i1!=cell1.id_list.cend(); ++i1) {
        const int id1 = *i1;
        RealType sigma = sg[id1];
        // --- Verlet lists for particles in the same cell
        for (auto i2=i1+1; i2!=cell1.id_list.cend(); ++i2) {
          // Get the id of the second particle
          int id2 = *i2;
          // Get the displacement - we never need to wrap, since its the same cell (and particle positions don't wrap 
          // until remake_verlet is called again).
          subtractVec(Base::simData->X(id1), Base::simData->X(id2), displacement, sim_dimensions);
          // Check the displacement
          if (sqr(displacement, sim_dimensions) < sqr(sigma + sg[id2] + skin_depth))
            pair_interaction(id1, id2);
        }

        // Loop through adjacent cells
        for (auto cell_ptr : cell1.adjacent_cells) {
          // Loop through the particles in the adjacent cells
          for (auto id2 : cell_ptr->id_list) {
            // Get the displacement
            getDisplacement(Base::simData->X(id1), Base::simData->X(id2), displacement, bounds, bcs, sim_dimensions); // Faster to just do this
            // Check the displacement
            if (sqr(displacement, sim_dimensions) < sqr(sigma + sg[id2] + skin_depth))
              pair_interaction(id1, id2);
          }
        }
      }
    }
    // Close all
    Base::forceMaster->close();
    // Clean up 
    delete [] displacement;
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
    auto product = [&] (int n) -> int {
      int total = 1;
      for (int i=n; i<sim_dimensions; ++i) total *= dims[i];
      return total;
    };

    linear = 0;
    for (int d=0; d<sim_dimensions; ++d)
      linear += tuple[d]*product(d+1);
  }

  inline void Domain::find_adjacent_cells(int *index, bool is_special, vector<Cell*>& neighbors) {
    // Helper arrays
    int *d_index = new int[sim_dimensions], *other_index = new int[sim_dimensions], linear;

    // --- Look in the first half (rounded down) sectors only
    for (uint c=0; c<floor(pow(3,sim_dimensions)/2); ++c) {
      // Convert to a base 3 number (-1, 0, 1)
      int c0=c;
      for (uint d=0; d<sim_dimensions; ++d) {
        d_index[d] = (c0%3) - 1;
        c0 /= 3;
      }
      // Look at that neighboring sector
      addVec(index, d_index, other_index, sim_dimensions);

      // --- Check that the other cell actually exists

      bool good = true;
      if (use_halo_cells) {
        // We are using harmonic cells for wrapping, so if an index is negative or to large, we are done
        for (int d=0; d<sim_dimensions && good; ++d) {
          if (other_index[d]<0) {
            good = false;
            break;
          }
          if (other_index[d]>=dims[d]) {
            good = false;
            break;
          }
        }
      }
      else { // There are no halo cells, and neighbors could be cells that we have to wrap to get to 
        const BCFlag *bcs = gflow->getBCs();
        for (int d=0; d<sim_dimensions && good; ++d) {
          // Ghost cells will not have wrapping neighbors, hence the check for is_special
          if (other_index[d]<0 && !is_special) {
            if (bcs[d]==BCFlag::WRAP && dims[d]>1) other_index[d] += dims[d];
            else {
              good = false;
              break; 
            }
          }
          // Ghost cells will not have wrapping neighbors, hence the check for is_special
          if (other_index[d]>=dims[d] && !is_special) {
            if (bcs[d]==BCFlag::WRAP && dims[d]>1) other_index[d] -= dims[d];
            else {
              good = false;
              break;
            }
          }
        }
      }
      // Do not evaluate if we shouldn't
      if (!good) continue;

      // The cell at other_index
      tuple_to_linear(linear, other_index);
      // A ghost or halo cell can only be neighbors with a central cell
      if (is_special && cells[linear].cellType!=CellType::Central) continue;
      // Include the cell
      neighbors.push_back(&cells[linear]);
    }
    // Clean up
    delete [] d_index;
    delete [] other_index;
  }

  inline void Domain::create_cells() {
    // Correct estimate so sectors are no smaller than our estimation
    for (int d=0; d<sim_dimensions; ++d) {
      widths[d] = cutoff;
      dims[d] = static_cast<int>(max(RealType(1.), bounds.wd(d)/widths[d]));
      // To solve (temporarily?) the "two sector" problem - see comments in "makeVerletLists"
      if (dims[d]==2) dims[d] = 1; 
      // Set widths and inverse widths
      widths[d] = bounds.wd(d)/dims[d];
      inverseW[d] = 1./widths[d];
    }

    // --- Communicate with other domains so we have common cell dimensions
    // @todo Implement this. ---> This actually probably isn't necessary.

    // --- Assign neighbor cells to every cell    
    const BCFlag *bcs = Base::gflow->getBCs(); // Get the boundary condition flags
    // Figure out what types of "extra" cells we'll need (other than normal central cells)
    for (int d=0; d<sim_dimensions; ++d) {
      // Upper cells
      if (domains_up[d]==-1) {
        if (bcs[d]==BCFlag::WRAP) {
          halo_up[d]    = use_halo_cells;
          halo_down[d]  = use_halo_cells;
          ghost_up[d]   = false;
          ghost_down[d] = false;
        }
        else {
          halo_up[d]    = false;
          halo_down[d]  = false;
          ghost_up[d]   = false;
          ghost_down[d] = false;
        }
      }
      else {
        halo_up[d]    = false;
        halo_down[d]  = false;
        ghost_up[d]   = true;
        ghost_down[d] = true;
      }

      // Lower cells
      if (domains_down[d]==-1) {
        if (bcs[d]==BCFlag::WRAP) {
          halo_up[d]    = use_halo_cells;
          halo_down[d]  = use_halo_cells;
          ghost_up[d]   = false;
          ghost_down[d] = false;
        }
        else {
          halo_up[d]    = false;
          halo_down[d]  = false;
          ghost_up[d]   = false;
          ghost_down[d] = false;
        }
      }
      else {
        halo_up[d]    = false;
        halo_down[d]  = false;
        ghost_up[d]   = true;
        ghost_down[d] = true;
      }
    }

    // Correct dims to reflect the new sectors, adjust extended_bounds
    extended_domain_bounds = domain_bounds;
    for (int d=0; d<sim_dimensions; ++d) {
      if (halo_down[d] || ghost_down[d]) {
        ++dims[d];
        extended_domain_bounds.min[d] -= widths[d];
      }
      if (halo_up[d] || ghost_up[d]) {
        ++dims[d];
        extended_domain_bounds.max[d] += widths[d];
      }
    }

    // --- Create the cells
    const int size = getNumCells();
    cells = vector<Cell>(size, Cell());

    // Assign cell types
    int *index = new int[sim_dimensions];
    for (int ind=0; ind<size; ++ind) {
      linear_to_tuple(ind, index);
      // Default is that this is a central cell
      cells[ind].cellType = CellType::Central;
      for (int d=0; d<sim_dimensions; ++d) {
        // Check if this cell is harmonic or a ghost
        if (index[d]==0) {
          if (!ghost_down[d] && halo_down[d])
            cells[ind].cellType = CellType::Halo;
          else if (ghost_down[d])
            cells[ind].cellType = CellType::Ghost;
        }
        if (index[d]==dims[d]-1) {
          if (!ghost_up[d] && halo_up[d]) 
            cells[ind].cellType = CellType::Halo;
          else if (ghost_up[d])
            cells[ind].cellType = CellType::Ghost;
        }
      }
      // Do this after this cell's type has been assigned
      if (cells[ind].cellType==CellType::Central) 
        for (int d=0; d<sim_dimensions; ++d) {
          // Check if this cell is at the edge of the regular domain and has harmonic image cells
          if ((index[d]==1 && halo_down[d]) || (index[d]==dims[d]-2 && halo_up[d])) 
            cells[ind].is_boundary_cell = true;
        }

    }

    // Find neighbors - must be done after cell types are assigned
    for (int ind=0; ind<size; ++ind) {
      linear_to_tuple(ind, index);
      bool is_special = ( cells[ind].cellType==CellType::Halo || cells[ind].cellType==CellType::Ghost );
      // Find cell's neighbors
      find_adjacent_cells(index, is_special, cells[ind].adjacent_cells);
    }
    // Clean up 
    delete [] index;
  }

  inline void Domain::fill_cells() {
    // Clear cells of old data
    for (int i=0; i<cells.size(); ++i) 
      cells[i].clear();

    // We should have just done a particle removal, so we can use number, not size (since all arrays are compressed
    int number = Base::simData->number(), linear;
    // Place particles in their cells
    for (int n=0; n<number; ++n) {
      // Calculate the (DIMENSION)-tuple cell coordinate
      int linear = getCellIndex(Base::simData->X(n));
      cells[linear].add(n);
    }
  }

  inline bool Domain::correct_index(int *index, bool is_special) {
    bool good = true;
    if (use_halo_cells) {
      // We are using harmonic cells for wrapping, so if an index is negative or to large, we are done
      for (int d=0; d<sim_dimensions && good; ++d) {
        if (index[d]<0) {
          good = false;
          break;
        }
        if (index[d]>=dims[d]) {
          good = false;
          break;
        }
      }
    }
    else { // There are no halo cells, and neighbors could be cells that we have to wrap to get to 
      const BCFlag *bcs = gflow->getBCs();
      for (int d=0; d<sim_dimensions && good; ++d) {
        // Ghost cells will not have wrapping neighbors, hence the check for is_special
        if (index[d]<0 && !is_special) {
          if (bcs[d]==BCFlag::WRAP && dims[d]>1) index[d] += dims[d];
          else {
            good = false;
            break; 
          }
        }
        // Ghost cells will not have wrapping neighbors, hence the check for is_special
        if (index[d]>=dims[d] && !is_special) {
          if (bcs[d]==BCFlag::WRAP && dims[d]>1) index[d] -= dims[d];
          else {
            good = false;
            break;
          }
        }
      }
    }
    return good;
  }

}
