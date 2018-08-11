#include "domain.hpp"
// Other files
#include "gflow.hpp"
#include "simdata.hpp"
#include "vectormath.hpp"
#include "verletlist.hpp"
#include "forcemaster.hpp"

namespace GFlowSimulation {

  Domain::Domain(GFlow *gflow) : DomainBase(gflow), use_halo_cells(false) {
    // --- Find what the domain bounds are
    // For now, we are only running on one core, this is the only domain
    domain_bounds = gflow->getBounds();

    // Set some default values
    setVec(domains_up, -1);
    setVec(domains_down, -1);
    setVec(domain_dims, 0);
    setVec(domain_index, 0);
    setVec(halo_up, false);
    setVec(halo_down, false);
    setVec(ghost_up, false);
    setVec(ghost_down, false);
  };

  void Domain::initialize() {
    // First call Base's initialize function
    Base::initialize();

    // Get the simulation bounds from gflow
    bounds = Base::gflow->getBounds();

    // Calculate data about the domain decomposition, and this domain's place in it
    parallel_assignments();

    // If bounds are unset, then don't make sectors
    if (domain_bounds.wd(0)==0) return;

    // --- Calculate cutoff
    RealType *sg = simData->sg;
    int number = simData->number;

    // We may have already set the cutoff, and don't want it to be overwritten
    if (cutoff <= minCutoff) {
      // Largest radius
      RealType maxSigma(0);
      // Look for largest and second largest radii
      for (int i=0; i<number; ++i)
        if (sg[i]>maxSigma) maxSigma = sg[i];
      max_small_sigma = maxSigma;
      minCutoff = 2*max_small_sigma + skin_depth; // Cutoff radius

      // The actual cutoff distance is some multiple of the minimum cutoff
      cutoff = minCutoff*cutoffFactor;

      // Create the cells
      create_cells();
    }

    // Make verlet lists
    remake_verlet();
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
    RealType *X = Base::simData->x[id];
    RealType displacement[DIMENSIONS];
    int linear = getCellIndex(X), lin;
    int index[DIMENSIONS], d_index[DIMENSIONS], other_index[DIMENSIONS];
    linear_to_tuple(linear, index);
    const BCFlag *bcs = Base::gflow->getBCs();

    bool is_special = ( cells[linear].cellType==CellType::Halo || cells[linear].cellType==CellType::Ghost );
    // --- DO THIS FOR NOW
    // Go through all adjacent cells
    for (int c=0; c<pow(3, DIMENSIONS); ++c) {
      // Convert to a base 3 number (-1, 0, 1)
      int c0=c;
      for (uint d=0; d<DIMENSIONS; ++d) {
        d_index[d] = (c0%3) - 1;
        c0 /= 3;
      }
      // Add the index displacement to the index to get a (potentially) adjacent index
      addVec(index, d_index, other_index);
      // Check if this cell really should be adjacent and if its a real cell
      if (correct_index(other_index, is_special)) {
        tuple_to_linear(lin, other_index);
        for (const auto n_id : cells[lin].id_list) {
          if (id==n_id) continue;
          getDisplacement(Base::simData->x[id], Base::simData->x[n_id], displacement, bounds, bcs);
          // Check that the particle is close enough
          if (sqr(displacement)<sqr(radius))
            neighbors.push_back(n_id);
        }
      }
    }
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
    int index[DIMENSIONS], linear;
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) {
      index[d] = static_cast<int>((x[d] - extended_domain_bounds.min[d])*inverseW[d]);
      // Even when wrapping, rounding errors (I assume) can cause index to be too large.
      // When not wrapping, particles could be outside the sectorization grid
      // This also may be useful later for parallel things, if a particle should interact with
      // particles in this domain, but is so big that it is outside of the ghost cells.
      if (index[d]>=dims[d]) index[d] = dims[d]-1; 
      else if (index[d]<0)   index[d] = 0;
    }
    // Add particle to the appropriate sector
    tuple_to_linear(linear, index);
    return linear;
  }

  void Domain::parallel_assignments() {
    // Right now, this function doesn't actually work with MPI or calculate anything about
    // parallel decomposition. It only works with single processor simulations.

    // Calculate how the divison of the simulation into domains
    // @todo Actually calculate this
    setVec(domain_dims, 1); // For now ...

    // Calculate which domain we are
    // @todo Actually calculate this
    zeroVec(domain_index); // For now ...

    // Calculate which domains are adjacent to us
    // @todo Actually calculate these
    setVec(domains_up, -1);   // For now ...
    setVec(domains_down, -1); // For now ...

    // Calculate the bounds of this domain
    domain_bounds = bounds; // For now ...
  }

  void Domain::remake_verlet() {
    // Setup and common tasks
    DomainBase::remake_verlet();

    // Fill the linked cells
    fill_cells();

    // Get the boundary conditions
    const BCFlag *bcs = gflow->getBCs();

    // A simple displacement function that we might get to use
    auto Displacement_Halo = [] (const RealType x1[DIMENSIONS], const RealType x2[DIMENSIONS], 
      RealType d[DIMENSIONS], const Bounds bounds, const BCFlag *bcs) 
    {
      for (int i=0; i<DIMENSIONS; ++i)
        d[i] = x1[i] - x2[i];
      // subtractVec(x1, x2, d);
    };

    // Decide which if we can use a simple displacement function
    bool no_wrap = true;
    for (int d=0; d<DIMENSIONS; ++d) 
      if (bcs[d]==BCFlag::WRAP) no_wrap = false;
    // Choose displacement function
    void (*Displacement) (const RealType*, const RealType*, RealType*, const Bounds, const BCFlag*); 
    Displacement = (use_halo_cells || no_wrap) ? Displacement_Halo : getDisplacement;

    // --- Calculate verlet lists
    RealType displacement[DIMENSIONS];
    for (const auto &cell1 : cells) {
      // Loop through every particle in the cell
      for (auto i1 = cell1.id_list.cbegin(); i1!=cell1.id_list.cend(); ++i1) {
        const int id1 = *i1;
        RealType sigma = Base::simData->sg[id1];
        // --- Verlet lists for particles in the same cell
        for (auto i2=i1+1; i2!=cell1.id_list.cend(); ++i2) {
          // Get the id of the second particle
          int id2 = *i2;
          // Get the displacement - we never need to wrap, since its the same cell (and particle positions dont wrap 
          // until remake_verlet is called again).
          subtractVec(Base::simData->x[id1], Base::simData->x[id2], displacement);
          // Check the displacement
          if (sqr(displacement) < sqr(sigma + Base::simData->sg[id2] + skin_depth))
            pair_interaction(id1, id2);
        }

        // Loop through adjacent cells
        for (auto cell_ptr : cell1.adjacent_cells) {
          // Loop through the particles in the adjacent cells
          for (auto id2 : cell_ptr->id_list) {
            // Get the displacement
            // Displacement(Base::simData->x[id1], Base::simData->x[id2], displacement, bounds, bcs);
            getDisplacement(Base::simData->x[id1], Base::simData->x[id2], displacement, bounds, bcs); // Faster to just do this
            // Check the displacement
            if (sqr(displacement) < sqr(sigma + Base::simData->sg[id2] + skin_depth))
              pair_interaction(id1, id2);
          }
        }
      }
    }
  }

  // Turns a linear cell index into a (DIMENSIONS)-dimensional index
  inline void Domain::linear_to_tuple(int linear, int *tuple) {
    // Same as the get address function in
    getAddress(linear, dims, tuple);
  }

  // Turns a (DIMENSIONS)-dimensional index into a linear cell index
  inline void Domain::tuple_to_linear(int &linear, const int *tuple) {
    // Product lambda
    auto product = [&] (int n) -> int {
      int total = 1;
      for (int i=n; i<DIMENSIONS; ++i) total *= dims[i];
      return total;
    };

    linear = 0;
    for (int d=0; d<DIMENSIONS; ++d)
      linear += tuple[d]*product(d+1);
  }

  inline void Domain::find_adjacent_cells(int index[DIMENSIONS], bool is_special, vector<Cell*>& neighbors) {
    // Helper arrays
    int d_index[DIMENSIONS], other_index[DIMENSIONS], linear;

    // --- Look in the first half (rounded down) sectors only
    for (uint c=0; c<floor(pow(3,DIMENSIONS)/2); ++c) {
      // Convert to a base 3 number (-1, 0, 1)
      int c0=c;
      for (uint d=0; d<DIMENSIONS; ++d) {
        d_index[d] = (c0%3) - 1;
        c0 /= 3;
      }
      // Look at that neighboring sector
      addVec(index, d_index, other_index);

      // --- Check that the other cell actually exists

      bool good = true;
      if (use_halo_cells) {
        // We are using harmonic cells for wrapping, so if an index is negative or to large, we are done
        for (int d=0; d<DIMENSIONS && good; ++d) {
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
        for (int d=0; d<DIMENSIONS && good; ++d) {
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
  }

  inline void Domain::create_cells() {
    // Correct estimate so sectors are no smaller than our estimation
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) {
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
    for (int d=0; d<DIMENSIONS; ++d) {
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
    for (int d=0; d<DIMENSIONS; ++d) {
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
    int index[DIMENSIONS];
    for (int ind=0; ind<size; ++ind) {
      linear_to_tuple(ind, index);
      // Default is that this is a central cell
      cells[ind].cellType = CellType::Central;
      for (int d=0; d<DIMENSIONS; ++d) {
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

        /// @todo Maybe do something about checking if cells are adjacent to 
      }
      // Do this after this cell's type has been assigned
      if (cells[ind].cellType==CellType::Central) 
        for (int d=0; d<DIMENSIONS; ++d) {
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
  }

  inline void Domain::fill_cells() {
    // Clear cells of old data
    for (int i=0; i<cells.size(); ++i) 
      cells[i].clear();

    // @todo Clear ghost index mapping in simData

    // --- Place particles in their cells
    int linear;
    for (int n=0; n<Base::simData->number; ++n) {
      // Calculate the (DIMENSION)-tuple cell coordinate

      int linear = getCellIndex(Base::simData->x[n]);

      /*
      #if _INTEL_ == 1
      #pragma unroll(DIMENSIONS)
      #endif 
      for (int d=0; d<DIMENSIONS; ++d) {
        index[d] = static_cast<int>((Base::simData->X(n, d) - extended_domain_bounds.min[d])*inverseW[d]);
        // Even when wrapping, rounding errors (I assume) can cause index to be too large.
        // When not wrapping, particles could be outside the sectorization grid
        // This also may be useful later for parallel things, if a particle should interact with
        // particles in this domain, but is so big that it is outside of the ghost cells.
        if (index[d]>=dims[d]) index[d] = dims[d]-1; 
        else if (index[d]<0)   index[d] = 0;
      }

      // Add particle to the appropriate sector
      tuple_to_linear(linear, index);
      */
      cells[linear].add(n);

      /*
      // Add cell to halo cells if neccessary
      if (cells[linear].is_boundary_cell) {
        // In what dimensions do we need to make halo particles
        vector<int> changes;
        vector<RealType> displacement;
        for (int d=0; d<DIMENSIONS; ++d) {
          if (halo_down[d]) { // If halo_down[d], then halo_up[d]
            // Upwards displacent
            if (index[d]==1) {
              displacement.push_back(domain_bounds.wd(d));
              changes.push_back(d);
            }
            else if (index[d]==dims[d]-2) {
              displacement.push_back(-domain_bounds.wd(d));
              changes.push_back(d);
            }
          }
        }
        // Create halo particles
        for (int i=0; i<pow(2, changes.size())-1; ++i) {
          RealType Xd[DIMENSIONS];
          copyVec(Base::simData->x[n], Xd);

        }
        
        // Done making halo particles
      }
      */
    }
  }

  inline bool Domain::correct_index(int index[DIMENSIONS], bool is_special) {
    bool good = true;
    if (use_halo_cells) {
      // We are using harmonic cells for wrapping, so if an index is negative or to large, we are done
      for (int d=0; d<DIMENSIONS && good; ++d) {
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
      for (int d=0; d<DIMENSIONS && good; ++d) {
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
