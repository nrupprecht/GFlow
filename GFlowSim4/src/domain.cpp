#include "domain.hpp"
// Other files
#include "gflow.hpp"
#include "simdata.hpp"
#include "vectormath.hpp"
#include "verletlist.hpp"
#include "forcemaster.hpp"
#include "force.hpp"

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
    if (domain_bounds.wd(0) == 0) return;

    // --- Calculate cutoff

    // Largest and second largest radii
    RealType maxCutR(0), secCutR(0);
    // Look for largest and second largest radii
    for (int i=0; i<simData->number; ++i) {
      if (Base::simData->sg[i]>maxCutR) maxCutR = Base::simData->sg[i];
      else if (Base::simData->sg[i]>secCutR) secCutR = Base::simData->sg[i];
    }
    minCutoff = maxCutR + secCutR + skinDepth; // Cutoff radius

    // The actual cutoff is some multiple of the minimum cutoff
    cutoff = minCutoff*cutoffFactor;

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
    // @todo Implement this

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
      find_adjacent_cells(index, is_special, cells[ind].adjacent_cell_id);
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
    // @todo Implement this
  }

  void Domain::getAllWithin(int id, RealType radius, vector<int>& neighbors) {
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

      addVec(index, d_index, other_index);

      if (correct_index(index, is_special)) {
        tuple_to_linear(lin, other_index);
        for (const auto n_id : cells[lin].id_list) {
          getDisplacement(Base::simData->x[id], Base::simData->x[n_id], displacement, bounds, bcs);
          // Check that the particle is close enough
          if (sqr(displacement)<sqr(radius))
            neighbors.push_back(n_id);
        }
      }
    }
  }

  void Domain::addToCell(int id) {
    int index = getCellIndex(simData->X(id));
    Cell &cell = cells[index];
    cell.add(id);
    // @todo Create halo particles if necessary
  }

  int Domain::getCellIndex(RealType *x) {
    int index[DIMENSIONS];
    for (int d=0; d<DIMENSIONS; ++d)
      index[d] = static_cast<int>((x[d] - extended_domain_bounds.min[d])*inverseW[d]);
    // Return the linear index
    int linear;
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
    // Increment counter
    ++number_of_remakes;
    // Record where the particles were
    fillXVL();
    // Clear old verlet lists
    Base::forceMaster->clearVerletLists();

    // Fill the linked cells
    fill_cells();

    // Get the boundary conditions
    const BCFlag *bcs = gflow->getBCs();

    // Pick the appropriate displacement function
    auto Displacement_Halo = [] (const RealType x1[DIMENSIONS], const RealType x2[DIMENSIONS], 
      RealType d[DIMENSIONS], const Bounds bounds, const BCFlag *bcs) 
    {
      subtractVec(x1, x2, d);
    };
    auto Displacement_No_Halo = [] (const RealType x1[DIMENSIONS], const RealType x2[DIMENSIONS], 
      RealType d[DIMENSIONS], const Bounds bounds, const BCFlag *bcs) 
    {
      getDisplacement(x1, x2, d, bounds, bcs);
    };
    bool no_wrap = true;
    for (int d=0; d<DIMENSIONS; ++d) 
      if (bcs[d]==BCFlag::WRAP) no_wrap = false;
    auto Displacement = (use_halo_cells || no_wrap) ? Displacement_Halo : Displacement_No_Halo;

    // --- Calculate verlet lists
    RealType displacement[DIMENSIONS];
    for (const auto &cell1 : cells) {
      // Loop through every particle in the cell
      for (auto i1 = cell1.id_list.cbegin(); i1!=cell1.id_list.cend(); ++i1) {
        const int id1 = *i1;
        RealType sigma = Base::simData->sg[id1];
        // Verlet lists for particles in the same cell
        for (auto i2 = i1+1; i2!=cell1.id_list.cend(); ++i2) {
          // Get the id of the second particle
          int id2 = *i2;
          // Get the displacement
          Displacement(Base::simData->x[id1], Base::simData->x[id2], displacement, bounds, bcs);
          // Check the displacement
          if (sqr(displacement) < sqr(sigma + Base::simData->sg[id2] + skinDepth))
            pair_interaction(id1, id2);
        }

        // Loop through adjacent cells
        for (auto cell_id : cell1.adjacent_cell_id) {
          const Cell &cell2 = cells[cell_id];
          // Loop through the particles in the adjacent cells
          for (auto id2 : cell2.id_list) {
            // Get the displacement
            Displacement(Base::simData->x[id1], Base::simData->x[id2], displacement, bounds, bcs);
            // Check the displacement
            if (sqr(displacement) < sqr(sigma + Base::simData->sg[id2] + skinDepth))
              pair_interaction(id1, id2);
          }
        }
      }
    }
  }

  // Turns a linear cell index into a (DIMENSIONS)-dimensional index
  void Domain::linear_to_tuple(int linear, int *tuple) {
    // Same as the get address function in
    getAddress(linear, dims, tuple);
  }

  // Turns a (DIMENSIONS)-dimensional index into a linear cell index
  void Domain::tuple_to_linear(int &linear, const int *tuple) {
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

  void Domain::find_adjacent_cells(int index[DIMENSIONS], bool is_special, vector<int>& neighbors) {
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
      neighbors.push_back(linear);
    }
  }

  inline void Domain::fill_cells() {
    // Clear cells of old data
    for (int i=0; i<cells.size(); ++i) 
      cells[i].clear();

    // @todo Clear ghost index mapping in simData

    // --- Place particles in their cells
    int index[DIMENSIONS], linear;
    for (int n=0; n<Base::simData->number; ++n) {
      // Calculate the (DIMENSION)-tuple cell coordinate
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
      cells[linear].add(n);

      // @todo Add cell to halo cells if neccessary

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
