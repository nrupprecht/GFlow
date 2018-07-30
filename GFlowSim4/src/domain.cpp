#include "domain.hpp"
// Other files
#include "gflow.hpp"
#include "simdata.hpp"
#include "vectormath.hpp"
#include "verletlist.hpp"
#include "forcemaster.hpp"
#include "force.hpp"

namespace GFlowSimulation {

  Domain::Domain(GFlow *gflow) : DomainBase(gflow) {
    // --- Find what the domain bounds are
    // For now, we are only running on one core, this is the only domain
    domain_bounds = gflow->getBounds();

    // --- Find neighboring domains
    // For now, we are only running on one core
    for (int d=0; d<DIMENSIONS; ++d) {
      domains_up[d] = -1;
      domains_down[d] = -1;
    }

  };

  void Domain::initialize() {
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
    
    // --- Assign neighbor cells to every cell

    // Get the boundary condition flags
    const BCFlag *bcs = Base::gflow->getBCs();

    // Figure out what types of "extra" cells we'll need (other than normal central cells)
    for (int d=0; d<DIMENSIONS; ++d) {
      // Upper cells
      if (domains_up[d]==-1) {
        if (bcs[d]==BCFlag::WRAP) {
          harmonic_up[d]   = true;
          harmonic_down[d] = true;
          ghost_up[d]      = false;
          ghost_down[d]    = false;
        }
        else {
          harmonic_up[d]   = false;
          harmonic_down[d] = false;
          ghost_up[d]      = false;
          ghost_down[d]    = false;
        }
      }
      else {
        harmonic_up[d]   = false;
        harmonic_down[d] = false;
        ghost_up[d]      = true;
        ghost_down[d]    = true;
      }

      // Lower cells
      if (domains_down[d]==-1) {
        if (bcs[d]==BCFlag::WRAP) {
          harmonic_up[d]   = true;
          harmonic_down[d] = true;
          ghost_up[d]      = false;
          ghost_down[d]    = false;
        }
        else {
          harmonic_up[d]   = false;
          harmonic_down[d] = false;
          ghost_up[d]      = false;
          ghost_down[d]    = false;
        }
      }
      else {
        harmonic_up[d]   = false;
        harmonic_down[d] = false;
        ghost_up[d]      = true;
        ghost_down[d]    = true;
      }
    }

    // Correct dims to reflect the new sectors, adjust extended_bounds
    extended_domain_bounds = domain_bounds;
    for (int d=0; d<DIMENSIONS; ++d) {
      if (harmonic_down[d] || ghost_down[d]) {
        ++dims[d];
        extended_domain_bounds.min[d] -= widths[d];
      }
      if (harmonic_up[d] || ghost_up[d]) {
        ++dims[d];
        extended_domain_bounds.max[d] += widths[d];
      }
    }

    // --- Create the cells
    const int size = getNumCells();
    cells = vector<Cell>(size);

    int index[DIMENSIONS];
    for (int ind=0; ind<size; ++ind) {
      linear_to_tuple(ind, index);
      vector<int> neighbors = find_adjacent_cells(index);
      // Set cell's neighbors
      if (cells[ind].adjacent_cell_id) delete [] cells[ind].adjacent_cell_id;
      cells[ind].adjacent_cell_id = new int[neighbors.size()];
      for (int i=0; i<neighbors.size(); ++i) 
        cells[ind].adjacent_cell_id[i] = neighbors[i];
      cells[ind].adjacent_cell_id_size = neighbors.size();
      // Assign this cell's type
      cells[ind].cellType = CellType::Central;
      for (int d=0; d<DIMENSIONS; ++d) {
        // Check if this cell is harmonic or a ghost
        if (index[d]==0) {
          if (!ghost_down[d] && bcs[d]==BCFlag::WRAP)
            cells[ind].cellType = CellType::Harmonic;
          else if (ghost_down[d])
            cells[ind].cellType = CellType::Ghost;
        }
        if (index[d]==dims[d]-1) {
          if (!ghost_up[d] && bcs[d]==BCFlag::WRAP) 
            cells[ind].cellType = CellType::Harmonic;
          else if (ghost_up[d])
            cells[ind].cellType = CellType::Ghost;
        }
        // Check if this cell is at the edge of the regular domain and has harmonic image cells
      }

    }

  }

  void Domain::pre_integrate() {
    // Reset time points
    lastCheck  = -1.;
    lastUpdate = -1.;
    updateDelay = 1.0e-4;
  }

  void Domain::addToCell(int id) {
    int index = getCellIndex(simData->X(id));
    Cell &cell = cells[index];
    cell.add(id);
    if (cell.cellType==CellType::Harmonic) { 
      // Create ghost particles in other cells
    }
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
    for (int d=0; d<DIMENSIONS; ++d) domain_dims[d] = 1;

    // Calculate which domain we are
    for (int d=0; d<DIMENSIONS; ++d) domain_index[d] = 0;

    // Calculate which domains are adjacent to us
    for (int d=0; d<DIMENSIONS; ++d) {
      domains_up[d] = -1;
      domains_down[d] = -1;
    }

    // Calculate the bounds of this domain
    domain_bounds = bounds;
  }

  void Domain::remake_verlet() {
    // @todo Implement this
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

  vector<int> Domain::find_adjacent_cells(int index[DIMENSIONS]) {
    // Helper arrays
    int d_index[DIMENSIONS], other_index[DIMENSIONS], linear;
    // A vector for the (linear) indices of adjacent cells
    vector<int> neighbors;

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

      // We are using harmonic cells for wrapping, so if an index is negative, we are done
      bool good = true;
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

      // Do not evaluate if we shouldn't
      if (!good) continue;
      // The cell at other_index is a neighbor that we should include
      tuple_to_linear(linear, other_index);
      neighbors.push_back(linear);
    }
    // Return
    return neighbors;
  }
}
