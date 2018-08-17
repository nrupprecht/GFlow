#include "domain2d_test.hpp"
// Other files
#include "simdata.hpp"
#include "vectormath.hpp"
#include "memoryoptimizer.hpp"

namespace GFlowSimulation {

  Domain2d::Domain2d(GFlow *gflow) : DomainBase(gflow) {
    // Make sure the simulation really is 2D
    if (DIMENSIONS!=2) {
      cout << "Test class \"Domain2D\" only works in ... you guessed it ... 2D! Fail!\n";
      exit(1);
    }

    cout << "Using Domain2d, an experimental class, as the DomainBase object.\n";
  }

  void Domain2d::initialize() {
    // First call Base's initialize function
    Base::initialize();

    // Get the simulation bounds from gflow
    bounds = Base::gflow->getBounds();
    domain_bounds = bounds;

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

    // Here is where we would usually call construct. There is no need to.
  }

  void Domain2d::pre_integrate() {
    // Reset time points
    lastCheck  = -1.;
    lastUpdate = -1.;
    updateDelay = 1.0e-4;
  }

  void Domain2d::calculateForces() {
    if (dims[0]<2 || dims[1]<2) {
      cout << "Dims[0] or Dims[1] are < 2. This is an uncovered case.\n";
      exit(1);
    }

    // --- Initialize
    
    // Load the initial cell
    load_cell(0, 0);

    // Go through all the cells
    for (int x=0; x<dims[0]; ++x) {
      for (int y=0; y<dims[1]; ++y) {
        // Load the required cell data
        if (y!=dims[1]-1) { // We don't have to load the row above if there is no row above
          int address = (y*dims[0] + x) + dims[0] + 1;
          release_cell(address);
        }

        // Get the current cell
        Cell2d &cell = getCell(x, y);
        int size = cell.size();
        // Check if there are any particles in ths cell
        if (size==0) {
          // Release the old data
          if (y!=0) { // We don't have to release the row below if there is no row below
            int address = (y*dims[0] + x) - dims[0] - 1;
            release_cell(address);
          }
          continue;
        }

        cout << "Size: " << size << endl;

        // Do interactions for this cell
        for (int i=0; i<size; ++i) 
          for (int j=i+1; j<size; ++j) {

            // Find displacements
            subtractVec(cell.x[i], cell.x[j], cell.dx[j]);
            RealType dsqr    = sqr(cell.dx[j]);
            RealType dsg2    = sqr( cell.sg[i] + cell.sg[j] );
            RealType distance = sqrt(dsqr);
            // The 1./distance will turn dx into a normal vector in the scalarMultVec multiplication
            RealType c1      = dsqr - dsg2 < 0. ? 1./distance : 0.;
            RealType magnitude = c1*(dsg2 - dsqr);
            // Find the (vectorial) force, place in cell2.dx array
            scalarMultVec(magnitude, cell.dx[j]);
            // Add force to accumulators
            plusEqVec (cell.f[i], cell.dx[j]);
            minusEqVec(cell.f[j], cell.dx[j]);
          }
        // Done with the interactions for just this cell

        // Look at all surrounding cells
        for (int y1=-1; y1<=1; ++y1) {
          // Keep cell indices in bounds
          y1 = y<0 ? y1+dims[1] : y1;
          y1 = dims[1]<=y1 ? y1-dims[1] : y1;
          for (int x1=-1; x1<=1; ++x1) {
            // Keep cell indices in bounds
            x1 = x1<0 ? x1+dims[0] : x1;
            x1 = dims[0]<=x1 ? x1-dims[0] : x1;

            // Get cell

            cout << x1 << " " << y1 << endl;

            Cell2d &cell2 = getCell(x1, y1);

            if (!cell2.loaded) load_cell(cell);

            int size2 = cell2.size();
            if (size2==0) continue;

            // Interactions:
            for (int i=0; i<size; ++i)
              for (int j=0; j<size2; ++j) {
                subtractVec(cell.x[i], cell2.x[j], cell2.dx[j]);
                RealType dsqr    = sqr(cell2.dx[j]);
                RealType dsg2    = sqr( cell.sg[i] + cell2.sg[j] );
                RealType distance = sqrt(dsqr);
                // The 1./distance will turn dx into a normal vector in the scalarMultVec multiplication
                RealType c1      = dsqr - dsg2 < 0. ? 1./distance : 0.;
                RealType magnitude = c1*(dsg2 - dsqr);
                // Find the (vectorial) force, place in cell2.dx array
                scalarMultVec(magnitude, cell2.dx[j]);
                // Add force to accumulators
                plusEqVec (cell.f[i],  cell2.dx[j]);
                minusEqVec(cell2.f[j], cell2.dx[j]);
              }
            // Done with interactions between cell, cell2
          }
        }

        // Release the old data
        if (y!=0) { // We don't have to release the row below if there is no row below
          int address = (y*dims[0] + x) - dims[0] - 1;
          release_cell(address);
        }

      }
    }

    // Finish - release the remaining cells
    for (int y=dims[1]-2; y<dims[1]; ++y)
      for (int x=0; x<dims[0]; ++x) {
        release_cell(x, y);
      }

  }

  void Domain2d::construct() {
    fill_cells();
  }

  inline void Domain2d::create_cells() {
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

    // --- Create the cells
    const int size = dims[0]*dims[1];
    cells = vector<Cell2d>(size, Cell2d());
  }

  inline void Domain2d::fill_cells() {
    // Clear cells of old data
    for (int i=0; i<cells.size(); ++i) 
      cells[i].clear();

    // --- Place particles in their cells
    int linear;
    for (int n=0; n<Base::simData->number; ++n) {
      // Calculate the (DIMENSION)-tuple cell coordinate
      int linear = getCellIndex(Base::simData->x[n]);
      cells[linear].add(n);
    }
  }

  inline void Domain2d::load_cell(int x, int y) {
    // Fetch cell
    Cell2d &cell = getCell(x, y);
    int size = cell.size();
    // Potentially allocate data arrays
    if (size>0 && !cell.loaded) {
      cell.x = alloc_array_2d<RealType>(size, DIMENSIONS);
      cell.dx = alloc_array_2d<RealType>(size, DIMENSIONS);
      cell.f = alloc_array_2d<RealType>(size, DIMENSIONS);
      cell.sg = new RealType[size];
      for (int i=0; i<size; ++i) {
        int id = cell.id_list[i];
        copyVec(Base::simData->x[id], cell.x[i]);
        zeroVec(cell.f[i]);
        cell.sg[i] = Base::simData->sg[id];
      }
    }
    // Set loaded to true
    cell.loaded = true;
  }

  inline void Domain2d::load_cell(int linear) {
    // Fetch cell
    Cell2d &cell = cells[linear];
    int size = cell.size();
    // Potentially allocate data arrays
    if (size>0 && !cell.loaded) {
      cell.x = alloc_array_2d<RealType>(size, DIMENSIONS);
      cell.dx = alloc_array_2d<RealType>(size, DIMENSIONS);
      cell.f = alloc_array_2d<RealType>(size, DIMENSIONS);
      cell.sg = new RealType[size];
      // Load data from simData to local copies in the cell
      for (int i=0; i<size; ++i) {
        int id = cell.id_list[i];
        copyVec(Base::simData->x[id], cell.x[i]);
        zeroVec(cell.f[i]);
        cell.sg[i] = Base::simData->sg[id];
      }
    }
    // Set loaded to true
    cell.loaded = true;
  }

  inline void Domain2d::load_cell(Cell2d& cell) {
    int size = cell.size();
    // Potentially allocate data arrays
    if (size>0 && !cell.loaded) {
      cell.x = alloc_array_2d<RealType>(size, DIMENSIONS);
      cell.dx = alloc_array_2d<RealType>(size, DIMENSIONS);
      cell.f = alloc_array_2d<RealType>(size, DIMENSIONS);
      cell.sg = new RealType[size];
      // Load data from simData to local copies in the cell
      for (int i=0; i<size; ++i) {
        int id = cell.id_list[i];
        copyVec(Base::simData->x[id], cell.x[i]);
        zeroVec(cell.f[i]);
        cell.sg[i] = Base::simData->sg[id];
      }
    }
    // Set loaded to true
    cell.loaded = true;
  }

  inline void Domain2d::release_cell(int x, int y) {
    // Fetch cell
    Cell2d &cell = getCell(x, y);
    int size = cell.size();
    // Potentially transfer and clear data arrays
    if (size>0 && cell.loaded) {
      dealloc_array_2d(cell.x);
      dealloc_array_2d(cell.dx);
      // Transfer forces back to simdata
      for (int i=0; i<size; ++i) {
        int id = cell.id_list[i];
        plusEqVec(Base::simData->f[id], cell.f[i]);
      }
      dealloc_array_2d(cell.f);
      delete [] cell.sg;
      // Set to null
      cell.x = cell.dx = cell.f = nullptr;
      cell.sg = nullptr;
    }
    // Set loaded to false
    cell.loaded = false;
  }

  inline void Domain2d::release_cell(int linear) {
    Cell2d &cell = cells[linear];
    int size = cell.size();
    // Potentially transfer and clear data arrays
    if (size>0 && cell.loaded) {
      dealloc_array_2d(cell.x);
      dealloc_array_2d(cell.dx);
      // Transfer forces back to simdata
      for (int i=0; i<size; ++i) {
        int id = cell.id_list[i];
        plusEqVec(Base::simData->f[id], cell.f[i]);
      }
      dealloc_array_2d(cell.f);
      delete [] cell.sg;
      // Set to null
      cell.x = cell.dx = cell.f = nullptr;
      cell.sg = nullptr;
    }
    // Set loaded to false
    cell.loaded = false;
  }

  inline int Domain2d::getCellIndex(RealType *x) {
    int index[DIMENSIONS], linear;
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) {
      index[d] = static_cast<int>((x[d] - domain_bounds.min[d])*inverseW[d]);
      if (index[d]>=dims[d]) index[d] = dims[d]-1; 
      else if (index[d]<0)   index[d] = 0;
    }
    // Convert to linear index
    linear = dims[0]*index[1] + index[0];
    // Return
    return linear;
  }

  inline Cell2d& Domain2d::getCell(int x, int y) {
    return cells[ dims[0]*y + x ];
  }

}