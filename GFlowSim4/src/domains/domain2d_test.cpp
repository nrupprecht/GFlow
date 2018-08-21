#include "domain2d_test.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../utility/memoryoptimizer.hpp"
#include "../utility/printingutility.hpp" // For debugging

#include <assert.h>

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

    // Construct
    construct();
  }

  void Domain2d::pre_integrate() {
    // Reset time points
    lastCheck  = -1.;
    lastUpdate = -1.;
    updateDelay = 1.0e-4;
  }

  void Domain2d::calculateForces() {
    // Get the boundary conditions
    const BCFlag *bcs = gflow->getBCs();
    
    // Load the initial cell
    load_cell(0, 0);

    // Go through all the cells
    for (int x1=0; x1<dims[0]; ++x1) {
      for (int y1=0; y1<dims[1]; ++y1) {
        // Load the required cell data
        if (y1!=dims[1]-1) { // We don't have to load the row above if there is no row above
          int address = (y1*dims[0] + x1) + dims[0] + 1;
          load_cell(address);
        }

        // Get the current cell
        Cell2d &cell1 = getCell(x1, y1);
        int size1 = cell1.size();
        // Check if there are any particles in ths cell
        if (size1==0) {
          // Release the old data
          if (y1!=0) { // We don't have to release the row below if there is no row below
            int address = (y1*dims[0] + x1) - dims[0] - 1;
            release_cell(address);
          }
          continue;
        }

        // Load cell 1 if neccessary
	      if (!cell1.loaded) load_cell(cell1);

        // Do interactions for this cell
        for (int j=0; j<size1; ++j) {
          for (int i=j+1; i<size1; ++i) {
            // x[i] - x[j] --> dx[j]
      	    //subtractVec(cell1.x[i], cell1.x[j], cell1.dx[j]);

            //**
            getDisplacement(cell1.x[i], cell1.x[j], cell1.dx[j], bounds, bcs);
            //**

      	    RealType dsqr = sqr(cell1.dx[j]);
      	    RealType dsg  = cell1.sg[i] + cell1.sg[j];
      	    // The 1./distance will turn dx into a normal vector in the scalarMultVec multiplication
      	    RealType magnitude = dsqr < sqr(dsg) ? DEFAULT_HARD_SPHERE_REPULSION*(dsg - sqrt(dsqr))/sqrt(dsqr) : 0.;
      	    // Find the (vectorial) force, place in cell2.dx array
      	    scalarMultVec(magnitude, cell1.dx[j]);
      	    // Add force to accumulators
      	    plusEqVec (cell1.f[i], cell1.dx[j]);
      	    minusEqVec(cell1.f[j], cell1.dx[j]);
          }
        }
        // Done with the interactions for just this cell

        // If only one was ==1, this would be a problem.
        if (dims[0]==1 || dims[1]==1) continue;

        // Look at all surrounding cells
        for (int dy2=-1; dy2<=1; ++dy2) {
	         int y2 = y1+dy2;
          // Keep cell indices in bounds
          y2 = (y2<0 ? y2+dims[1] : y2);
          y2 = (dims[1]<=y2 ? y2-dims[1] : y2);
	  
          for (int dx2=-1; dx2<=1; ++dx2) {
	          int x2 = x1+dx2;
            // Keep cell indices in bounds
            x2 = x2<0 ? x2+dims[0] : x2;
            x2 = dims[0]<=x2 ? x2-dims[0] : x2;

      	    // Do not interact with yourself here - we already did that.
      	    if (dy2==0 && dx2==0) continue;

            // Get cell
            Cell2d &cell2 = getCell(x2, y2);
            if (!cell2.loaded) load_cell(cell2);

            // Check if there are any particles in cell 2
            int size2 = cell2.size();
            if (size2==0) continue;

            // Interactions:
            for (int i=0; i<size1; ++i)
              for (int j=0; j<size2; ++j) {

                //subtractVec(cell1.x[i], cell2.x[j], cell2.dx[j]);

                //**
                getDisplacement(cell1.x[i], cell2.x[j], cell2.dx[j], bounds, bcs);
                //**

                RealType dsqr = sqr(cell2.dx[j]);
                RealType dsg  = cell1.sg[i] + cell2.sg[j];
                RealType distance = sqrt(dsqr);

                // The 1./distance will turn dx into a normal vector in the scalarMultVec multiplication
                RealType c1      = distance - dsg < 0. ? 1./distance : 0.;
                RealType magnitude = DEFAULT_HARD_SPHERE_REPULSION*c1*(dsg - distance);
                // Find the (vectorial) force, place in cell2.dx array
                scalarMultVec(magnitude, cell2.dx[j]);
                // Add force to accumulators
                plusEqVec (cell1.f[i], cell2.dx[j]);
                minusEqVec(cell2.f[j], cell2.dx[j]);
              }
            // Done with interactions between cell, cell2
          }
        }

        // Release the old data
        if (y1!=0) { // We don't have to release the row below if there is no row below
          int address = (y1*dims[0] + x1) - dims[0] - 1;
          release_cell(address);
        }

      }
    }

    // Finish - release the remaining cells
    /*
    for (int y=dims[1]-2; y<dims[1]; ++y)
      for (int x=0; x<dims[0]; ++x) {
        release_cell(x, y);
      }
    */

    // Release *all* the cells, just in case.
    for (int i=0; i<dims[0]*dims[1]; ++i) release_cell(i);

  }

  void Domain2d::construct() {
    fill_cells();
  }

  bool Domain2d::check_needs_remake() { 
    return true; 
  }

  inline void Domain2d::create_cells() {

    //**
    cutoff = 0.25*max(bounds.wd(0), bounds.wd(1));
    //**

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
    load_cell(cell);
  }

  inline void Domain2d::load_cell(int linear) {
    // Fetch cell
    Cell2d &cell = cells[linear];
    load_cell(cell);
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
