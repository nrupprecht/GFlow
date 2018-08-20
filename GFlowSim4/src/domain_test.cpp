#include "domain.hpp"
// Other files
#include "gflow.hpp"
#include "simdata.hpp"
#include "vectormath.hpp"
#include "verletlist.hpp"
#include "forcemaster.hpp"

namespace GFlowSimulation {

  inline void Domain::load_cell(Cell& cell) {
    // Check if anything needs to be done
    if (cell.loaded || cell.size()==0) return;

    // Always have a multiple of simd_data_size worth of space for particles
    cell.capacity = simd_data_size * ceil( static_cast<float>(cell.id_list.size())/simd_data_size );

    // We will store a "structure of array," { x1, x2, ... , y1, y2, ..., etc }
    cell.x    = new RealType[cell.capacity*DIMENSIONS];
    cell.f    = new RealType[cell.capacity*DIMENSIONS];
    cell.mask = new int[cell.capacity];

    // Get the coordinates
    for (int i=0; i<cell.size(); ++i) {
      int id = cell.id_list[i];
      for (int d=0; d<DIMENSIONS; ++d) {
        cell.x[i+d*cell.capacity] = simData->x[id][d];
        cell.f[i+d*cell.capacity] = 0.;
      }
    }

    // Set the mask and sigma
    int i;
    for (i=0; i<cell.size(); ++i) {
      cell.mask[i] = 0x1;
      cell.sg[i] = simData->sg[cell.id_list[i]];
    }
    // Set the rest of the mask
    for (; i<cell.capacity; ++i) cell.mask[i] = 0x0;

    // Set loaded to true
    cell.loaded = true;
  }

  inline void Domain::release_cell(Cell& cell) {
    // Check if anything needs to be done
    if (!cell.loaded || cell.size()==0) return;

    // Transfer forces back to particles
    for (int i=0; i<cell.id_list.size(); ++i) {
      int id = cell.id_list[i];
      
    }

    // Clean up arrays
    delete [] cell.x;
    cell.x = nullptr;
    delete [] cell.f;
    cell.f = nullptr;
    delete [] cell.sg;
    cell.sg = nullptr;
    delete [] cell.mask;
    cell.mask = nullptr;

    // Set loaded to false
    cell.loaded = false;
  }

}