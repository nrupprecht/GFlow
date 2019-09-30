#include "domain_2d_cells.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../base/forcemaster.hpp"

namespace GFlowSimulation {

  Domain2DCells::Domain2DCells(GFlow *gflow) : DomainBase(gflow) {
    // Initialize to 0.
    dim_shift_down[0] = 1;
    dim_shift_down[1] = 0;
    dim_shift_up[0] = 1;
    dim_shift_up[1] = 1;
  };

  void Domain2DCells::traversePairs(PairFunction body) {
    auto x = simData->X();
    auto sg = simData->Sg();
    RealType dx[2], rsqr;
    auto type = simData->Type();

    for (const auto &c1 : cells)
      for (auto p=c1.particle_ids.begin(); p!=c1.particle_ids.end(); ++p) {
        // The id of the first particle.
        int id1 = *p;
        RealType radius1 = sg[id1]*max_cutoffs[type[id1]] + skin_depth;
        // Search through the same cell
        for (auto q = p+1; q!=c1.particle_ids.end(); ++q) {
          // The id of the second particle.
          int id2 = *q;
          // Get distance between particles.
          dx[0] = x[id1][0] - x[id2][0];
          dx[1] = x[id1][1] - x[id2][1];
          rsqr = sqr(dx[0]) + sqr(dx[1]);
          // If close enough, check if they interact.
          if (rsqr < sqr(radius1 + sg[id2] + skin_depth))
            body(id1, id2, 0, radius1, sg[id2], rsqr);
        }

        // Search through adjacent cells.
        if (sg[id1] < max_small_sigma) {
          for (auto c2 : c1.adjacent) check_cell(id1, c2, body);
        }
        else {
          // Large particle.
        }
      }
  }

  void Domain2DCells::calculate_domain_cell_dimensions() {
    // sim_dimensions = 2
    for (int d=0; d<sim_dimensions; ++d) {
      dims[d] = static_cast<int>(process_bounds.wd(d)/target_cell_size);
      // Check that the bounds are good.
      if (dims[d]<=0) throw BadBounds();
      widths[d] = process_bounds.wd(d)/dims[d];
      inverseW[d] = 1./widths[d];

      // Add halo or ghost cells.
      dims[d] += (dim_shift_down[d] + dim_shift_up[d]);
    }
  }

  void Domain2DCells::construct_halo_particles() {
    // No particles.
    if (simData->size()==0) return;
    // Get wrapping information.
    bool wrapX = gflow->getBC(0)==BCFlag::WRAP;
    bool wrapY = gflow->getBC(1)==BCFlag::WRAP;
    // Check if anything needs to be done.
    if (!wrapX && !wrapY) return;
    // Set up. We don't try to get the position pointer from simdata, because it may be invalidated by adding particles, which can
    // cause a resizing event.
    RealType *x = nullptr; 
    RealType marginX = widths[0];
    RealType marginY = widths[1];
    RealType left   = simulation_bounds.min[0] + marginX;
    RealType right  = simulation_bounds.max[0] - marginX; 
    RealType bottom = simulation_bounds.min[1] + marginY;
    RealType top    = simulation_bounds.max[1] - marginX;
    Vec dx(2); dx[0] = simulation_bounds.wd(0);
    Vec dy(2); dy[1] = simulation_bounds.wd(1);
    Vec dz(2); dz[0] = dx[0]; dz[1] = dy[1];
    // We actually need to put this here, not in the for loop, because the loop changes simData->_size by adding halo particles.
    int size = simData->size();
    // Check each particle.
    for (int i=0; i<size; ++i) {
      // Get the position of the particle
      x = simData->X(i);
      // Halo in the (+1, 0) direction.
      if (wrapX &&x[0]<=left) simData->createHaloOf(i, dx);
      // Check in Y direction
      if (wrapY) { 
        if (x[1]<=bottom) {
          // Halo in the (0, +1) direction.
          simData->createHaloOf(i, dy);

          // Special cases: corner cells
          if (wrapX) {
            // Also needs a halo in the (+1, +1) direction.
            if(x[0]<=left) simData->createHaloOf(i, dz);
            // Also needs a halo in the (-1, +1) direction.
            else if (right<x[0]) {
              dz[0] = -dz[0];
              simData->createHaloOf(i, dz);
              dz[0] = -dz[0];
            }
          }
        }
      }
    }

    // Tell the location of the first halo particle.
    simData->setFirstHalo(size);
    simData->setFirstGhost(simData->size());
  }

  void Domain2DCells::structure_updates() {
    // Fill up cells.
    update_cells();
  }

  void Domain2DCells::create_cells() {
    // Compute number of cells.
    const int size = dims[0]*dims[1];
    if (size==cells.size()) return;
    cells = vector<Cell>(size, Cell());

    // Create stencil of offsets.
    const int stencil_size = 4;
    int stencil[stencil_size];
    stencil[0] = -dims[0]-1; // Bottom left
    stencil[1] = -1;         // Left
    stencil[2] = dims[0]-1;  // Top left
    stencil[3] = -dims[0];   // Bottom

    for (int y=0; y<dims[1]; ++y) {
      for (int x=0; x<dims[0]; ++x) {
        // Linear index of this cell.
        int c = y*dims[0] + x;
        // Add adjacent cells.
        if (0<=x-1) {
          if (0<=y-1) cells[c].adjacent.push_back(&cells[c + stencil[0]]);
          cells[c].adjacent.push_back(&cells[c + stencil[1]]);
          if (y+1<dims[1]) cells[c].adjacent.push_back(&cells[c + stencil[2]]);
        }
        if (0<=y-1) cells[c].adjacent.push_back(&cells[c + stencil[3]]); 
      }
    }
  }

  inline void Domain2DCells::update_cells() {
    // Get the size.
    int size = simData->size();

    // First "clear" cells
    for (auto &cell : cells) cell.particle_ids.clear();

    // Get data arrays
    RealType **x = simData->X();

    // Add particles to cells
    RealType shiftX = process_bounds.min[0]*inverseW[0], shiftY = process_bounds.min[1]*inverseW[1];
    for (int i=0; i<size; ++i) {
      // Compute the cell index of the particle
      int cx = static_cast<int>(x[i][0]*inverseW[0] - shiftX) + dim_shift_down[0];
      int cy = static_cast<int>(x[i][1]*inverseW[1] - shiftY) + dim_shift_down[1];
      
      // Checks -> these may not be necessary.
      if (cx<0) cx = 0;
      if (dims[0]<=cx) cx = dims[0]-1;
      if (cy<0) cy = 0;
      if (dims[1]<=cy) cy = dims[1]-1;

      // Compute the linear cell index
      int index = cy*dims[0] + cx;
      cells[index].particle_ids.push_back(i);
    }
  }

  inline void Domain2DCells::check_cell(int id1, const Cell *cell, PairFunction body) {
    // Set up.
    RealType dx[2], radius1 = simData->Sg(id1) + skin_depth, rsqr = 0;
    auto x = simData->X();
    auto sg = simData->Sg();

    for (const auto id2 : cell->particle_ids) {
      // If the other particle is a large particle, it will take care of this interaction
      if (sg[id2]>max_small_sigma) continue;
      // Get distance between particles.
      dx[0] = x[id1][0] - x[id2][0];
      dx[1] = x[id1][1] - x[id2][1];
      rsqr = sqr(dx[0]) + sqr(dx[1]);
      // If close enough, check if they interact.
      if (rsqr < sqr(radius1 + sg[id2] + skin_depth))
        body(id1, id2, 0, radius1, sg[id2], rsqr);
      else if (rsqr > sqr(0.9*simulation_bounds.wd(0))) // *max_cutoffs[simData->Type(id2)]
        body(id1, id2, 1, radius1, sg[id2], rsqr);
    }
  }

  inline void Domain2DCells::check_cell_large(int id1, int id2, PairFunction) {
    if (id2<0) return;

    // Set up.
    RealType dx[2], radius1 = simData->Sg(id1)*max_cutoffs[simData->Type(id1)] + skin_depth, rsqr = 0;
    RealType **x = simData->X();
    RealType *sg = simData->Sg();

    /*
    // Go through particles in second cell, starting with id2.
    while (-1<id2) {
      // If the other particle is large, continue.
      if (simData->Sg(id2)>simData->Sg(id1) || (simData->Sg(id2)==simData->Sg(id1) && id1<id2) || id1==id2) {
        id2 = linked_cells[id2];
        continue;
      };
      // Get distance between particles.
      dx[0] = x[id1][0] - x[id2][0];
      dx[1] = x[id1][1] - x[id2][1];
      rsqr = sqr(dx[0]) + sqr(dx[1]);
      // If close enough, check if they interact.
      if (rsqr < sqr(radius1 + sg[id2]) || rsqr > sqr(0.9*simulation_bounds.wd(0))) // *max_cutoffs[simData->Type(id2)]
        pair_interaction(id1, id2);

      // "Increment" id2
      id2 = linked_cells[id2];
    }
    */
  }

}