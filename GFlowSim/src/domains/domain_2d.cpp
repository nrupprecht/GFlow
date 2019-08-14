#include "domain_2d.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../base/forcemaster.hpp"

namespace GFlowSimulation {

  Domain2D::Domain2D(GFlow *gflow) : DomainBase(gflow) {
    // Initialize to 0.
    min_side_edge_cells[0] = 1;
    min_side_edge_cells[1] = 0;
    max_side_edge_cells[0] = 1;
    max_side_edge_cells[1] = 1;
    // Initialize to false.
    halo_cells[0] = halo_cells[1] = false;
  };

  Domain2D::~Domain2D() {
    if (cell_pointers) delete [] cell_pointers;
    if (linked_cells) delete [] linked_cells;
  }

  void Domain2D::initialize() {
    // Common initialization tasks
    DomainBase::initialize();

    // If bounds are unset, then don't make sectors. We cannot initialize if simdata is null
    if (process_bounds.vol()<=0 || simData==nullptr || simData->size()==0) return; 

    // Calculate skin depth
    RealType rho = simData->size() / process_bounds.vol();
    RealType candidate = inv_sphere_volume((target_list_size)/rho + 0.5*sphere_volume(max_small_sigma, sim_dimensions), sim_dimensions) - 2*max_small_sigma;
    skin_depth = max(static_cast<RealType>(0.5 * max_small_sigma), candidate);

    // Use max_small_sigma
    target_cell_size = min_small_cutoff = 2*max_small_sigma+skin_depth;

    // Calculate cell grid data
    calculate_domain_cell_dimensions();

    // Construct the interaction handlers for the forces
    construct();

    // The domain is now initialized.
    initialized = true;
  }

  void Domain2D::getAllWithin(int id1, vector<int>& neighbors, RealType distance) {
    // STUB
  }

  void Domain2D::removeOverlapping(RealType) {
    // STUB
  }

  void Domain2D::setSkinDepth(RealType) {

  }

  void Domain2D::setCellSize(RealType) {

  }

  void Domain2D::construct() {
    DomainBase::construct();

    // Only bother constructing under some circumstances.
    if (forceMaster==nullptr || !forceMaster->interactingTypes()) return;

    // Start timer.
    timer.start();

    // Fill up cells.
    update_cells();

    // Create stencil of offsets.
    const int stencil_size = 4;
    int stencil[stencil_size];
    stencil[0] = -dims[0]-1; // Bottom left
    stencil[1] = -1;         // Left
    stencil[2] = dims[0]-1;  // Top left
    stencil[3] = -dims[0];   // Bottom
    
    // Iterate through cells.
    for (int c=0; c<cells_size; ++c) {
      // Get the first position in the array.
      int id1 = cell_pointers[c];
      // Check if this is a real particle.
      if (id1<0) continue;
      // Go through all id1's in the first cell.
      while (id1!=-1) {
        // Is this particle "small"?
        bool is_small = (simData->Sg(id1)<=max_small_sigma);
        // Go through all the particles in the same cells.
        int id2 = linked_cells[id1];
        // Check this cell.
        check_cell(id1, id2);
        // Check surrounding cells.
        if (is_small) {
          // Go through particles in adjacent cells.
          for (int s=0; s<stencil_size; ++s) {
            // Get the number of the other cell.
            int c2 = c + stencil[s];
            // Make sure the cell is valid.
            if (c2<0 || cells_size<=c2) continue;
            // Go through all particles in the cell.
            id2 = cell_pointers[c2];
            // Check cells.
            check_cell(id1, id2);
          }
        }
        else { // "Large" particle.
          // Find which cell this is.
          int cx1 = c % dims[0];
          int cy1 = c / dims[0];
          RealType search_width = 2*simData->Sg(id1)*max_cutoffs[simData->Type(id1)]+skin_depth;
          int sx = static_cast<int>(ceil(search_width/widths[0]));
          int sy = static_cast<int>(ceil(search_width/widths[1]));
          bool wrapX = gflow->getBC(0)==BCFlag::WRAP;
          bool wrapY = gflow->getBC(1)==BCFlag::WRAP;
          // Search around the particle.
          for (int dy=-sy; dy<=sy; ++dy)
            for (int dx=-sx; dx<=sx; ++dx) {
              // Cell coords
              int cx2 = cx1 + dx, cy2 = cy1 + dy;
              // Potentially wrap cells
              if (wrapX) {
                if (cx2<min_side_edge_cells[0]) cx2 += dims[0] - min_side_edge_cells[0] - max_side_edge_cells[0];
                else if (dims[0]+min_side_edge_cells[0]<cx2) cx2 -= dims[0] - min_side_edge_cells[0] - max_side_edge_cells[0];
              }
              if (wrapY) {
                if (cy2<min_side_edge_cells[1]) cy2 += dims[1] - min_side_edge_cells[1] - max_side_edge_cells[1];
                else if (dims[1]+min_side_edge_cells[0]<cy2) cy2 -= dims[1] - min_side_edge_cells[1] - max_side_edge_cells[1];
              }
              // If a good cell, search through it.
              if ( min_side_edge_cells[0]<=cx2 && cx2<dims[0]+min_side_edge_cells[0] 
                && min_side_edge_cells[1]<=cy2 && cy2<dims[1]+min_side_edge_cells[1]) {
                int c2 = cy2*dims[0] + cx2;
                id2 = cell_pointers[c2];
                check_cell_large(id1, id2);
              }
            }
        }

        // "Increment" the cell id.
        id1 = linked_cells[id1];
      }
    }

    forceMaster->close();

    // Start timer.
    timer.stop();
  }

  void Domain2D::calculate_domain_cell_dimensions() {
    // sim_dimensions = 2
    for (int d=0; d<sim_dimensions; ++d) {
      dims[d] = static_cast<int>(process_bounds.wd(d)/target_cell_size);
      // Check that the bounds are good.
      if (dims[d]<=0) throw BadBounds();
      widths[d] = process_bounds.wd(d)/dims[d];
      inverseW[d] = 1./widths[d];

      // Add halo or ghost cells.
      dims[d] += (min_side_edge_cells[d] + max_side_edge_cells[d]);
    }

    // Compute number of cells. If we need to change it, allocate a new array.
    if (cells_size!=dims[0]*dims[1]) {
      if (cell_pointers) delete [] cell_pointers;
      cells_size = dims[0]*dims[1];
    }

    // Create the cell pointers array.
    if (cells_size>0) cell_pointers = new int[cells_size];
    else cell_pointers = nullptr;
  }


  void Domain2D::migrate_particles() {

  }

  void Domain2D::construct_halo_particles() {
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

  void Domain2D::construct_ghost_particles() {

  }

  inline void Domain2D::update_cells() {
    // Get the size.
    int size = simData->size();

    // Check if we have to remake the linked cells array.
    if (size>list_size) {
      if (linked_cells) delete [] linked_cells;
      if (size>0) linked_cells = new int[size];
      else linked_cells = nullptr;
      list_size = size;
    }

    // First "clear" cells
    for (int i=0; i<cells_size; ++i) cell_pointers[i] = -1;
    for (int i=0; i<list_size; ++i)  linked_cells [i] = -1;

    // Get data arrays
    RealType **x = simData->X();

    // Add particles to cells
    RealType shiftX = process_bounds.min[0]*inverseW[0], shiftY = process_bounds.min[1]*inverseW[1];
    for (int i=0; i<size; ++i) {
      // Compute the cell index of the particle
      int cx = static_cast<int>(x[i][0]*inverseW[0] - shiftX) + min_side_edge_cells[0];
      int cy = static_cast<int>(x[i][1]*inverseW[1] - shiftY) + min_side_edge_cells[1];
      
      // Checks -> these may not be necessary.
      if (cx<0) cx = 0;
      if (dims[0]<=cx) cx = dims[0]-1;
      if (cy<0) cy = 0;
      if (dims[1]<=cy) cy = dims[1]-1;

      // Compute the linear cell index
      int index = cy*dims[0] + cx;
      linked_cells[i] = cell_pointers[index];
      cell_pointers[index] = i;
    }
  }

  inline void Domain2D::check_cell(int id1, int id2) {
    if (id2<0) return;

    // Set up.
    RealType dx[2], radius1 = simData->Sg(id1)*max_cutoffs[simData->Type(id1)] + skin_depth, rsqr = 0;
    RealType **x = simData->X();
    RealType *sg = simData->Sg();
    // Go through particles in second cell, starting with id2.
    while (-1<id2) {
      // If the other particle is large, continue.
      if (simData->Sg(id2)<max_small_sigma) {
        // Get distance between particles.
        dx[0] = x[id1][0] - x[id2][0];
        dx[1] = x[id1][1] - x[id2][1];
        rsqr = sqr(dx[0]) + sqr(dx[1]);
        // If close enough, check if they interact.
        if (rsqr < sqr(radius1 + sg[id2])) // *max_cutoffs[simData->Type(id2)]
          pair_interaction_nw(id1, id2);
      }

      // "Increment" id2
      id2 = linked_cells[id2];
    }
  }

  inline void Domain2D::check_cell_large(int id1, int id2) {
    if (id2<0) return;

    // Set up.
    RealType dx[2], radius1 = simData->Sg(id1)*max_cutoffs[simData->Type(id1)] + skin_depth, rsqr = 0;
    RealType **x = simData->X();
    RealType *sg = simData->Sg();
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
  }

}