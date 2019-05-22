#include "domain_2d.hpp"
// Other files
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  Domain2D::Domain2D(GFlow *gflow) : DomainBase(gflow) {
    // Initialize to 0.
    max_side_edge_cells[0] = max_side_edge_cells[1] = 0;
    min_side_edge_cells[0] = min_side_edge_cells[1] = 0;
    // Initialize to false.
    halo_cells[0] = halo_cells[1] = false;
  };

  void Domain2D::initialize() {
    // Common initialization tasks
    Base::initialize();

    // Get the simulation bounds from gflow
    bounds = Base::gflow->getBounds();
    // Get the bounds from the gflow object - for now assumes this is the only domain, so bounds==domain_bounds
    domain_bounds = gflow->getBounds();

    // If bounds are unset, then don't make sectors. We cannot initialize if simdata is null
    if (domain_bounds.vol()<=0 || simData==nullptr) return; 

    // Calculate the maxiumu "small sigma"
    calculate_max_small_sigma();

    // Use max_small_sigma
    target_cell_size = min_small_cutoff = 2*max_small_sigma+skin_depth;

    // Calculate cell grid data
    calculate_domain_cell_dimensions();

    // Construct the interaction handlers for the forces
    construct();
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
    make_cells();

    // Create stencil of offsets.
    const int stencil_size = 4;
    int stencil[stencil_size];
    stencil[0] = -dims[0]-1; // Bottom left
    stencil[1] = -1;         // Left
    stencil[2] = dims[0]-1;  // Top left
    stencil[3] = -dims[0];   // Bottom

    // Get data arrays.
    RealType **x = simData->X();
    RealType *sg = simData->Sg();
    RealType dx, dy, rsqr;

    // Iterate through cells.
    for (int c=0; c<cells_size; ++c) {
      int cell_start = cell_pointers[c];
      // Check if the cell is empty
      if (cell_start<0) continue;

      // Go through all particles in the cell
      for (int p=linked_cells[cell_start]; p!=-1; p=linked_cells[p]) {
        // If the particle is a small particle
        if (sg[p]<=max_small_sigma) {
          // Look at the other particles in the same cell
          for (int q = linked_cells[p]; q!=-1; q=linked_cells[q]) {
            // If the other particle is a big particle, continue.
            if (sg[q]>max_small_sigma) continue;
            dx = x[p][0] - x[q][0];
            dy = x[p][1] - x[q][1];
            rsqr = dx*dx + dy*dy;
            // If close enough, check if they interact.
            if (rsqr < sqr(sg[p] + sg[q] + skin_depth))
              pair_interaction(p, q);
          }

          // Look at the surrounding cells
          for (int s=0; s<stencil_size; ++s) {
            // Get neighbor cell.
            int d = c + stencil[s];
            if (d<0 || cells_size<=d) continue;
            // Look at all particles in the cell
            int cell_start_2 = cell_pointers[d];
            // If the cell is empty
            if (cell_start_2<0) continue;
            // Go through all the particles in the cell.
            for (int q = linked_cells[cell_start_2]; q!=-1; q=linked_cells[q]) {
              // If the other particle is a big particle, continue.
              if (sg[q]>max_small_sigma) continue;
              dx = x[p][0] - x[q][0];
              dy = x[p][1] - x[q][1];
              rsqr = dx*dx + dy*dy;
              // If close enough, check if they interact.
              if (rsqr < sqr(sg[p] + sg[q] + skin_depth))
                pair_interaction(p, q);
            }
          }

        }
        // The particle is a large particle, and needs to look in a larger region.
        else {
          // Need a special stencil for the large particle
          
        }
      }



    }


  }

  inline void Domain2D::calculate_domain_cell_dimensions() {
    // sim_dimensions = 2
    for (int d=0; d<sim_dimensions; ++d) {
      dims[d] = static_cast<int>(domain_bounds.wd(d)/target_cell_size);
      // Check that the bounds are good.
      if (dims[d]<=0) throw BadBounds();
      widths[d] = domain_bounds.wd(d)/dims[d];
      inverseW[d] = 1./widths[d];

      // Add halo or ghost cells.
      dims[d] += (min_side_edge_cells[d] + max_side_edge_cells[d]);
    }

    // Compute number of cells
    cells_size = dims[0]*dims[1];

    // Create the cell pointers array.
    cell_pointers = new int[cells_size];
  }

  inline void Domain2D::make_cells() {
    // Get the size.
    int size = simData->size();
    // Check if we have to remake the linked cells array.
    if (size>list_size) {
      if (linked_cells) delete [] linked_cells;
      linked_cells = new int[size];
    }

    // First "clear" cells
    for (int i=0; i<cells_size; ++i) cell_pointers[i] = -1;
    for (int i=0; i<list_size; ++i) linked_cells[i]  = -1;

    // Get data arrays
    RealType **x = simData->X();
    RealType *sg = simData->Sg();

    // Add particles to cells
    for (int i=0; i<size; ++i) {
      // Compute the cell index of the particle
      int cx = (x[i][0] - domain_bounds.min[0])*inverseW[0];
      int cy = (x[i][1] - domain_bounds.min[1])*inverseW[1];
      cx += min_side_edge_cells[0];
      cy += min_side_edge_cells[1];
      // Compute the linear cell index
      int linear = cy*dims[0] + cx;

      linked_cells[i] = cell_pointers[linear];
      cell_pointers[linear] = i;

      // Create halo cell if needed
      if (halo_cells[0]) {

        // If a "small particle," then we only need to check if the particle falls the boundary.
        if (sg[i]<=max_small_sigma) {
          if (cx==0 && cy==0) {
            // Add three halo particles, top, right, and top right
            RealType width = domain_bounds.wd(0);
            RealType height = domain_bounds.wd(1);

          }
          else if (cx==0) {

          }
          else if (cy==0) {

          }
        }
        // Have to check more things for large particles
        else {

        }
      }
      
    }
  }


}