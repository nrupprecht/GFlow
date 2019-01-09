#include "linkedcells.hpp"
// Other files
#include "../gflow.hpp"
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../interactionhandlers/verletlist.hpp"
#include "../base/forcemaster.hpp"

namespace GFlowSimulation {

  LinkedCells::LinkedCells(GFlow *gflow) : DomainBase(gflow), heads(nullptr), list(nullptr), ncells(0), nlist(0) {};

  void LinkedCells::initialize() {
    // This initialization is borrowed from the DomainTest object.

    // Common initialization tasks
    Base::initialize();

    // Get the simulation bounds from gflow
    bounds = Base::gflow->getBounds();
    // Get the bounds from the gflow object - for now assumes this is the only domain, so bounds==domain_bounds
    domain_bounds = gflow->getBounds();

    // If bounds are unset, then don't make sectors
    if (domain_bounds.vol()<=0) return;

    // We cannot initialize if simdata is null
    if (simData==nullptr) return;

    // Find average sigma
    RealType sigma = 0, max_sigma = 0;
    for (int n=0; n<Base::simData->size(); ++n) {
      if (Base::simData->Type()<0) continue;
      RealType s = Base::simData->Sg(n);
      sigma += s;
      if (s>max_sigma) max_sigma = s;
    }
    sigma /= Base::simData->number();
    // Threshhold sigma is between average and maximum sigma
    RealType threshold = 0.5*(sigma + max_sigma), max_under = sigma;
    if (threshold!=sigma) {
      for (int n=0; n<Base::simData->size(); ++n) {
        if (Base::simData->Type()<0) continue;
        RealType s = Base::simData->Sg(n);
        if (s<threshold && max_under<s) max_under = s;
      }
    }
    max_small_sigma = 1.025*max_under;
    RealType target_width = (2*max_small_sigma+skin_depth);

    // Use max_small_sigma
    for (int d=0; d<sim_dimensions; ++d) {
      dims[d] = static_cast<int>(domain_bounds.wd(d)/target_width);
      if (dims[d]<=0) throw BadBounds(); // Make sure bin numbers are positive
      widths[d] = domain_bounds.wd(d)/dims[d];
      inverseW[d] = 1./widths[d];
    }
    minCutoff = 2*max_small_sigma + skin_depth;

    //! @todo Processors might need to communicate with one another about what they chose at this point

    // Create the cells
    create_cells();

    // Construct the interaction handlers for the forces
    construct();
  }

  void LinkedCells::exchange_particles() {
    // Stub
  }

  void LinkedCells::getAllWithin(int, RealType, vector<int>&) {
    // Stub
  }

  void LinkedCells::removeOverlapping(RealType) {
    // Stub
  }

  void LinkedCells::construct() {
    // Domain base common tasks
    DomainBase::construct();

    // Update particles in the cells
    update_cells();

    // Set a "maximum reasonable distance"
    RealType max_reasonable = sqr(0.9*bounds.wd(0));

    // Get particle data
    RealType *sg = Base::simData->Sg();
    RealType **X = Base::simData->X();
    // Get the boundary conditions
    const BCFlag *bcs = gflow->getBCs();

    // Go through all cells
    for (int y=0; y<dims[1]; ++y) {
      for (int x=0; x<dims[0]; ++x) {
        int index = y*dims[0] + x; // Index of the cell
        int id1 = heads[index];    // ID of the particle
        int id2;
        while (id1!=-1) {

          if (sg[id1]<=max_small_sigma) {
            // Check particles in the same cell
            id2 = list[id1];
            while (id2!=-1) {
              // Check distance
              RealType r2 = getDistanceSqrNoWrap(X[id1], X[id2], sim_dimensions);
              if (r2 < sqr(sg[id1] + sg[id2] + skin_depth))
                pair_interaction(id1, id2);

              // Increment
              id2 = list[id2];
            }

            // Look at surrounding sectors
            int ix = x-1, iy = y-1;
            if (correct_index(ix, iy))
              check_cell(iy*dims[0]+ix, id1);
            ix = x-1; iy = y;
            if (correct_index(ix, iy))
              check_cell(iy*dims[0]+ix, id1);
            ix = x-1; iy = y+1;
            if (correct_index(ix, iy))
              check_cell(iy*dims[0]+ix, id1);
            ix = x; iy = y-1;
            if (correct_index(ix, iy))
              check_cell(iy*dims[0]+ix, id1);
          }
          else {
            // We have to check more than just the surrounding sectors.

            // Calculate sweep "radius"
            RealType search_width = 2*sg[id1]+skin_depth;
            int sx = static_cast<int>(ceil(search_width/widths[0]));
            int sy = static_cast<int>(ceil(search_width/widths[1]));

            // Sweep over the square
            for (int dy=-sy; dy<=sy; ++sy)
              for (int dx=-sx; dx<=sx; ++sx) {
                int x2 = x+dx, y2 = y+dy;
                if (correct_index(x2, y2)) {
                  int index = y2*dims[0] + x2;
                  id2 = heads[index];
                  while (id2!=-1) {
                    if (id1==id2 || sg[id2]>sg[id1]) continue;
                    RealType r2 = getDistanceSqrNoWrap(X[id1], X[id2], sim_dimensions);
                    if (r2 < sqr(sg[id1] + sg[id2] + skin_depth) || max_reasonable<r2)
                      pair_interaction(id1, id2);
                    // Increment
                    id2 = list[id2];
                  }
                }
              }

          }
          // Increment
          id1 = list[id1];
        }

      }
    }

  }

  void LinkedCells::setCellSize(RealType) {

  }

  inline void LinkedCells::create_cells() {
    if (heads) delete [] heads;
    if (list)  delete [] list;
    // Calculate number of cells
    ncells = 1;
    for (int d=0; d<sim_dimensions; ++d) ncells *= dims[d];
    // Create head list
    heads = new int[ncells];
    // Create list
    nlist = Base::simData->size();
    list = new int[nlist];
  }

  inline void LinkedCells::update_cells() {
    // Reset heads, list
    for (int i=0; i<ncells; ++i) heads[i] = -1;
    for (int i=0; i<nlist; ++i)  list [i] = -1;

    // Assign particles to sectors
    RealType **X = simData->X();
    for (int i=0; i<simData->size(); ++i) {
      // Assumes two dimensions
      int x = static_cast<int>((X[i][0] - domain_bounds.min[0])*inverseW[0]);
      int y = static_cast<int>((X[i][1] - domain_bounds.min[1])*inverseW[1]);
      int index = y*dims[0] + x;
      // Add particle to its sector
      list[i] = heads[index];
      heads[index] = i;
    }

  }

  inline void LinkedCells::check_cell(int index, int id1) {
    // Get data
    RealType **x = simData->X();
    RealType *sg = simData->Sg();
    // Set a "maximum reasonable distance"
    RealType max_reasonable = sqr(0.9*bounds.wd(0));
    
    int id2 = heads[index];
    while (id2!=-1) {
      // If the other particle is a large particle, it will take care of this interaction
      if (sg[id2]>max_small_sigma) continue;
      // Look for distance between particles
      RealType r2 = getDistanceSqrNoWrap(x[id1], x[id2], sim_dimensions);
      // Possibly interact
      if (r2 < sqr(sg[id1] + sg[id2] + skin_depth) || max_reasonable<r2)
        pair_interaction(id1, id2);
      // Increment
      id2 = list[id2];
    }
  }

  inline bool LinkedCells::correct_index(int& ix, int& iy) {
    // Get the boundary conditions
    const BCFlag *bcs = gflow->getBCs();
    // Correct X
    if (ix<0) {
      if (bcs[0]==BCFlag::WRAP) ix += dims[0];
      else return false;
    }
    else if (dims[0]<=ix) {
      if (bcs[0]==BCFlag::WRAP) ix -= dims[0];
      else return false;
    }
    // Correct Y
    if (iy<0) {
      if (bcs[1]==BCFlag::WRAP) iy += dims[1];
      else return false;
    }
    else if (dims[1]<=iy) {
      if (bcs[1]==BCFlag::WRAP) iy -= dims[1];
      else return false;
    }
    // Return success
    return true;
  }

}