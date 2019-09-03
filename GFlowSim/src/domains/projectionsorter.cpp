#include "projectionsorter.hpp"

namespace GFlowSimulation {

  ProjectionSorter::ProjectionSorter(GFlow *gflow) : InteractionHandler(gflow), axis(gflow->getSimDimensions()) {
    auto_record_positions = false;
  };

  void ProjectionSorter::construct() {
    // Common tasks
    InteractionHandler::construct();

    // Start timer.
    timer.start();

    // Sort the particles
    simData->sortParticles(axis);

    RealType **x = simData->X();
    RealType *sg = simData->Sg();
    int *type = simData->Type();
    int size = simData->size();
    RealType cutoff = 2*max_small_sigma + skin_depth;
    RealType width = simulation_bounds.wd(0);
    bool wrap = true;
    RealType dX[2];
    
    for (int id1=0; id1<size; ++id1) {
      int id2 = id1+1;
      RealType radius = sg[id1]; // \todo Not necessarily the cutoff.
      RealType end_value = max_small_sigma + sg[id1] + skin_depth + x[id1][0];
      while (id2<size && x[id2][0] < end_value) {
        RealType dsqr = getDistanceSqrNoWrap(x[id1], x[id2], sim_dimensions);
        // If the particles are close enough, add them to the interaction.
        if (dsqr<sqr(sg[id1] + sg[id2] + skin_depth)) pair_interaction(id1, id2);
        // Increment.
        ++id2;
      }
      
      if (wrap && id2==size) {
        id2 = 0;
        end_value = max_small_sigma + sg[id1] + skin_depth + x[id1][0] - width;

        while (id2<size && x[id2][0] < x[id2][0] < end_value) {
          gflow->getDisplacement(x[id1], x[id2], dX);
          RealType dsqr = dX[0]*dX[0] + dX[1]*dX[1];
          if (dsqr<sqr(sg[id1] + sg[id2] + skin_depth)) pair_interaction(id1, id2, 1);
          ++id2;
        }
      }
    }

    // Record positions.
    record_positions();

    // End timer.
    timer.stop();
  }

  void ProjectionSorter::getAllWithin(int id, vector<int>& neighbors, RealType distance) {

  }

  void ProjectionSorter::getAllWithin(Vec position, vector<int>& neighbors, RealType distance) {

  }

  void ProjectionSorter::removeOverlapping(RealType factor) {

  }

}