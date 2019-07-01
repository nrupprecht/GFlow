#include "projectionsorter.hpp"

namespace GFlowSimulation {

  ProjectionSorter::ProjectionSorter(GFlow *gflow) : InteractionHandler(gflow), axis(gflow->getSimDimensions()) {};

  void ProjectionSorter::construct() {
    // Common tasks
    InteractionHandler::construct();

    // Sort the particles
    simData->sortParticles(axis);

    RealType **x = simData->X();
    RealType *sg = simData->Sg();
    int *type = simData->Type();
    int size = simData->size();
    RealType cutoff = 2*max_small_sigma + skin_depth;

    Vec dX(sim_dimensions);

    for (int id1=0; id1<size; ++id1) {
      int id2 = id1+1;
      while (id2<size && x[id2][0] - x[id1][0] < cutoff) {
        RealType dsqr = getDistanceSqrNoWrap(x[id1], x[id2], sim_dimensions);

        // If the particles are close enough, add them to the interaction.
        if (dsqr<sqr(sg[id1] + sg[id2] + skin_depth)) pair_interaction(id1, id2);

        // Increment.
        ++id2;
      }

    }

  }

}