#include "boxcreator.hpp"

namespace GFlowSimulation {

  BoxCreator::BoxCreator(int argc, char **argv) : Creator(argc, argv) {};

  GFlow* BoxCreator::createSimulation() {
    // Create a new gflow object
    GFlow *gflow = new GFlow;

    // Create an integrator
    gflow->integrator = new VelocityVerlet(gflow);

    // Set the bounds of the gflow object --- for now, just make it [0,1] in each dimension
    for (int d=0; d<DIMENSIONS; ++d) {
      gflow->bounds.min[d] = 0.;
      gflow->bounds.max[d] = 1.;
    }

    // Set wrapping
    gflow->setAllWrap(true);

    // Add some objects
    int number = 10;
    gflow->simData->reserve(number);

    // Give some random positions and velocities
    RealType **x = gflow->simData->x;
    RealType **v = gflow->simData->v;
    for (int n=0; n<number; ++n)
      for (int d=0; d<DIMENSIONS; ++d) {
        x[n][d] = drand48();
        v[n][d] = drand48();
      }

    // Set the correct number of particles
    gflow->simData->number = number;

    // Make sure all forces are zero
    gflow->simData->clearF();

    // Make sure everything is initialized (this isn't really necessary - 
    // initialization happens at the beginning of every run)
    gflow->initialize();

    // Give data objects to data master
    gflow->dataMaster->addDataObject(new PositionData(gflow));

    return gflow;
  }

}