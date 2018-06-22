#include "boxcreator.hpp"

namespace GFlowSimulation {

  BoxCreator::BoxCreator(int argc, char **argv) : Creator(argc, argv) {};

  GFlow* BoxCreator::createSimulation() {
    // Seed random number generators
    srand48(time(0));

    // Values
    RealType time = 10;
    int number = 10;
    RealType radius = 0.05;
    bool animate = false;;

    // Gather command line arguments
    ArgParse parser(argc, argv);
    parser.get("time", time);
    parser.get("number", number);
    parser.get("radius", radius);
    parser.get("animate", animate);

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
    gflow->simData->reserve(number);

    // --- Set all particle data

    // Get pointers to particle data
    RealType **x = gflow->simData->x;
    RealType **v = gflow->simData->v;
    RealType *sg = gflow->simData->sg;
    int *type    = gflow->simData->type;
    for (int n=0; n<number; ++n) {
      // Give some random positions and velocities
      for (int d=0; d<DIMENSIONS; ++d) {
        x[n][d] = drand48();
        v[n][d] = 0.25*(drand48() - 0.5);
      }
      sg[n] = radius;
      type[n] = 0;
    }

    // Set the correct number of particles
    gflow->simData->number = number;

    // Make sure all forces are zero
    gflow->simData->clearF();

    // Request some amount of time to run
    gflow->requestTime(time); // -- 10 is a stub value

    // Give data objects to data master
    if (animate) 
      gflow->dataMaster->addDataObject(new PositionData(gflow));

    return gflow;
  }

}