#include "debugcreator.hpp"

namespace GFlowSimulation {

  // Constructor
  DebugCreator::DebugCreator(int argc, char** argv) : Creator(argc, argv) {};

  DebugCreator::DebugCreator(ArgParse *p) : Creator(p) {};

  // Create simulation
  GFlow* DebugCreator::createSimulation() {
    // Seed random number generators
    srand48(time(0));

    // Values
    RealType time = 3.;
    RealType radius = 0.05;
    RealType velocity = 0.25;
    RealType skinDepth = -1.;

    // Gather command line arguments
    parserPtr->get("time", time);
    parserPtr->get("radius", radius);
    parserPtr->get("skinDepth", skinDepth);

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
    gflow->simData->reserve(2);

    // Get pointers to particle data
    RealType **x   = gflow->simData->x;
    RealType **v   = gflow->simData->v;
    RealType *sg   = gflow->simData->sg;
    RealType *im   = gflow->simData->im;
    int      *type = gflow->simData->type;

    // Rightwards ball
    x[0][0] = 0.3; x[1][0] = 0.7;
    v[0][0] = velocity;  v[1][0] = -velocity;
    x[0][1] = 0.5*(1.-radius); x[1][1] = 0.5*(1.+radius);
    for (int d=2; d<DIMENSIONS; ++d) {
      x[0][d] = x[1][d] = 0.5;
      v[0][d] = v[1][d] = 0;
    }
    for (int n=0; n<2; ++n) {
      sg[n] = radius;
      im[n] = 1.0 / (1.0 * PI*sqr(radius)); // Density of 1
      type[n] = 0;
    }

    // Set the correct number of particles
    gflow->simData->number = 2;

    // --- Handle forces
    gflow->forceMaster->setNTypes(1);
    Force *hard_sphere = new HardSphere(gflow);
    gflow->forceMaster->setForce(0, 0, hard_sphere);

    // Set skin depth
    if (skinDepth>0) gflow->sectorization->setSkinDepth(skinDepth);

    // Make sure all forces are zero
    gflow->simData->clearF();

    // Request some amount of time to run
    gflow->requestTime(time); // -- 10 is a stub value

    return gflow;
  }

}
