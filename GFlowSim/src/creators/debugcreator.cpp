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
    bool lj = false;

    // Gather command line arguments
    parserPtr->get("time", time);
    parserPtr->get("radius", radius);
    parserPtr->get("skinDepth", skinDepth);
    parserPtr->get("lj", lj);

    // Create a new gflow object
    GFlow *gflow = new GFlow(sim_dimensions);
    gflow->setAllBCs(bcFlag);

    // Create an integrator
    gflow->integrator = new VelocityVerlet(gflow);

    // Set the bounds of the gflow object --- for now, just make it [0,1] in each dimension
    for (int d=0; d<sim_dimensions; ++d) {
      gflow->bounds.min[d] = 0.;
      gflow->bounds.max[d] = 1.;
    }

    // Set wrapping
    gflow->setAllBCs(BCFlag::WRAP);

    // Add some objects
    gflow->simData->reserve(2);

    // Get pointers to particle data
    SimData *simData = gflow->simData;

    // Rightwards ball
    simData->X(0, 0) = 0.3; simData->X(1, 0) = 0.7;
    simData->V(0, 0) = velocity;  simData->V(1, 0) = -velocity;
    simData->X(0, 1) = 0.5*(1.-radius); simData->X(1,1) = 0.5*(1.+radius);
    for (int d=2; d<sim_dimensions; ++d) {
      simData->X(0, d) = simData->X(1, 1) = 0.5;
      simData->V(0, d) = simData->V(1, d) = 0;
    }
    for (int n=0; n<2; ++n) {
      simData->Sg(n) = radius;
      simData->Im(n) = 1.0 / (1.0 * PI*sqr(radius)); // Density of 1
      simData->Type(n) = 0;
    }

    // Set the correct number of particles
    gflow->simData->number = 2;

    // --- Handle forces
    gflow->forceMaster->setNTypes(1);
    Interaction *force;
    if (lj) force = new LennardJones(gflow);
    else force = new HardSphere(gflow);
    gflow->forceMaster->setInteraction(0, 0, force);

    // Set skin depth
    if (skinDepth>0) gflow->domain->setSkinDepth(skinDepth);

    // Make sure all forces are zero
    gflow->simData->clearF();

    // Request some amount of time to run
    gflow->requestTime(time); // -- 10 is a stub value

    return gflow;
  }

}
