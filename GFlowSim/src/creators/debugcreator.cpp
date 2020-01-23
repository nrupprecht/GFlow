#include "debugcreator.hpp"
// Other files
#include "../interactions/interaction-choice.hpp"

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
    gflow->integrator = choose_velocity_verlet(gflow, sim_dimensions);

    // Set the bounds of the gflow object --- for now, just make it [0,1] in each dimension
    for (int d=0; d<sim_dimensions; ++d) {
      gflow->bounds.min[d] = -0.5;
      gflow->bounds.max[d] = 0.5;
    }

    // Set wrapping
    gflow->setAllBCs(BCFlag::WRAP);

    // Add some objects
    gflow->simData->reserve(2);
    gflow->simData->addParticle(2);

    // Get pointers to particle data
    auto simData = gflow->simData;

    // Rightwards ball
    zeroVec(simData->X(0), sim_dimensions); 
    zeroVec(simData->X(1), sim_dimensions);
    simData->X(0, 0) = -0.3; 
    simData->X(1, 0) =  0.3;
    simData->X(0, 1) = -0.5*radius;
    simData->X(1, 1) =  0.5*radius;
    zeroVec(simData->V(0), sim_dimensions); 
    zeroVec(simData->V(1), sim_dimensions);
    simData->V(0, 0) = velocity;  simData->V(1, 0) = -velocity;
    for (int n=0; n<2; ++n) {
      simData->Sg(n) = radius;
      simData->Im(n) = 1./sphere_volume(radius, sim_dimensions);
      simData->Type(n) = 0;
    }

    // --- Handle forces
    gflow->forceMaster->setNTypes(1);
    shared_ptr<Interaction> force;
    if (lj) force = InteractionChoice::choose(gflow, LennardJonesToken, sim_dimensions);
    else    force = InteractionChoice::choose(gflow, HardSphereToken, sim_dimensions);
    gflow->forceMaster->setInteraction(0, 0, force);

    // Set skin depth
    if (skinDepth>0) gflow->handler->setSkinDepth(skinDepth);

    // Make sure all forces are zero
    gflow->simData->clearF();

    // Request some amount of time to run
    gflow->requestTime(time);

    return gflow;
  }

}
