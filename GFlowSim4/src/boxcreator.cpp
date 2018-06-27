#include "boxcreator.hpp"
#include "overdampedintegrator.hpp"

namespace GFlowSimulation {

  BoxCreator::BoxCreator(int argc, char **argv) : Creator(argc, argv) {};

  BoxCreator::BoxCreator(ArgParse *p) : Creator(p) {};

  GFlow* BoxCreator::createSimulation() {
    // Seed random number generators
    srand48(time(0));

    // Values
    RealType time = 10.;
    int number = 10;
    RealType radius = 0.05;
    bool animate = false;;
    bool over_damped_flag = false;

    // Gather command line arguments
    parserPtr->get("time", time);
    parserPtr->get("number", number);
    parserPtr->get("radius", radius);
    parserPtr->get("animate", animate);
    parserPtr->get("overdamped", over_damped_flag);

    // Create a new gflow object
    GFlow *gflow = new GFlow;

    // Create an integrator
    if (over_damped_flag) gflow->integrator = new OverdampedIntegrator(gflow);
    else gflow->integrator = new VelocityVerlet(gflow);

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
    RealType *im = gflow->simData->im;
    int *type    = gflow->simData->type;
    for (int n=0; n<number; ++n) {
      // Give some random positions and velocities
      for (int d=0; d<DIMENSIONS; ++d) {
        x[n][d] = drand48();
        v[n][d] = 0.25*(drand48() - 0.5);
      }
      sg[n] = radius;
      im[n] = 1.0 / (1.0 * PI*sqr(radius)); // Density of 1
      type[n] = 0;
    }

    // Set the correct number of particles
    gflow->simData->number = number;

    // --- Handle forces
    gflow->forceMaster->setNTypes(1);
    Force *hard_sphere = new HardSphere(gflow);
    gflow->forceMaster->setForce(0, 0, hard_sphere);

    // Make sure all forces are zero
    gflow->simData->clearF();

    // Request some amount of time to run
    gflow->requestTime(time); // -- 10 is a stub value

    return gflow;
  }

}