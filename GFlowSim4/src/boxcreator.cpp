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
    int number = -1; // -1 means use volume density
    RealType radius = 0.05;
    RealType phi = 0.5;
    RealType width = 4.;
    bool animate = false;;
    bool over_damped_flag = false;

    // Gather command line arguments
    parserPtr->get("time", time);
    parserPtr->get("number", number);
    parserPtr->get("radius", radius);
    parserPtr->get("phi", phi);
    parserPtr->get("width", width);
    parserPtr->get("animate", animate);
    parserPtr->get("overdamped", over_damped_flag);

    // Create a new gflow object
    GFlow *gflow = new GFlow;

    // Use overdamped integrator for relaxation
    gflow->integrator = new OverdampedIntegrator(gflow);

    // Set the bounds of the gflow object
    if (width<=0) width = 1.; // In case of bad argument
    for (int d=0; d<DIMENSIONS; ++d) {
      gflow->bounds.min[d] = -0.5*width;
      gflow->bounds.max[d] =  0.5*width;
    }

    // Set wrapping
    gflow->setAllWrap(true);

    // --- Set initial particle data
    // Find how many objects to use
    if (number<0) {
      RealType vol = pow(width, DIMENSIONS);  
      number = phi*vol/sphere_volume(radius);
    }
    // Add some objects
    gflow->simData->reserve(number);
    gflow->simData->number = number;
    // Get pointers to particle data
    RealType **x = gflow->simData->x;
    RealType **v = gflow->simData->v;
    RealType *sg = gflow->simData->sg;
    RealType *im = gflow->simData->im;
    int *type    = gflow->simData->type;
    for (int n=0; n<number; ++n) {
      // Give some random initial positions - we will allow these to relax
      for (int d=0; d<DIMENSIONS; ++d) x[n][d] = (drand48()-0.5)*width;
      sg[n] = radius;
      im[n] = 1.0 / (1.0 * PI*sqr(radius)); // Density of 1
      type[n] = 0;
    }

    // --- Handle forces
    gflow->forceMaster->setNTypes(1); // Only one type of particle
    Force *hard_sphere = new HardSphere(gflow);
    gflow->forceMaster->setForce(0, 0, hard_sphere);

    // Make sure all forces are zero
    gflow->simData->clearF();

    // Relax simulation
    gflow->requestTime(0.25); // Should be long enough
    gflow->run();

    // Set new integrator
    delete [] gflow->integrator;
    if (over_damped_flag) gflow->integrator = new OverdampedIntegrator(gflow);
    else gflow->integrator = new VelocityVerlet(gflow);

    // --- Set velocities
    for (int n=0; n<number; ++n) {
      // Give some random positions and velocities
      for (int d=0; d<DIMENSIONS; ++d) v[n][d] = 0.25*(drand48() - 0.5);
      sg[n] = radius;
    }

    // Request time
    gflow->resetAllTimes();
    gflow->requestTime(time);

    return gflow;
  }

}