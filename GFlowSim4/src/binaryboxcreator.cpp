#include "binaryboxcreator.hpp"
#include "gflow.hpp"

namespace GFlowSimulation {

  BinaryBoxCreator::BinaryBoxCreator(int argc, char **argv) : Creator(argc, argv) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
    real_dist = std::uniform_real_distribution<RealType>(0., 1.);
  };

  BinaryBoxCreator::BinaryBoxCreator(ArgParse *p) : Creator(p) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
    real_dist = std::uniform_real_distribution<RealType>(0., 1.);
  };

  void BinaryBoxCreator::seedGenerator(uint s) {
    Creator::seedGenerator(s);
    generator = std::mt19937(seed);
  }

  GFlow* BinaryBoxCreator::createSimulation() {
    // Seed random number generators
    srand48(time(0));

    // Values
    RealType time = 10.;
    int number = -1; // -1 means use volume density
    RealType dt = 0.001;
    RealType radius = 0.05;
    RealType phi = 0.5;
    RealType width = 4.;
    RealType vsgma = 0.25;
    RealType portion = 0.5; // Portion of the particles that are type 0
    RealType skinDepth = -1.;
    bool animate = false;;
    bool over_damped_flag = false;

    // Gather command line arguments
    parserPtr->get("time", time);
    parserPtr->get("number", number);
    parserPtr->get("dt", dt);
    parserPtr->get("radius", radius);
    parserPtr->get("phi", phi);
    parserPtr->get("width", width);
    parserPtr->get("vsgma", vsgma);
    parserPtr->get("portion", portion);
    parserPtr->get("skinDepth", skinDepth);
    parserPtr->get("animate", animate);
    parserPtr->get("overdamped", over_damped_flag);

    // Create a new gflow object
    GFlow *gflow = new GFlow;

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
      // Choose particle type
      RealType randreal = real_dist(generator);
      if (randreal<portion) type[n] = 0;
      else                  type[n] = 1;
    }

    // --- Handle forces

    gflow->forceMaster->setNTypes(2); // Only one type of particle
    Force *hard_sphere = new HardSphere(gflow);
    gflow->forceMaster->setForce(0, 1, hard_sphere);
    gflow->forceMaster->setForce(1, 0, hard_sphere);

    // Set skin depth
    if (skinDepth>0) gflow->sectorization->setSkinDepth(skinDepth);

    // --- Create integrator

    // Relax the setup
    relax(gflow);
    if (over_damped_flag) gflow->integrator = new OverdampedIntegrator(gflow);
    else gflow->integrator = new VelocityVerlet(gflow);
    // Set integrator initial time step
    gflow->integrator->setDT(dt);

    // --- Set velocities

    RealType normal[DIMENSIONS];
    for (int n=0; n<number; ++n) {
      // Give some random velocities
      double ke = fabs(vsgma*normal_dist(generator));
      double velocity = sqrt(2*im[n]*ke/127.324);
      // Random normal vector
      randomNormalVec(normal);
      // Set the velocity
      scalarMultVec(velocity, normal, v[n]);
      sg[n] = radius;
    }

    // --- Request the desired amount of time

    gflow->resetAllTimes();
    gflow->requestTime(time);

    return gflow;
  }

}