#include "boxcreator.hpp"
#include "gflow.hpp"

namespace GFlowSimulation {

  BoxCreator::BoxCreator(int argc, char **argv) : Creator(argc, argv), phi(0.5), width(4.), radius(0.05) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  };

  BoxCreator::BoxCreator(ArgParse *p) : Creator(p), phi(0.5), width(4.), radius(0.05) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  };

  void BoxCreator::seedGenerator(uint s) {
    Creator::seedGenerator(s);
    generator = std::mt19937(seed);
  }

  GFlow* BoxCreator::createSimulation() {
    // Seed random number generators
    srand48(time(0));

    // Values
    int number = -1; // -1 means use volume density
    RealType dt = 0.001;
    RealType vsgma = 0.25;
    RealType skinDepth = -1.;
    RealType repulsion = 1.;
    bool animate = false;
    bool over_damped_flag = false;
    bool lj_flag = false;

    // Gather command line arguments
    if (parserPtr) {
      parserPtr->get("number", number);
      parserPtr->get("dt", dt);
      parserPtr->get("radius", radius);
      parserPtr->get("phi", phi);
      parserPtr->get("width", width);
      parserPtr->get("vsgma", vsgma);
      parserPtr->get("skinDepth", skinDepth);
      parserPtr->get("repulsion", repulsion);
      parserPtr->get("animate", animate);
      parserPtr->get("overdamped", over_damped_flag);
      parserPtr->get("lj", lj_flag);
    }

    // Create a new gflow object
    GFlow *gflow = new GFlow;
    gflow->setAllBCs(bcFlag);

    // Set the bounds of the gflow object
    if (width<=0) width = 1.; // In case of bad argument
    for (int d=0; d<DIMENSIONS; ++d) {
      gflow->bounds.min[d] = -0.5*width;
      gflow->bounds.max[d] =  0.5*width;
    }

    // --- Set initial particle data
    // Find how many objects to use
    if (number<0) {
      RealType vol = pow(width, DIMENSIONS);  
      number = phi*vol/sphere_volume(radius);
    }
    // Add some objects
    SimData *simData = gflow->simData;
    simData->reserve(number);
    simData->number = number;
    for (int n=0; n<number; ++n) {
      // Give some random initial positions - we will allow these to relax
      for (int d=0; d<DIMENSIONS; ++d) simData->X(n, d) = (drand48()-0.5)*width;
      simData->Sg(n) = radius;
      simData->Im(n) = 1.0 / (1.0 * PI*sqr(radius)); // Density of 1
      simData->Type(n) = 0;
    }

    // --- Handle forces
    gflow->forceMaster->setNTypes(1); // Only one type of particle
    Force *force;
    if(lj_flag) {
      force = new LennardJones(gflow);
      dynamic_cast<LennardJones*>(force)->setStrength(
        repulsion*DEFAULT_LENNARD_JONES_STRENGTH
      );
    }
    else {
      force = new HardSphere(gflow);
      dynamic_cast<HardSphere*>(force)->setRepulsion(
        repulsion*DEFAULT_HARD_SPHERE_REPULSION
      );
    }
    gflow->forceMaster->setForce(0, 0, force);

    // Set skin depth
    if (skinDepth>0) gflow->domain->setSkinDepth(skinDepth);

    // Relax the setup
    hs_relax(gflow);
    if (over_damped_flag) gflow->integrator = new OverdampedIntegrator(gflow);
    else gflow->integrator = new VelocityVerlet(gflow);
    // Set integrator initial time step
    gflow->integrator->setDT(dt);

    // --- Set velocities
    RealType normal[DIMENSIONS];
    for (int n=0; n<number; ++n) {
      // Give some random velocities
      double ke = fabs(vsgma*normal_dist(generator));
      double velocity = sqrt(2*simData->Im(n)*ke/127.324);
      // Random normal vector
      randomNormalVec(normal);
      // Set the velocity
      scalarMultVec(velocity, normal, simData->V(n));
      simData->Sg(n) = radius;
    }

    return gflow;
  }

}
