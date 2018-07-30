#include "bond_boxcreator.hpp"
// Other files
#include "harmonicbond.hpp"
#include "gflow.hpp"

namespace GFlowSimulation {

  BondBoxCreator::BondBoxCreator(int argc, char **argv) : Creator(argc, argv) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  };

  BondBoxCreator::BondBoxCreator(ArgParse *p) : Creator(p) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  };

  void BondBoxCreator::seedGenerator(uint s) {
    Creator::seedGenerator(s);
    generator = std::mt19937(seed);
  }

  GFlow* BondBoxCreator::createSimulation() {
    // Seed random number generators
    srand48(time(0));

    // Values
    int number = -1; // -1 means use volume density
    RealType dt = 0.001;
    RealType radius = 0.05;
    RealType phi = 0.5;
    RealType width = 4.;
    RealType vsgma = 0.25;
    RealType skinDepth = -1.;
    bool animate = false;;
    bool over_damped_flag = false;

    // Gather command line arguments
    parserPtr->get("number", number);
    parserPtr->get("dt", dt);
    parserPtr->get("radius", radius);
    parserPtr->get("phi", phi);
    parserPtr->get("width", width);
    parserPtr->get("vsgma", vsgma);
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
    gflow->setAllBCs(BCFlag::WRAP);

    // --- Set initial particle data
    // Find how many objects to use
    if (number<0) {
      RealType vol = pow(width, DIMENSIONS);  
      number = phi*vol/sphere_volume(radius);
    }
    // Add some objects
    RealType eqLength = 2.1*radius; // Equilibrium bond length
    number = 2*(number/2);
    SimData *simData = gflow->simData;
    simData->reserve(number);
    simData->number = number;
    // For choosing the center of a particle pair, and their displacement from the center
    RealType center[DIMENSIONS], normal[DIMENSIONS]; 
    for (int n=0; n<number/2; ++n) {
      int id1 = 2*n, id2 = 2*n+1;
      // --- Give the pair some random initial positions near each other - we will allow these to relax
      for (int d=0; d<DIMENSIONS; ++d) center[d] = (drand48()-0.5)*width;
      // Random normal vector, to displace the particle along
      randomNormalVec(normal);
      scalarMultVec(0.5*eqLength, normal);

      // Set positions
      addVec(center, normal, simData->X(id1));
      subtractVec(center, normal, simData->X(id2));
      // Put a bond between them
      HarmonicBond *bond = new HarmonicBond(gflow, id1, id2);
      bond->setSpringK(1.);
      gflow->addModifier(bond);
      // Give a radius, mass, and type
      simData->Sg(id1) = simData->Sg(id2) = radius;
      simData->Im(id1) = simData->Im(id2) = 1.0 / (1.0 * PI*sqr(radius)); // Density of 1
      simData->Type(id1) = simData->Type(id2) = 0;
    }
    // Wrap positions
    gflow->wrapPositions();

    // --- Handle forces
    gflow->forceMaster->setNTypes(1); // Only one type of particle
    Force *hard_sphere = new HardSphere(gflow);
    gflow->forceMaster->setForce(0, 0, hard_sphere);

    // Set skin depth
    if (skinDepth>0) gflow->domain->setSkinDepth(skinDepth);

    // Relax the setup
    hs_relax(gflow);
    if (over_damped_flag) gflow->integrator = new OverdampedIntegrator(gflow);
    else gflow->integrator = new VelocityVerlet(gflow);
    // Set integrator initial time step
    gflow->integrator->setDT(dt);

    // --- Set velocities - give pairs the same velocity
    for (int n=0; n<number/2; ++n) {
      // Give some random velocities
      double ke = fabs(vsgma*normal_dist(generator));
      double velocity = sqrt(2*simData->Im(n)*ke/127.324);
      // Random normal vector
      randomNormalVec(normal);
      // Set the velocity
      scalarMultVec(velocity, normal, simData->V(2*n));
      scalarMultVec(velocity, normal, simData->V(2*n+1));
    }

    return gflow;
  }

}
