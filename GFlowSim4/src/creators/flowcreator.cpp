#include "flowcreator.hpp"

namespace GFlowSimulation {

  FlowCreator::FlowCreator(int argc, char **argv) : Creator(argc, argv), phi(0.5), width(4.), radius(0.05) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  };

  FlowCreator::FlowCreator(ArgParse *p) : Creator(p), phi(0.5), width(4.), radius(0.05) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  };

  void FlowCreator::seedGenerator(uint s) {
    Creator::seedGenerator(s);
    generator = std::mt19937(seed);
  }

  GFlow* FlowCreator::createSimulation() {
    // Seed random number generators
    srand48(time(0));

    // Values
    int number = -1; // -1 means use volume density
    int sample = 0;
    RealType dt = 0.001;
    RealType length = 4.;
    RealType vsgma = 0.25;
    RealType skinDepth = -1.;
    RealType cell_size = -1;
    RealType repulsion = 1.;
    bool animate = false;
    bool over_damped_flag = false;
    RealType langevin_temp = -1.;
    bool lj_flag = false;

    // Gather command line arguments
    if (parserPtr) {
      parserPtr->get("number", number);
      parserPtr->get("sample", sample);
      parserPtr->get("dt", dt);
      parserPtr->get("radius", radius);
      parserPtr->get("phi", phi);
      parserPtr->get("width", width);
      parserPtr->get("vsgma", vsgma);
      parserPtr->get("skinDepth", skinDepth);
      parserPtr->get("cell_size", cell_size);
      parserPtr->get("repulsion", repulsion);
      parserPtr->get("animate", animate);
      parserPtr->get("overdamped", over_damped_flag);
      parserPtr->get("langevin", langevin_temp);
      parserPtr->get("lj", lj_flag);
    }

    // Create a new gflow object
    GFlow *gflow = new GFlow;
    gflow->setAllBCs(BCFlag::REPL); // All boundaries but the first are repulsive
    // Make the first dimension wrapping. This is the direction of flow.
    gflow->setBC(0, BCFlag::WRAP);

    // Set the bounds of the gflow object
    if (width<=0)  width = 1.; // In case of bad argument
    if (length<=0) length = 1.;
    // Set length
    gflow->bounds.min[0] = -0.5*length;
    gflow->bounds.max[0] =  0.5*length;
    // Set remaining widths
    for (int d=1; d<DIMENSIONS; ++d) {
      gflow->bounds.min[d] = -0.5*width;
      gflow->bounds.max[d] =  0.5*width;
    }

    // --- Set initial particle data
    // Find how many objects to use
    if (number<0) {
      RealType vol = gflow->bounds.vol();
      number = phi*vol/sphere_volume(radius);
    }
    // Add some objects
    SimData *simData = gflow->simData;
    simData->reserve(number);
    //-- simData->number = number;
    RealType X[DIMENSIONS], V[DIMENSIONS]; 
    zeroVec(V); // Give them zero velocity (for now)
    for (int n=0; n<number; ++n) {
      // Give some random initial positions - we will allow these to relax
      X[0] = (drand48()-0.5)*length;
      for (int d=1; d<DIMENSIONS; ++d) {
        X[d] = (drand48()-0.5)*width;
      }
      RealType im = 1.0 / (1.0 * PI*sqr(radius));
      simData->addParticle(X, V, radius, im, 0);
    }
    // --- Initialize domain
    gflow->domain->initialize();

    // --- Handle forces
    gflow->forceMaster->setNTypes(1); // Only one type of particle
    Interaction *force;
    if(lj_flag) {
      auto *LJ = new LennardJones(gflow);
      LJ->setStrength(repulsion*DEFAULT_LENNARD_JONES_STRENGTH);
      force = LJ;
    }
    else {
      auto *HS = new HardSphere(gflow);
      HS->setRepulsion(repulsion*DEFAULT_HARD_SPHERE_REPULSION);
      force = HS;
    }
    gflow->forceMaster->setInteraction(0, 0, force);

    // --- Set some parameters
    if (skinDepth>0) gflow->domain->setSkinDepth(skinDepth);
    if (cell_size>0) gflow->domain->setCellSize(cell_size);
    if (sample>0) gflow->domain->setSampleSize(sample);

    // Relax the setup in two steps
    hs_relax(gflow, 0.1); // 1) To make sure particles don't stop on top of one another
    relax(gflow, 0.15);   // 2) To let particles relax naturally
    // Set the integrator
    if (over_damped_flag) gflow->integrator = new OverdampedIntegrator(gflow);
    else if (langevin_temp>=0) gflow->integrator = new LangevinIntegrator(gflow, langevin_temp);
    else gflow->integrator = new VelocityVerlet(gflow);
    // Set integrator initial time step
    gflow->integrator->setDT(dt);

    // --- Set velocities
    /*
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
    */

    // --- Add a flow object
    gflow->addModifier(new Flow(gflow));

    return gflow;
  }

}