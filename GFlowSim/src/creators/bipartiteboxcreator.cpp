#include "bipartiteboxcreator.hpp"
#include "../gflow.hpp"

namespace GFlowSimulation {

  BipartiteBoxCreator::BipartiteBoxCreator(int argc, char **argv) : Creator(argc, argv), phi(0.5), width(4.), radius(0.05) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
    real_dist = std::uniform_real_distribution<RealType>(0., 1.);
  };

  BipartiteBoxCreator::BipartiteBoxCreator(ArgParse *p) : Creator(p), phi(0.5), width(4.), radius(0.05) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
    real_dist = std::uniform_real_distribution<RealType>(0., 1.);
  };

  void BipartiteBoxCreator::seedGenerator(uint s) {
    Creator::seedGenerator(s);
    generator = std::mt19937(seed);
  }

  GFlow* BipartiteBoxCreator::createSimulation() {
    // Checks
    if (radius<=0 || phi<=0 || width<=0) return nullptr;

    // Seed random number generators
    srand48(time(0));

    // Values
    int number = -1; // -1 means use volume density
    int ntypes = 2;  
    RealType dt = 0.001;
    RealType vsgma = 0.25;
    RealType portion = 0.5; // Portion of the particles that are type 0
    RealType skinDepth = -1.;
    RealType repulsion = 1.;
    bool animate = false;;
    bool over_damped_flag = false;
    bool lj_flag = false;

    // Gather command line arguments
    if (parserPtr) {
      parserPtr->get("number", number);
      parserPtr->get("ntypes", ntypes);
      parserPtr->get("dt", dt);
      parserPtr->get("radius", radius);
      parserPtr->get("phi", phi);
      parserPtr->get("width", width);
      parserPtr->get("vsgma", vsgma);
      parserPtr->get("portion", portion);
      parserPtr->get("skinDepth", skinDepth);
      parserPtr->get("repulsion", repulsion);
      parserPtr->get("animate", animate);
      parserPtr->get("overdamped", over_damped_flag);
      parserPtr->get("lj", lj_flag);
    }

    // Create a new gflow object
    GFlow *gflow = new GFlow(sim_dimensions);
    gflow->setAllBCs(bcFlag);

    // Set the bounds of the gflow object
    if (width<=0) width = 1.; // In case of bad argument
    for (int d=0; d<sim_dimensions; ++d) {
      gflow->bounds.min[d] = -0.5*width;
      gflow->bounds.max[d] =  0.5*width;
    }

    // --- Set initial particle data

    // Find how many objects to use
    if (number<0) {
      RealType vol = pow(width, sim_dimensions);  
      number = phi*vol/sphere_volume(radius, sim_dimensions);
    }
    // Add some objects
    gflow->simData->reserve(number);
    gflow->simData->addParticle(number);
    // Get pointers to particle data
    SimData *simData = gflow->simData;
    for (int n=0; n<number; ++n) {
      // Give some random initial positions - we will allow these to relax
      for (int d=0; d<sim_dimensions; ++d) simData->X(n, d) = (drand48()-0.5)*width;
      simData->Sg(n) = radius;
      simData->Im(n) = 1.0 / (1.0 * PI*sqr(radius)); // Density of 1
      // Choose particle type
      RealType randreal = real_dist(generator);
      if (ntypes>2) {
        simData->Type(n) = static_cast<int>(randreal*ntypes);
      }
      else {
        if (randreal<portion) simData->Type(n) = 0;
        else                  simData->Type(n) = 1;
      }
    }

    // --- Handle forces

    gflow->forceMaster->setNTypes(ntypes);
    Interaction *force;
    if (lj_flag) {
      if (sim_dimensions==2) force = new LennardJones_VerletPairs_2d(gflow);
      else throw false;
    }
    else {
      HardSphere *hs_force;
      if (sim_dimensions==2) hs_force = new HardSphere_VerletPairs_2d(gflow);
      else if (sim_dimensions==3) hs_force = new HardSphere_VerletPairs_3d(gflow);
      else throw false;
      hs_force->setRepulsion(repulsion*DEFAULT_HARD_SPHERE_REPULSION);
      force = hs_force;
    }
    for (int t1=0; t1<ntypes; ++t1) 
      for (int t2 = 0; t2<ntypes; ++t2) 
        if (t1!=t2)
          gflow->forceMaster->setInteraction(t1, t2, force);

    // Set skin depth
    if (skinDepth>0) gflow->handler->setSkinDepth(skinDepth);

    // --- Create integrator

    // Relax the setup
    relax(gflow);
    if (over_damped_flag) gflow->integrator = new OverdampedIntegrator(gflow);
    else gflow->integrator = new VelocityVerlet(gflow);
    // Set integrator initial time step
    gflow->integrator->setDT(dt);

    // --- Set velocities
    number = simData->number();
    for (int n=0; n<number; ++n) {
      // Give some random velocities
      double ke = fabs(vsgma*normal_dist(generator));
      double velocity = sqrt(2*simData->Im(n)*ke/127.324);
      // Random normal vector
      randomNormalVec(simData->V(n), sim_dimensions);
      // Set the velocity
      scalarMultVec(velocity, simData->V(n), sim_dimensions);
      simData->Sg(n) = radius;
    }

    return gflow;
  }

}
