#include "linecrystalcreator.hpp"
// Other files
#include "../gflow.hpp"

namespace GFlowSimulation {

  LineCrystalCreator::LineCrystalCreator(int argc, char **argv) : Creator(argc, argv), radius(0.05), lambda(0.08), width(4.), height(4.), N(10), M(10) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  };

  LineCrystalCreator::LineCrystalCreator(ArgParse *p) : Creator(p), radius(0.05), lambda(0.08), width(4.), height(4.), N(10), M(10) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  };

  void LineCrystalCreator::seedGenerator(uint s) {
    Creator::seedGenerator(s);
    generator = std::mt19937(seed);
  }

  GFlow* LineCrystalCreator::createSimulation() {
    // Check dimensions
    if (DIMENSIONS!=2) {
      cout << "LineCrystalCreator expects DIMENSIONS=2. DIMENSIONS=" << DIMENSIONS << ". Exiting.";
      exit(0);
    }

    // --- Command line parameters
    RealType vsgma = 0.25;

    // Gather command line arguments
    if (parserPtr) {
      parserPtr->get("vsgma", vsgma);
      parserPtr->get("lambda", lambda);
    }

    // Create a new gflow object
    GFlow *gflow = new GFlow;
    gflow->setAllBCs(bcFlag);

    // Find the proper width and height
    width = 1.1 * N*(2. + gamma(lambda/radius))*radius;
    height = lambda*M;

    // Set bounds
    gflow->bounds.min[0] = 0;
    gflow->bounds.max[0] = width;
    gflow->bounds.min[1] = 0;
    gflow->bounds.max[1] = height;

    // --- Set initial particle data
    // Find how many objects to use
    int number = N*M;

    // Add some objects
    SimData *simData = gflow->simData;
    simData->reserve(number);
    //-- simData->number = number;
    RealType X[DIMENSIONS], V[DIMENSIONS];
    zeroVec(V); 
    // --- Set up a crystal
    RealType im = 1.0 / (1.0 * PI*sqr(radius));
    RealType dx = width / N;
    RealType dy = lambda;
    V[1] = 0.;
    for (int m=0; m<M; ++m) {
      RealType offset = m%2==0 ? 0. : 0.5*dx;
      X[1] = m*dy;
      for (int n=0; n<N; ++n) {
        X[0] = offset + n*dx;
        // Random velocity
        double v = (drand48() > 0.5 ? 1. : -1) * sqrt(2*im*vsgma*fabs(normal_dist(generator))/127.324);
        V[0] = v;
        // Insert particle
        simData->addParticle(X, V, radius, im, 0); //--
      }
    }

    // --- Initialize domain
    gflow->domain->initialize();

    // --- Handle forces
    gflow->forceMaster->setNTypes(1); // Only one type of particle
    Interaction *force = new HardSphere(gflow);
    gflow->forceMaster->setInteraction(0, 0, force);

    // We have to create an integrator
    gflow->integrator = new VelocityVerlet(gflow);

    return gflow;
  }

  RealType LineCrystalCreator::gamma(RealType y) {
    return y<2.*radius ? sqrt(1.-sqr(0.5*y/radius)) : 0;
  }


}