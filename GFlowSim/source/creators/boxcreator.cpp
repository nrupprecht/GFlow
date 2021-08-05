#include <creators/boxcreator.hpp>
// Other files
#include <gflow.hpp>
#include <interactions/interaction-choice.hpp>

using namespace GFlowSimulation;

BoxCreator::BoxCreator(int argc, char **argv)
    : Creator(argc, argv), phi(0.5), radius(0.05) {
  seed = std::chrono::system_clock::now().time_since_epoch().count();
  seedGenerator(seed);
  normal_dist = std::normal_distribution<RealType>(0., 1.);
};

BoxCreator::BoxCreator(ArgParse *p)
    : Creator(p), phi(0.5), radius(0.05) {
  seed = std::chrono::system_clock::now().time_since_epoch().count();
  seedGenerator(seed);
  normal_dist = std::normal_distribution<RealType>(0., 1.);
};

void BoxCreator::seedGenerator(uint s) {
  Creator::seedGenerator(s);
  generator = std::mt19937(seed);
}

GFlow *BoxCreator::createSimulation() {
  // Seed random number generators
  srand48(time(0));

  // Values
  int number = -1; // -1 means use volume density
  int sample = 0;
  RealType width = 4.;
  RealType height = -1.;
  RealType dt = 0.001;
  RealType vsgma = 0.25;
  RealType skinDepth = -1.;
  RealType cell_size = -1;
  RealType repulsion = 1.;
  bool over_damped_flag = false;
  RealType langevin_temp = -1.;
  bool lj_flag = false;
  bool interact_flag = true;

  // Gather command line arguments
  if (parserPtr) {
    parserPtr->get("number", number);
    parserPtr->get("sample", sample);
    parserPtr->get("dt", dt);
    parserPtr->get("radius", radius);
    parserPtr->get("phi", phi);
    parserPtr->get("width", width);
    parserPtr->get("height", height);
    parserPtr->get("vsgma", vsgma);
    parserPtr->get("skinDepth", skinDepth);
    parserPtr->get("cell_size", cell_size);
    parserPtr->get("repulsion", repulsion);
    parserPtr->get("overdamped", over_damped_flag);
    parserPtr->get("langevin", langevin_temp);
    parserPtr->get("lj", lj_flag);
    parserPtr->get("interact", interact_flag);
  }

  // Create a new gflow object
  GFlow *gflow = new GFlow(sim_dimensions);
  gflow->setAllBCs(bcFlag);

  // Set the bounds of the gflow object
  Bounds bnds = gflow->getBounds();
  if (width <= 0) {
    width = 1.;
  } // In case of bad argument
  if (height <= 0) {
    height = width;
  }
  for (int d = 0; d < sim_dimensions; ++d) {
    if (d != 1) {
      bnds.min[d] = -0.5 * width;
      bnds.max[d] = 0.5 * width;
    }
    else {
      bnds.min[d] = -0.5 * height;
      bnds.max[d] = 0.5 * height;
    }
  }
  gflow->setBounds(bnds);

  // --- Set initial particle data
  // Find how many objects to use
  if (number < 0) {
    RealType vol = pow(width, sim_dimensions);
    number = phi * vol / sphere_volume(radius, sim_dimensions);
  }

  // The simdata object
  auto simData = gflow->simData;

  // Add some objects
  simData->reserve(number);

  // Holders
  Vec X(sim_dimensions), V(sim_dimensions);
  X.zero();
  V.zero();
  //RealType *X = new RealType[sim_dimensions], *V = new RealType[sim_dimensions];
  //zeroVec(V, sim_dimensions);
  // Calculate the inverse mass
  RealType im = 1.0 / (1.0 * PI * sqr(radius));
  for (int n = 0; n < number; ++n) {
    // Give some random initial positions - we will allow these to relax
    for (int d = 0; d < sim_dimensions; ++d) {
      X[d] = (drand48() - 0.5) * width;
    }
    simData->addParticle(X.data, V.data, radius, im, 0);
  }

  // --- Initialize handler
  gflow->initialize();

  // --- Handle forces
  gflow->forceMaster->setNTypes(1); // Only one type of particle
  shared_ptr<Interaction> force;
  if (lj_flag) {
    force = InteractionChoice::choose(gflow, LennardJonesToken, sim_dimensions);
  }
  else if (!interact_flag) {
    force = nullptr;
  }
  else {
    // RealType strength = repulsion * DEFAULT_HARD_SPHERE_REPULSION;
    force = InteractionChoice::choose(gflow, HardSphereToken, sim_dimensions);
    // std::dynamic_pointer_cast<Interaction>(force); //->setRepulsion(strength);
  }
  gflow->forceMaster->setInteraction(0, 0, force);

  // --- Set some parameters
  if (skinDepth > 0) {
    gflow->handler->setSkinDepth(skinDepth);
  }
  if (sample > 0) {
    gflow->handler->setSampleSize(sample);
  }

  // Relax the setup in two steps
  hs_relax(gflow, 0.1); // Make sure particles don't stop on top of one another

  // Set the integrator
  if (over_damped_flag) {
    gflow->integrator = choose_overdamped_integrator(gflow, sim_dimensions);
  }
  else if (langevin_temp >= 0) {
    gflow->integrator = new LangevinIntegrator(gflow, langevin_temp);
  }
  else {
    gflow->integrator = choose_velocity_verlet(gflow, sim_dimensions);
  }
  // Set integrator initial time step
  gflow->integrator->setDT(dt);

  // --- Set velocities
  // In case the number of particles has changed
  number = gflow->simData->number();
  for (int n = 0; n < number; ++n) {
    // Give some random velocities
    double ke = fabs(vsgma * normal_dist(generator));
    double velocity = sqrt(2 * simData->Im(n) * ke / 127.324);
    // Random normal vector
    randomNormalVec(simData->V(n), sim_dimensions);
    // Set the velocity
    scalarMultVec(velocity, simData->V(n), sim_dimensions);
    simData->Sg(n) = radius;
  }
  // Clean up
  //delete [] X;
  //delete [] V;
  // Return
  return gflow;
}
