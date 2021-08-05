#include <base/creator.hpp>
// Other files
#include <interactions/interaction-choice.hpp>

using namespace GFlowSimulation;

Creator::Creator(int ac, char **av)
    : argc(ac), argv(av), simBounds(Bounds(2)), ourParser(true), sim_dimensions(2) {
  parserPtr = new ArgParse(ac, av);
};

Creator::Creator(ArgParse *p)
    : argc(p->getArgc()),
      argv(p->getArgv()),
      simBounds(Bounds(2)),
      parserPtr(p),
      ourParser(false),
      sim_dimensions(2) {};

Creator::~Creator() {
  if (ourParser) {
    delete parserPtr;
  }
}

void Creator::setSeed(uint s) {
  seed = s;
}

unsigned Creator::getSeed() {
  return seed;
}

void Creator::seedGenerator(uint s) {
  seed = s;
};

void Creator::setBCFlag(BCFlag b) {
  bcFlag = b;
}

void Creator::setDimensions(int d) {
  sim_dimensions = d;
}

void Creator::hs_relax(GFlow *gflow, RealType time, bool relax_integrator, HeadNode *head) {
  // Check for valid object
  if (gflow == nullptr) {
    return;
  }

  // Get sim dimensions
  int sim_dimensions = gflow->getSimDimensions();

  // Save gflow's data
  auto integrator = gflow->integrator; // Save integrator
  auto master = gflow->forceMaster; // Save old master
  auto interactions = gflow->interactions; // Save old forces

  // Use hard sphere forces
  int ntypes = gflow->getNTypes();
  auto hsForce = InteractionChoice::choose(gflow, HardSphereToken, sim_dimensions);

  // New force master - has only hard sphere forces.
  ForceMaster forceMaster(gflow, ntypes);
  Integrator *rx_integrator = nullptr;
  if (relax_integrator) {
    rx_integrator = choose_relax_integrator(gflow, sim_dimensions);
  }
  else {
    rx_integrator = choose_overdamped_integrator(gflow, sim_dimensions);
  }

  // Give gflow new data
  gflow->interactions.clear(); // Clear old forces
  gflow->integrator = rx_integrator;
  gflow->forceMaster = &forceMaster; // Give it to gflow

  // All particles interact as hard spheres
  forceMaster.setInteraction(hsForce);

  // Make sure all forces are zero
  gflow->simData->clearF();
  // Relax simulation
  gflow->requestTime(time);
  gflow->setRunMode(RunMode::INIT);
  gflow->run(); // GFlow run calls initialize, so all the base objects' pointers will be correct.
  // Reset times
  gflow->dataMaster->resetTimer(); // So the setup does not count towards the run time.
  gflow->resetAllTimes();

  // Reset GFlow
  gflow->integrator = integrator;
  gflow->forceMaster = master;
  gflow->interactions = interactions;
  gflow->initialize(); // Correct all the base objects' pointers.

  // Clear forces
  gflow->simData->clearF();

  // Clean up
  delete rx_integrator;
}

void Creator::relax(class GFlow *gflow, RealType time, HeadNode *head) {
  // Check for valid object
  if (gflow == nullptr) {
    return;
  }

  // Get sim dimensions
  int sim_dimensions = gflow->getSimDimensions();

  // Use overdamped integrator for relaxation
  Integrator *integrator = gflow->integrator; // Save integrator
  OverdampedIntegratorBase *rx_integrator = choose_relax_integrator(gflow, sim_dimensions);
  gflow->integrator = rx_integrator;

  // Parse a head node to set up the relaxation.
  if (head) {
    // Create a parser
    TreeParser parser(head);
    // Add a heading.
    parser.addHeadingOptional("StepDelay");
    parser.addHeadingOptional("DampingConstant");
    // Gather parameters
    int step_delay;
    if (parser.firstArg("StepDelay", step_delay)) {
      rx_integrator->setStepDelay(step_delay);
    }
    RealType x = 0;
    if (parser.firstArg("DampingConstant", x)) {
      rx_integrator->setDamping(x);
    }
    if (parser.firstArg("MinDT", x)) {
      rx_integrator->setMinDT(x);
    }
    if (parser.firstArg("MaxDT", x)) {
      rx_integrator->setMaxDT(x);
    }
  }

  // Make sure all forces are zero
  gflow->simData->clearF();

  // Relax simulation
  gflow->requestTime(time); // Should be long enough
  gflow->setRunMode(RunMode::INIT);
  gflow->run();
  gflow->dataMaster->resetTimer(); // So the setup does not count towards the run time
  gflow->resetAllTimes();

  // Reset integrator
  gflow->integrator = integrator;
  // Destroy relax integrator.
  delete rx_integrator;
}

void Creator::correct_global_ids(GFlow *gflow) {
  auto topology = gflow->getTopology();
  int size = gflow->simData->size_owned();
  int rank = topology->getRank();
  auto simData = gflow->getSimData();

  // Gather the number of particles on each processor.
  vector<int> array(topology->getNumProc());
  MPIObject::mpi_allgather(size, array);

  // Rank 0 never needs to correct global ids. It does need to have its next_global_id updated though.
  if (rank > 0) {
    int shift = std::accumulate(array.begin(), array.begin() + rank, 0);
    simData->shift_global_ids(shift);
    // Correct global indices for groups.
    for (auto group : gflow->global_id_reliant) {
      group->shift_global_ids(shift);
    }
  }

  // Since any processor can create particles, we update next_global_id so that it points past the global
  // IDs of current particles, and each processor has a different next_global_id. From now on, when a
  // processor creates a particle, it increment next_global_id by num_proc, so that we never have global id
  // conflicts. This may lead to certain global ids being unused, but that doesn't actually matter.
  int total = std::accumulate(array.begin(), array.end(), 0);
  simData->next_global_id = total + 1 + rank;
  simData->d_global_id = topology->getNumProc();;
}

void Creator::fix_particle_velocities(const shared_ptr<SimData> &simData) {
  for (auto &fix : particle_fixers) {
    // Get local id of the particle
    int id = simData->getLocalID(fix.global_id);
    if (0 <= id && id < simData->size_owned()) {
      // Set the initial velocity of the particle.
      copyVec(fix.velocity, simData->V(id));
    }
  }
}

void Creator::clear_particle_fixers() {
  particle_fixers.clear();
}

