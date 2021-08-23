#include "gflow.hpp"
// Other files
#include "allbaseobjects.hpp"
#include "allbodies.hpp"
#include "alldomains.hpp"
#include "allparallelobjects.hpp"

using namespace GFlowSimulation;

GFlow::GFlow(int dims)
    : simulation_bounds(dims), sim_dimensions(dims) {
  // Set up basic objects. The integrator will be created by the creator.
  simData = std::make_shared<SimData>(this);
  integrator = nullptr;
  handler = new Domain(this);
  dataMaster = new DataMaster(this);
  forceMaster = new ForceMaster(this);
  topology = new KDTreeTopology(this);
  // Set up boundary conditions
  boundaryConditions = new BCFlag[sim_dimensions];
  // Set wrapping to true by default
  setAllBCs(BCFlag::WRAP);
  // Initialize all base objects, so they have valid pointers to one another
  initialize();
}

GFlow::~GFlow() {
  delete integrator;
  delete handler;
  delete dataMaster;
  delete forceMaster;
  delete topology;
  integrator = nullptr;
  handler = nullptr;
  dataMaster = nullptr;
  forceMaster = nullptr;
  topology = nullptr;
  delete[] boundaryConditions;
  boundaryConditions = nullptr;
}

bool GFlow::initialize() {
  bool non_null = true;
  // --- Initialize all the subobjects
  simData->initialize();

  if (topology) {
    topology->initialize();
  }
  else {
    non_null = false;
  }

  if (integrator) {
    integrator->initialize();
  }
  else {
    non_null = false;
  }

  if (forceMaster) {
    forceMaster->initialize(); // Should go before handler
  }
  else {
    non_null = false;
  }

  if (handler) {
    handler->initialize();
  }
  else {
    non_null = false;
  }

  if (dataMaster) {
    dataMaster->initialize();
  }
  else {
    non_null = false;
  }

  for (auto &md : modifiers) {
    if (md) {
      md->initialize();
    }
    else {
      non_null = false;
    }
  }

  for (auto &it: interactions) {
    if (it) {
      it->initialize();
    }
    else {
      non_null = false;
    }
  }
  // Return whether pointers were non-null
  return non_null;
}

void GFlow::initializeBase(Base *base) {
  // Make sure we aren't handed a null pointer
  if (base == nullptr) {
    return;
  }
  // Set the number of dimensions
  base->sim_dimensions = sim_dimensions;
  // Give pointer to this GFlow object
  base->gflow = this;
  // Set other objects
  base->simData = simData;
  base->integrator = integrator;
  base->handler = handler;
  base->dataMaster = dataMaster;
  base->forceMaster = forceMaster;
  base->topology = topology;
  // Set vectors
  // base->modifiersPtr = &modifiers;
  // base->interactionsPtr = &interactions;
}

void GFlow::run(long double rt) {
  // If a parameter was passed in, it is the requested time
  if (rt > 0) {
    requested_time = rt;
  }

  // Record this request
  total_requested_time += requested_time;

  // If there are no particles, we are done
  if (simData->number() == 0) {
    elapsed_time += requested_time;
    total_time += requested_time;
    return;
  }

  // Make sure we have initialized everything
  if (!initialize()) {
    // Some object was null
    throw UnexpectedNullPointer("Error: Some object was null at GFlow initialization.");
  }
  // Check that simdata has good arrays.
  if (simData->X().isnull() || simData->V().isnull() || simData->F().isnull()) {
    throw UnexpectedNullPointer("Some array in simdata was null that shouldn't be.");
  }

  // Set run mode, if it is currently "idle"
  if (runMode == RunMode::IDLE) {
    setRunMode(RunMode::SIM);
  }

  // GFlow local timers. Timed objects' timers are cleared by call to pre_integrate.
  bonded_timer.clear_timer();
  body_timer.clear_timer();
  mpi_exchange_timer.clear_timer();
  mpi_ghost_timer.clear_timer();

  // --> Pre-integrate
  running_ = true;
  terminate_ = false;
  elapsed_time = 0;
  iter = 0;
  // Set up all objects.
  // Put modifiers first so if anything is modified, the other objects can react to it.
  for (const auto& m : modifiers) {
    m->pre_integrate();
  }
  for (const auto& b : bodies) {
    b->pre_integrate();
  }
  simData->pre_integrate();
  integrator->pre_integrate();
  for (const auto& it : additional_integrators) {
    it->pre_integrate();
  }
  handler->pre_integrate();
  dataMaster->pre_integrate();
  forceMaster->pre_integrate();

  // Run Timer - for printing updates. We use a timer, not a timed object, so that if timed objects are turned off,
  // we still can collect timing data.
  Timer timer;
  timer.start();

  // Do integration for the requested amount of time
  while (running_ && requested_time > 0) {

    // --> Pre-step
    if (!modifiers.empty()) {
      modifier_timer.start_timer();
      for (const auto& m : modifiers) {
        m->pre_step();
      }
      modifier_timer.stop_timer();
    }
    integrator->pre_step();
    for (const auto& it : additional_integrators) {
      it->pre_step();
    }
    dataMaster->pre_step();
    handler->pre_step();

    // --> Pre-force
    if (!modifiers.empty()) {
      modifier_timer.start_timer();
      for (const auto& m : modifiers) {
        m->pre_forces();
      }
      modifier_timer.stop_timer();
    }
    integrator->pre_forces(); // -- This is where VV first half kick happens (if applicable)
    for (const auto& it : additional_integrators) {
      it->pre_forces();
    } // -- First half kick could also happen here.
    dataMaster->pre_forces();
    handler->pre_forces();   // -- This is where resectorization / verlet list creation might happen. Particle exchange can also happen here.

    // --- Do interactions

    // Clear force buffers
    clearForces();

    // This involves velocities, so it should be done before ghost particles, but can be done before or after clear forces.
    reflectPositions();
    repulsePositions(); // This needs to be done after clear forces.
    attractPositions(); // This does too.

    // Update ghost particles. This the positions of particles on this processor that are ghosts on other processors back to
    // the other processors. This should be done after VV second half kick happens.
    if (1 < topology->getNumProc() && !_handler_remade && use_ghosts_) {
      topology->send_ghost_updates();
      topology->recv_ghost_updates();
    }

    // Calculate interactions and forces.
    if (use_forces_) {
      // Calculate interactions.
      forceMaster->interact();
      forceMaster->interact_ghosts();

      // Calculate bonded interactions.
      if (!bondedInteractions.empty()) {
        bonded_timer.start_timer();
        for (const auto& bd : bondedInteractions)
          bd->interact();
        bonded_timer.stop_timer();
      }
    }

    /*
    // Finish updating ghost particles.
    if (!_handler_remade && use_ghosts_) simData->finishGhostParticleUpdates();
    // Calculate ghost particle forces
    if (use_forces_) forceMaster->interact_ghosts();
    */

    // Update bodies - this may include forces.
    if (!bodies.empty()) {
      body_timer.start_timer();
      for (const auto& b : bodies) {
        b->correct();
      }
      body_timer.stop_timer();
    }

    // Do modifier removal
    handleModifiers();

    // --> Post-forces
    // This is where modifiers should do forces (if they need to)
    if (!modifiers.empty()) {
      modifier_timer.start_timer();
      for (const auto& m : modifiers) {
        m->post_forces();
      }
      modifier_timer.stop_timer();
    }
    // Continue with normal order of updates.
    simData->post_forces();
    integrator->post_forces();                 // -- This is where VV second half kick happens (if applicable)
    for (const auto& it : additional_integrators) {
      it->post_forces();
    }
    dataMaster->post_forces();
    handler->post_forces();

    // --> Post-step
    if (requested_time <= elapsed_time)
      terminate();
    if (!modifiers.empty()) {
      modifier_timer.start_timer();
      for (const auto& m : modifiers) {
        m->post_step();
      }
      modifier_timer.stop_timer();
    }
    integrator->post_step();
    for (const auto& it : additional_integrators) {
      it->post_step();
    }
    dataMaster->post_step();
    handler->post_step();

    // Timer updates
    ++iter;
    RealType dt = integrator->getTimeStep();
    elapsed_time += dt;
    total_time += dt;
    // Check for bad numerical precision
    if (total_time - dt == total_time) {
      std::cout << "Loss of precision. Stopping simulation.\n";
      terminate();
    }
    // Possibly print updates to the screen or to a file.
    if (print_updates && topology->getRank() == 0 && runMode == RunMode::SIM &&
        static_cast<int>((elapsed_time - dt) / update_interval_) < static_cast<int>((elapsed_time) / update_interval_)) {
      RealType live_ratio = elapsed_time / timer.current();
      (*monitor_) << "Simulation time: " << static_cast<int>(elapsed_time) << "\t";
      (*monitor_) << "Ratio: " << elapsed_time / timer.current() << "\t";
      (*monitor_) << "Est. time: " << (requested_time - elapsed_time) / live_ratio << "\t";
      (*monitor_) << endl;
    }

    // Reset flags.
    simData->setNeedsRemake(false);
    _simdata_remade = false;
    _handler_remade = false;
  }

  // Possibly print a closing update
  if (print_updates && topology->getRank() == 0 && runMode == RunMode::SIM) {
    (*monitor_) << " ---- End of run ----\n\n";
  }

  // --> Post-integrate
  requested_time = 0;
  simData->post_integrate();
  integrator->post_integrate();
  for (const auto& it : additional_integrators) {
    it->post_integrate();
  }
  handler->post_integrate();
  dataMaster->post_integrate();
  forceMaster->post_integrate();
  for (const auto& m : modifiers) {
    m->post_integrate();
  }

  // Don't leave the timer hanging, even though it doesn't matter.
  timer.stop();

  // Back to idle mode
  setRunMode(RunMode::IDLE);

  // End of run barrier.
  MPIObject::barrier();
}

void GFlow::writeData(const string& dirName) {
  // Make sure data master is non-null, print a warning if any writes failed
  if (dataMaster && !dataMaster->writeToDirectory(dirName)) {
    cout << "Warning: Some writes failed.\n";
  }
}

long double GFlow::getRequestedTime() const {
  return requested_time;
}

long double GFlow::getTotalRequestedTime() const {
  return total_requested_time;
}

long double GFlow::getElapsedTime() const {
  return elapsed_time;
}

long double GFlow::getTotalTime() const {
  return total_time;
}

RealType GFlow::getBoundaryForce() const {
  return boundaryForce;
}

RealType GFlow::getBoundaryEnergy() const {
  return boundaryEnergy;
}

RealType GFlow::getDT() const {
  return integrator->getTimeStep();
}

long int GFlow::getIter() const {
  return iter;
}

int GFlow::getNumInteractions() const {
  return interactions.size();
}

int GFlow::getNumParticles() const {
  return simData->number_owned();
}

Bounds &GFlow::getBounds() {
  return simulation_bounds;
}

const Bounds &GFlow::getBounds() const {
  return simulation_bounds;
}

const BCFlag *GFlow::getBCs() const {
  return boundaryConditions;
}

BCFlag GFlow::getBC(int dim) const {
  if (dim < 0 || sim_dimensions <= dim) {
    throw BadDimension("Bad dim in get BC.");
  }
  return boundaryConditions[dim];
}

int GFlow::getNTypes() const {
  return forceMaster->getNTypes();
}

int GFlow::getSimDimensions() const {
  return sim_dimensions;
}

bool GFlow::getUseForces() const {
  return use_forces_;
}

pair<int, char **> GFlow::getCommand() const {
  return pair<int, char **>(argc, argv);
}

const vector<shared_ptr<class Interaction> > &GFlow::getInteractions() const {
  return interactions;
}

const vector<shared_ptr<class Bonded> > &GFlow::getBondedInteractions() const {
  return bondedInteractions;
}

shared_ptr<class SimData> GFlow::getSimData() {
  return simData;
}

DataMaster *GFlow::getDataMaster() {
  return dataMaster;
}

ForceMaster *GFlow::getForceMaster() {
  return forceMaster;
}

Integrator *GFlow::getIntegrator() {
  return integrator;
}

Topology *GFlow::getTopology() {
  return topology;
}

InteractionHandler *GFlow::getInteractionHandler() {
  return handler;
}

int GFlow::getNumIntegrators() const {
  return (integrator != nullptr ? 1 : 0) + additional_integrators.size();
}

void GFlow::getDisplacement(const RealType *x, const RealType *y, RealType *dis) {
  for (int d = 0; d < sim_dimensions; ++d) {
    dis[d] = x[d] - y[d];
    if (boundaryConditions[d] == BCFlag::WRAP) {
      RealType dx = simulation_bounds.max[d] - simulation_bounds.min[d] - fabs(dis[d]);
      if (dx < fabs(dis[d])) {
        dis[d] = dis[d] > 0 ? -dx : dx;
      }
    }
  }
}

void GFlow::minimumImage(RealType *dis) {
  for (int d = 0; d < sim_dimensions; ++d) {
    if (boundaryConditions[d] == BCFlag::WRAP) {
      RealType dx = simulation_bounds.max[d] - simulation_bounds.min[d] - fabs(dis[d]);
      if (dx < fabs(dis[d])) {
        dis[d] = dis[d] > 0 ? -dx : dx;
      }
    }
  }
}

void GFlow::minimumImage(RealType &dis, int d) {
  if (boundaryConditions[d] == BCFlag::WRAP) {
    RealType dx = simulation_bounds.wd(d) - fabs(dis);
    if (dx < fabs(dis)) {
      dis = dis > 0 ? -dx : dx;
    }
  }
}

RealType GFlow::getDistance(const RealType *x, const RealType *y) {
  return sqrt(getDistanceSqr(x, y));
}

RealType GFlow::getDistanceSqr(const RealType *x, const RealType *y) {
  RealType dist = 0;
  for (int d = 0; d < sim_dimensions; ++d) {
    RealType ds = x[d] - y[d];
    if (boundaryConditions[d] == BCFlag::WRAP) {
      RealType dx = simulation_bounds.max[d] - simulation_bounds.min[d] - fabs(ds);
      if (dx < fabs(ds)) {
        ds = ds > 0 ? -dx : dx;
      }
    }
    dist += sqr(ds);
  }
  return dist;
}

RunMode GFlow::getRunMode() {
  return runMode;
}

void GFlow::addInteraction(const std::shared_ptr<Interaction>& inter) {
  if (inter && !contains(interactions, inter)) {
    interactions.push_back(inter);
  }
}

void GFlow::addBonded(const std::shared_ptr<Bonded>& bnd) {
  if (bnd && !contains(bondedInteractions, bnd)) {
    bondedInteractions.push_back(bnd);
  }
}

void GFlow::addBody(const std::shared_ptr<Body>& bdy) {
  if (bdy && !contains(bodies, bdy)) {
    bodies.push_back(bdy);
  }
}

void GFlow::addIntegrator(const std::shared_ptr<Integrator>& it) {
  if (it) {
    additional_integrators.push_back(it);
  }
}

void GFlow::setCommand(int argc, char **argv) {
  if (argv) {
    this->argc = argc;
    this->argv = argv;
  }
}

void GFlow::setAllBCs(BCFlag type) {
  for (int d = 0; d < sim_dimensions; ++d) {
    boundaryConditions[d] = type;
  }
}

void GFlow::setBC(const int d, const BCFlag type) {
  boundaryConditions[d] = type;
}

void GFlow::setUseForces(bool f) {
  use_forces_ = f;
}

void GFlow::setBounds(const Bounds &bnds) {
  simulation_bounds = bnds;
  // Set topology's simulation_bounds first.
  if (topology) {
    topology->setSimulationBounds(bnds);
  }
  // Tell handler to check its simulation_bounds.
  if (handler) {
    handler->checkBounds();
  }
}

void GFlow::setRepulsion(RealType r) {
  if (r < 0) {
    return;
  }
  repulsion = r;
}

void GFlow::setDissipation(RealType d) {
  if (d < 0) {
    return;
  }
  dissipation = d;
}

void GFlow::setAttraction(RealType g) {
  center_attraction = g;
}

void GFlow::setPrintUpdates(bool p) {
  print_updates = p;
}

void GFlow::setUpdateInterval(RealType u) {
  if (u > 0) {
    update_interval_ = u;
  }
}

void GFlow::setRunMode(RunMode m) {
  runMode = m;
}

void GFlow::requestTime(RealType t) {
  if (t < 0) {
    t = 0;
  }
  requested_time = t;
}

void GFlow::setElapsedTime(RealType t) {
  elapsed_time = t;
}

void GFlow::wrapPositions() {
  // Start simdata timer.
  simData->start_timer();

  // Get a pointer to position data and the number of particles in simData
  auto x = simData->X();
  int size = simData->size_owned();
  for (int d = 0; d < sim_dimensions; ++d) {
    RealType min_bound = simulation_bounds.min[d], max_bound = simulation_bounds.max[d],
        width = simulation_bounds.wd(d);
    if (boundaryConditions[d] == BCFlag::WRAP) {
      for (int n = 0; n < size; ++n) {
        // Create a local copy
        RealType xlocal = x(n, d);
        // Periodic boundary correction.
        if (xlocal < min_bound) {
          xlocal = max_bound - fmod(min_bound - xlocal, width);
        }
        else if (max_bound <= xlocal) {
          xlocal = fmod(xlocal - min_bound, width) + min_bound;
        }
        // Set
        x(n, d) = xlocal;
      }
    }
  }

  // Start simdata timer.
  simData->stop_timer();
}

void GFlow::reflectPositions() {
  // Start simdata timer.
  simData->start_timer();

  // Get an accessor to position and velocity data and the number of particles in simData
  auto x = simData->X(), v = simData->V();
  int size = simData->size_owned();
  // Reflect all the particles
  for (int d = 0; d < sim_dimensions; ++d) {
    RealType min_bound = simulation_bounds.min[d], max_bound = simulation_bounds.max[d];
    if (boundaryConditions[d] == BCFlag::REFL) {
      for (int n = 0; n < size; ++n) {
        // Create a local copy
        RealType xlocal = x(n, d);
        if (xlocal < min_bound) {
          xlocal = 2 * min_bound - xlocal;
          v(n, d) = -v(n, d);
        }
        else if (max_bound < xlocal) {
          xlocal = 2 * max_bound - xlocal;
          v(n, d) = -v(n, d);
        }
        x(n, d) = xlocal;
      }
    }
  }

  // Start simdata timer.
  simData->stop_timer();
}

void GFlow::repulsePositions() {
  // Start simdata timer.
  simData->start_timer();

  // Get a pointer to position data and the number of particles in simData
  auto x = simData->X(), v = simData->V(), f = simData->F();
  int size = simData->size_owned();
  // Reset boundary force and energy
  boundaryForce = 0;
  boundaryEnergy = 0;
  // Reflect all the particles
  for (int d = 0; d < sim_dimensions; ++d) {
    RealType min_bound = simulation_bounds.min[d], max_bound = simulation_bounds.max[d];
    if (boundaryConditions[d] == BCFlag::REPL) {
      for (int n = 0; n < size; ++n) {
        // Create a local copy
        if (x(n, d) < min_bound) {
          RealType dx = min_bound - x(n, d);
          RealType F = 10 * repulsion * dx + dissipation * clamp(-v(n, d));
          f(n, d) += F;
          boundaryForce += F;
          boundaryEnergy += 0.5 * repulsion * sqr(dx);
        }
        else if (max_bound < x(n, d)) {
          RealType dx = (x(n, d) - max_bound);
          RealType F = 10 * repulsion * dx + dissipation * clamp(v(n, d));
          f(n, d) -= F;
          boundaryForce += F;
          boundaryEnergy += 0.5 * repulsion * sqr(dx);
        }
      }
    }
  }

  // Start simdata timer.
  simData->stop_timer();
}

void GFlow::attractPositions() {
  // Only do this if center_attraction is nonzero
  if (center_attraction == 0) {
    return;
  }

  // Start simdata timer.
  simData->start_timer();

  // Get a pointer to position data and the number of particles in simData
  auto x = simData->X(), f = simData->F();
  auto im = simData->Im();
  int size = simData->size_owned();
  // Find the center of the simulation
  auto center = new RealType[sim_dimensions];
  auto X = new RealType[sim_dimensions], dX = new RealType[sim_dimensions];
  simulation_bounds.center(center);

  // Attract particles towards center with constant acceleration
  for (int n = 0; n < size; ++n) {
    if (simData->Type(n) < 0) {
      continue;
    }
    copyVec(x(n), X, sim_dimensions);
    subtractVec(center, X, dX, sim_dimensions);
    normalizeVec(dX, sim_dimensions);
    scalarMultVec(center_attraction / im(n), dX, sim_dimensions);
    plusEqVec(f(n), dX, sim_dimensions);
  }
  // Clean up
  delete[] center;
  delete[] X;

  // Start simdata timer.
  simData->stop_timer();
}

void GFlow::removeOverlapping(RealType fraction) {
  if (handler) {
    handler->removeOverlapping(fraction);
  }
}

void GFlow::addDataObject(const shared_ptr<DataObject>& dob) {
  if (dob) {
    dataMaster->addDataObject(dob);
    // Check if this object inherits from group.
    auto group = std::dynamic_pointer_cast<Group>(dob);
    if (group != nullptr) {
      global_id_reliant.push_back(group);
    }
  }
}

void GFlow::addModifier(const shared_ptr<Modifier>& mod) {
  if (mod) {
    modifiers.push_back(mod);
    // Check if this object inherits from group.
    auto group = std::dynamic_pointer_cast<Group>(mod);
    if (group != nullptr) {
      global_id_reliant.push_back(group);
    }
  }
}

void GFlow::resetAllTimes() {
  requested_time = 0.;
  total_requested_time = 0.;
  elapsed_time = 0.;
  total_time = 0.;
  iter = 0;
}

void GFlow::setStartRecTime(RealType t) {
  if (dataMaster) {
    dataMaster->setStartRecTime(t);
  }
}

void GFlow::setFPS(RealType fps) {
  if (dataMaster) {
    dataMaster->setFPS(fps);
  }
}

void GFlow::setFPS(int dob_id, RealType fps) {
  if (dataMaster) {
    dataMaster->setFPS(dob_id, fps);
  }
}

void GFlow::setDT(RealType dt) {
  if (integrator) {
    integrator->setDT(dt);
  }
}

void GFlow::setMaxDT(RealType mdt) {
  if (integrator) {
    integrator->setMaxDT(mdt);
  }
}

void GFlow::setDMCmd(int argc, char **argv) {
  if (dataMaster) {
    dataMaster->setCommand(argc, argv);
  }
}

void GFlow::giveFileToDataMaster(const string& filename, const string& file_contents) {
  if (dataMaster) {
    dataMaster->giveFile(filename, file_contents);
  }
}

void GFlow::startMPIExchangeTimer() {
  mpi_exchange_timer.start_timer();
}

void GFlow::stopMPIExchangeTimer() {
  mpi_exchange_timer.stop_timer();
}

void GFlow::startMPIGhostTimer() {
  mpi_ghost_timer.start_timer();
}

void GFlow::stopMPIGhostTimer() {
  mpi_ghost_timer.stop_timer();
}

void GFlow::syncRunning() {
  // Sync amongst processors.
  MPIObject::mpi_or(terminate_);
  // Check if any processors called for program termination.
  if (terminate_) {
    running_ = false;
  }
}

void GFlow::terminate() {
  // If this is the only processor, we don't have to wait for anyone.
  if (topology->getNumProc() == 1) {
    running_ = false;
  }
  // By doing this, we only record the first termination time.
  if (!terminate_) {
    // Set terminate flag to true.
    terminate_ = true;
    // Set the terminate time
    termination_time = elapsed_time;
  }
}

void GFlow::registerException(Exception *ex) {
  error_handling.push_back(ex);
  terminate();
}

inline void GFlow::clearForces() {
  simData->clearF();
  // Clear torque buffer, if it exists.
  if (sim_dimensions == 2) {
    simData->clearScalar("Tq");
  }
}

inline void GFlow::handleModifiers() {
  // If there are no modifiers, there is nothing to do!
  if (modifiers.empty()) {
    return;
  }
  // List of modifiers to remove
  vector<decltype(modifiers)::iterator> remove;
  // Find modifiers that need to be removed
  for (auto it = modifiers.begin(); it != modifiers.end(); ++it) {
    if ((*it)->getRemove()) {
      remove.push_back(it);
    }
  }
  // Remove modifiers
  for (auto &m : remove) {
    modifiers.erase(m);
  }
}
