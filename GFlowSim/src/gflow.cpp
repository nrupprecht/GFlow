#include "gflow.hpp"
// Other files
#include "allbaseobjects.hpp"
#include "allbodies.hpp"
#include "alldomains.hpp"
#include "allparallelobjects.hpp"
#include "allmodifiers.hpp"
#include "alltopologies.hpp"
#include "allbonded.hpp"
#include "allbodies.hpp"

namespace GFlowSimulation {

  GFlow::GFlow(int dims) : requested_time(0), total_requested_time(0), elapsed_time(0), total_time(0), 
    iter(0), repulsion(DEFAULT_HARD_SPHERE_REPULSION), dissipation(0), center_attraction(0), sim_dimensions(dims),
    bounds(dims)
    
  {
    // Set up basic objects. The integrator will be created by the creator.
    simData = std::make_shared<SimData>(this);
    integrator   = nullptr;
    handler      = new Domain(this);
    dataMaster   = new DataMaster(this);
    forceMaster  = new ForceMaster(this);
    topology     = new KDTreeTopology(this);
    // Set up boundary conditions
    boundaryConditions = new BCFlag[sim_dimensions];
    // Set wrapping to true by default
    setAllBCs(BCFlag::WRAP);
    // Initialize all base objects, so they have valid pointers to one another
    initialize();
  }

  GFlow::~GFlow() {
    if (integrator)   delete integrator;
    if (handler)      delete handler;
    if (dataMaster)   delete dataMaster;
    if (forceMaster)  delete forceMaster;
    if (topology)     delete topology;
    integrator = nullptr;
    handler = nullptr;
    dataMaster = nullptr;
    forceMaster = nullptr;
    topology = nullptr;
    if (boundaryConditions) delete [] boundaryConditions;
    boundaryConditions = nullptr;
  }

  bool GFlow::initialize() {
    bool non_null = true;
    // --- Initialize all the subobjects
    simData->initialize();

    if (topology) topology->initialize();
    else non_null = false;

    if(integrator) integrator->initialize();
    else non_null = false;

    if (forceMaster) forceMaster->initialize(); // Should go before handler
    else non_null = false;

    if (handler) handler->initialize();
    else non_null = false;

    if (dataMaster) dataMaster->initialize();
    else non_null = false;

    for (auto &md : modifiers) {
      if (md) md->initialize();
      else non_null = false; 
    }

    for (auto &it: interactions) {
      if (it) it->initialize();
      else non_null = false;
    }
    // Return whether pointers were non-null
    return non_null;
  }

  void GFlow::initializeBase(Base *base) {
    // Make sure we aren't handed a null pointer
    if (base==nullptr) return;
    // Set the number of dimensions
    base->sim_dimensions = sim_dimensions;
    // Give pointer to this GFlow object
    base->gflow = this;
    // Set other objects
    base->simData      = simData;
    base->integrator   = integrator;
    base->handler      = handler;
    base->dataMaster   = dataMaster;
    base->forceMaster  = forceMaster;
    base->topology     = topology;
    // Set vectors
    // base->modifiersPtr = &modifiers;
    // base->interactionsPtr = &interactions;
  }

  void GFlow::run(long double rt) {
    // If a parameter was passed in, it is the requested time
    if (rt>0) requested_time = rt;

    // Record this request
    total_requested_time += requested_time;

    // If there are no particles, we are done
    if (simData->number()==0) {
      elapsed_time += requested_time;
      total_time   += requested_time;
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
    if (runMode==RunMode::IDLE) setRunMode(RunMode::SIM);

    // GFlow local timers. Timed objects' timers are cleared by call to pre_integrate.
    bonded_timer.clear_timer();
    body_timer.clear_timer();
    mpi_exchange_timer.clear_timer();
    mpi_ghost_timer.clear_timer();

    // --> Pre-integrate
    _running = true;
    _terminate = false;
    elapsed_time = 0;
    iter = 0;
    // Set up all objects.
    // Put modifiers first so if anything is modified, the other objects can react to it.
    for (auto m : modifiers) m->pre_integrate(); 
    for (auto b : bodies) b->pre_integrate();
    simData->pre_integrate();
    integrator->pre_integrate();
    for (auto it : additional_integrators) it->pre_integrate();
    handler->pre_integrate();
    dataMaster->pre_integrate();
    forceMaster->pre_integrate();

    // Run Timer - for printing updates. We use a timer, not a timed object, so that if timed objects are turned off,
    // we still can collect timing data.
    Timer timer;
    timer.start();

    // Do integration for the requested amount of time
    while (_running && requested_time>0) {

      // --> Pre-step
      for (auto m : modifiers) m->pre_step();
      integrator->pre_step();
      for (auto it : additional_integrators) it->pre_step();
      dataMaster->pre_step();
      handler->pre_step();

      // --> Pre-force
      for (auto m : modifiers) m->pre_forces();
      integrator->pre_forces(); // -- This is where VV first half kick happens (if applicable)
      for (auto it : additional_integrators) it->pre_forces(); // -- First half kick could also happen here.
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
      if (!_handler_remade && _use_ghosts) {
        topology->send_ghost_updates();
        topology->recv_ghost_updates();
      }

      // Calculate interactions and forces.
      if (_use_forces) {
        // Calculate interactions.
        forceMaster->interact();
        forceMaster->interact_ghosts();

        // Calculate bonded interactions.
        if (!bondedInteractions.empty()) {
          bonded_timer.start_timer();
          for (auto bd : bondedInteractions) bd->interact();
          bonded_timer.stop_timer();
        }
      }

      /*
      // Finish updating ghost particles.
      if (!_handler_remade && _use_ghosts) simData->finishGhostParticleUpdates();      
      // Calculate ghost particle forces
      if (_use_forces) forceMaster->interact_ghosts();
      */
      
      // Update bodies - this may include forces.
      if (!bodies.empty()) {
        body_timer.start_timer();
        for (auto b : bodies) b->correct();
        body_timer.stop_timer();
      }

      // Do modifier removal
      handleModifiers();

      // --> Post-forces
      for (auto m : modifiers) m->post_forces(); // -- This is where modifiers should do forces (if they need to)
      // Update halo particles. This should be done after all force calculations, but before the next integration steps.
      // For this reason, this occurs *after* modifiers do their post-force routine.
      simData->updateHaloParticles();
      // Continue with normal order of updates.
      simData->post_forces();
      integrator->post_forces();                 // -- This is where VV second half kick happens (if applicable)
      for (auto it : additional_integrators) it->post_forces();
      dataMaster->post_forces();
      handler->post_forces();

      // --> Post-step
      if (requested_time<=elapsed_time) terminate();
      for (auto m : modifiers) m->post_step();
      integrator->post_step();
      for (auto it : additional_integrators) it->post_step();
      dataMaster->post_step();
      handler->post_step();

      // Timer updates
      ++iter;
      RealType dt = integrator->getTimeStep();
      elapsed_time += dt;
      total_time += dt;
      // Check for bad numerical precision
      if (total_time - dt == total_time) {
        cout << "Loss of precision. Stopping simulation.\n";
	      terminate();
      }
      // Possibly print updates to the screen or to a file.
      if (print_updates && topology->getRank()==0 && runMode==RunMode::SIM && 
        static_cast<int>((elapsed_time-dt)/update_interval) < static_cast<int>((elapsed_time)/update_interval)) 
      {
        RealType live_ratio = elapsed_time / timer.current();
        (*monitor) << "Simulation time: " << static_cast<int>(elapsed_time) << "\t";
        (*monitor) << "Ratio: " << elapsed_time / timer.current() << "\t";
        (*monitor) << "Est. time: " << (requested_time - elapsed_time)/live_ratio << "\t";
	      (*monitor) << endl;
      }

      // Reset flags.
      simData->setNeedsRemake(false);
      _simdata_remade = false;
      _handler_remade = false;
    }
    
    // Possibly print a closing update
    if (print_updates && topology->getRank()==0 && runMode==RunMode::SIM) {
      (*monitor) << " ---- End of run ----\n\n";
    }

    // --> Post-integrate
    requested_time = 0;
    simData->post_integrate();
    integrator->post_integrate();
    for (auto it : additional_integrators) it->post_integrate();
    handler->post_integrate();
    dataMaster->post_integrate();
    forceMaster->post_integrate();
    for (auto m : modifiers) m->post_integrate();

    // Don't leave the timer hanging, even though it doesn't matter.
    timer.stop();

    // Back to idle mode
    setRunMode(RunMode::IDLE);

    // End of run barrier.
    MPIObject::barrier();
  }

  void GFlow::writeData(string dirName) {
    // Make sure data master is non-null, print a warning if any writes failed
    if (dataMaster && !dataMaster->writeToDirectory(dirName))
      cout << "Warning: Some writes failed.\n";
  }

  long double GFlow::getRequestedTime() const {
    return requested_time;
  }

  long double GFlow::getTotalRequestedTime() const {
    return total_requested_time;
  }

  long double GFlow::getElapsedTime() const{
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

  Bounds GFlow::getBounds() const {
    return bounds;
  }

  const BCFlag* GFlow::getBCs() const {
    return boundaryConditions;
  }

  BCFlag GFlow::getBC(int dim) const {
    if (dim<0 || sim_dimensions<=dim) 
      throw BadDimension("Bad dim in get BC.");
    return boundaryConditions[dim];
  }

  int GFlow::getNTypes() const {
    return forceMaster->getNTypes();
  }

  int GFlow::getSimDimensions() const {
    return sim_dimensions;
  }

   bool GFlow::getUseForces() const {
    return _use_forces;
  }

  pair<int, char**> GFlow::getCommand() const {
    return pair<int, char**>(argc, argv);
  }

  const vector<shared_ptr<class Interaction> >& GFlow::getInteractions() const {
    return interactions;
  }

  const vector<shared_ptr<class Bonded> >& GFlow::getBondedInteractions() const {
    return bondedInteractions;
  }

  shared_ptr<class SimData> GFlow::getSimData() {
    return simData;
  }

  DataMaster* GFlow::getDataMaster() {
    return dataMaster;
  }

  ForceMaster* GFlow::getForceMaster() {
    return forceMaster;
  }

  Integrator* GFlow::getIntegrator() {
    return integrator;
  }

  Topology* GFlow::getTopology() {
    return topology;
  }

  InteractionHandler* GFlow::getInteractionHandler() {
    return handler;
  }

  int GFlow::getNumIntegrators() const {
    return (integrator!=nullptr ? 1 : 0) + additional_integrators.size();
  }

  void GFlow::getDisplacement(const RealType *x, const RealType *y, RealType *dis) {
    for (int d=0; d<sim_dimensions; ++d) {
      dis[d] = x[d] - y[d];
      if (boundaryConditions[d]==BCFlag::WRAP) {
        RealType dx = bounds.max[d] - bounds.min[d] - fabs(dis[d]);
        if (dx<fabs(dis[d])) dis[d] = dis[d]>0 ? -dx : dx;
      }      
    }
  }

  void GFlow::minimumImage(RealType *dis) {
    for (int d=0; d<sim_dimensions; ++d) {
      if (boundaryConditions[d]==BCFlag::WRAP) {
        RealType dx = bounds.max[d] - bounds.min[d] - fabs(dis[d]);
        if (dx<fabs(dis[d])) dis[d] = dis[d]>0 ? -dx : dx;
      }      
    }
  }

  void GFlow::minimumImage(RealType &dis, int d) {
    if (boundaryConditions[d]==BCFlag::WRAP) {
      RealType dx = bounds.wd(d) - fabs(dis);
      if (dx<fabs(dis)) dis = dis>0 ? -dx : dx;
    }  
  }

  RealType GFlow::getDistance(const RealType *x, const RealType *y) {
    return sqrt(getDistanceSqr(x, y));
  }

  RealType GFlow::getDistanceSqr(const RealType *x, const RealType *y) {
    RealType dist = 0;
    for (int d=0; d<sim_dimensions; ++d) {
      RealType ds = x[d] - y[d];
      if (boundaryConditions[d]==BCFlag::WRAP) {
        RealType dx = bounds.max[d] - bounds.min[d] - fabs(ds);
        if (dx<fabs(ds)) ds = ds>0 ? -dx : dx;
      }
      dist += sqr(ds);
    }
    return dist;
  }

  RunMode GFlow::getRunMode() {
    return runMode;
  }

  void GFlow::addInteraction(shared_ptr<Interaction> inter) {
    if (inter && !contains(interactions, inter))
      interactions.push_back(inter);
  }

  void GFlow::addBonded(shared_ptr<Bonded> bnd) {
    if (bnd && !contains(bondedInteractions, bnd))
      bondedInteractions.push_back(bnd);
  }

  void GFlow::addBody(shared_ptr<Body> bdy) {
    if (bdy && !contains(bodies, bdy))
      bodies.push_back(bdy);
  }

  void GFlow::addIntegrator(shared_ptr<Integrator> it) {
    if (it)
      additional_integrators.push_back(it);
  }

  void GFlow::setCommand(int argc, char **argv) {
    if (argv) {
      this->argc = argc;
      this->argv = argv;
    }
  }

  void GFlow::setAllBCs(BCFlag type) {
    for (int d=0; d<sim_dimensions; ++d) 
      boundaryConditions[d] = type;
  }

  void GFlow::setBC(const int d, const BCFlag type) {
    boundaryConditions[d] = type;
  }

  void GFlow::setUseForces(bool f) {
    _use_forces = f;
  }

  void GFlow::setBounds(const Bounds& bnds) {
    bounds = bnds;
    // Set topology's bounds first.
    if (topology) topology->setSimulationBounds(bnds);
    // Tell handler to check its bounds.
    if (handler) handler->checkBounds();
  }

  void GFlow::setRepulsion(RealType r) {
    if (r<0) return;
    repulsion = r;
  }

  void GFlow::setDissipation(RealType d) {
    if (d<0) return;
    dissipation = d;
  }

  void GFlow::setAttraction(RealType g) {
    center_attraction = g;
  }

  void GFlow::setPrintUpdates(bool p) {
    print_updates = p;
  }

  void GFlow::setUpdateInterval(RealType u) {
    if (u>0) update_interval = u;
  }

  void GFlow::setRunMode(RunMode m) {
    runMode = m;
  }

  void GFlow::requestTime(RealType t) {
    if (t<0) t = 0;
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
    for (int d=0; d<sim_dimensions; ++d) {
      RealType min_bound = bounds.min[d], max_bound = bounds.max[d], width = bounds.wd(d);
      if (boundaryConditions[d]==BCFlag::WRAP) {
        for (int n=0; n<size; ++n) {
          // Create a local copy
          RealType xlocal = x(n, d);
          // Periodic boundary correction.
          if (xlocal<min_bound) xlocal = max_bound-fmod(min_bound-xlocal, width);
          else if (max_bound<=xlocal) xlocal = fmod(xlocal-min_bound, width)+min_bound;
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
    for (int d=0; d<sim_dimensions; ++d) {
      RealType min_bound = bounds.min[d], max_bound = bounds.max[d];
      if (boundaryConditions[d]==BCFlag::REFL) { 
        for (int n=0; n<size; ++n) {
          // Create a local copy
          RealType xlocal = x(n, d);
          if (xlocal<min_bound) {
            xlocal = 2*min_bound - xlocal ;
            v(n, d) = -v(n, d);
          }
          else if (max_bound<xlocal) {
            xlocal = 2*max_bound - xlocal;
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
    for (int d=0; d<sim_dimensions; ++d) {
      RealType min_bound = bounds.min[d], max_bound = bounds.max[d];
      if (boundaryConditions[d]==BCFlag::REPL) { 
        for (int n=0; n<size; ++n) {
          // Create a local copy
          if (x(n, d)<min_bound) {
            RealType dx = min_bound - x(n, d);
            RealType F = 10*repulsion*dx + dissipation*clamp(-v(n, d));
            f(n, d) += F;
            boundaryForce += F;
            boundaryEnergy += 0.5*repulsion*sqr(dx);
          }
          else if (max_bound<x(n, d)) {
            RealType dx = (x(n, d) - max_bound);
            RealType F = 10*repulsion*dx + dissipation*clamp(v(n, d));
            f(n, d) -= F;
            boundaryForce += F;
            boundaryEnergy += 0.5*repulsion*sqr(dx);
          }
        }
      }
    }

    // Start simdata timer.
    simData->stop_timer();
  }

  void GFlow::attractPositions() {
    // Start simdata timer.
    simData->start_timer();

    // Only do this if center_attraction is nonzero
    if (center_attraction==0) return;
    // Get a pointer to position data and the number of particles in simData
    auto x = simData->X(), f = simData->F();
    auto im = simData->Im();
    int size = simData->size_owned();
    // Find the center of the simulation
    RealType *center = new RealType[sim_dimensions];
    RealType *X = new RealType[sim_dimensions], *dX = new RealType[sim_dimensions];
    bounds.center(center);

    // Attract particles towards center with constant acceleration
    for (int n=0; n<size; ++n) {
      if (simData->Type(n)<0) continue;
      copyVec(x(n), X, sim_dimensions);
      subtractVec(center, X, dX, sim_dimensions);
      normalizeVec(dX, sim_dimensions);
      scalarMultVec(center_attraction/im(n), dX, sim_dimensions);
      plusEqVec(f(n), dX, sim_dimensions);
    }
    // Clean up
    delete [] center;
    delete [] X;

    // Start simdata timer.
    simData->stop_timer();
  }

  void GFlow::removeOverlapping(RealType fraction) {
    if (handler) handler->removeOverlapping(fraction);
  }

  void GFlow::addDataObject(shared_ptr<DataObject> dob) {
    if (dob) dataMaster->addDataObject(dob);
  }

  void GFlow::addModifier(shared_ptr<Modifier> mod) {
    if (mod) modifiers.push_back(mod);
  }

  void GFlow::resetAllTimes() {
    requested_time       = 0.;
    total_requested_time = 0.;
    elapsed_time         = 0.;
    total_time           = 0.;
    iter                 = 0 ;
  }

  void GFlow::setStartRecTime(RealType t) {
    if (dataMaster) dataMaster->setStartRecTime(t);
  }

  void GFlow::setFPS(RealType fps) {
    if (dataMaster) dataMaster->setFPS(fps);
  }

  void GFlow::setFPS(int dob_id, RealType fps) {
    if (dataMaster) dataMaster->setFPS(dob_id, fps);
  }

  void GFlow::setDT(RealType dt) {
    if (integrator) integrator->setDT(dt);
  }

  void GFlow::setMaxDT(RealType mdt) {
    if (integrator) integrator->setMaxDT(mdt);
  }

  void GFlow::setDMCmd(int argc, char** argv) {
    if (dataMaster) dataMaster->setCommand(argc, argv);
  }

  void GFlow::giveFileToDataMaster(string filename, string file_contents) {
    if (dataMaster) dataMaster->giveFile(filename, file_contents);
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
    MPIObject::mpi_or(_terminate);
    // Check if any processors called for program termination.
    if (_terminate) _running = false;
  }

  void GFlow::terminate() {
    // If this is the only processor, we don't have to wait for anyone.
    if (topology->getNumProc()==1) _running = false;
    // By doing this, we only record the first termination time.
    if (!_terminate) {
      // Set terminate flag to true.
      _terminate = true;
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
    if (sim_dimensions==2) simData->clearScalar("Tq");
  }

  inline void GFlow::handleModifiers() {
    // If there are no modifiers, there is nothing to do!
    if (modifiers.empty()) return;
    // List of modifiers to remove
    vector< decltype(modifiers)::iterator > remove;
    // Find modifiers that need to be removed
    for (auto it = modifiers.begin(); it!=modifiers.end(); ++it) {
      if ((*it)->getRemove()) remove.push_back(it);
    }
    // Remove modifiers 
    for (auto &m : remove) modifiers.erase(m);
  }

}
