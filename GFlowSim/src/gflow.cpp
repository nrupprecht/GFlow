#include "gflow.hpp"
// Other files
#include "allbaseobjects.hpp"
#include "allbodies.hpp"
#include "alldomains.hpp"
#include "allparallelobjects.hpp"
#include "allmodifiers.hpp"
#include "alltopologies.hpp"
#include "allbonded.hpp"
#include "base/dataobject.hpp"

namespace GFlowSimulation {

  GFlow::GFlow(int dims) : running(false), useForces(true), requested_time(0), total_requested_time(0), elapsed_time(0), total_time(0), 
    iter(0), argc(0), argv(nullptr), repulsion(DEFAULT_HARD_SPHERE_REPULSION), dissipation(0), center_attraction(0), sim_dimensions(dims),
    bounds(Bounds(2))
  {
    #if USE_MPI == 1
    #if _CLANG_ == 1
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    #else
    rank    = MPI::COMM_WORLD.Get_rank();
    numProc = MPI::COMM_WORLD.Get_size();
    #endif
    #endif
    // Set the static counter.
    //DataObject::resetTotalObjects();
    // Set up basic objects. The integrator will be created by the creator.
    simData      = new SimData(this);
    integrator   = nullptr;
    domain       = new Domain(this);
    dataMaster   = new DataMaster(this);
    forceMaster  = new ForceMaster(this);
    topology     = new GridTopology(sim_dimensions);
    // Set up boundary conditions
    boundaryConditions = new BCFlag[sim_dimensions];
    // Set up bounds to have the propper dimensions
    bounds = Bounds(sim_dimensions);
    // Set wrapping to true by default
    setAllBCs(BCFlag::WRAP);
    // Initialize all base objects, so they have valid pointers to one another
    initialize();
  }

  GFlow::~GFlow() {
    if (simData)      delete simData;
    if (integrator)   delete integrator;
    if (domain)       delete domain;
    if (dataMaster)   delete dataMaster;
    if (forceMaster)  delete forceMaster;
    if (topology)     delete topology;
    simData = nullptr;
    integrator = nullptr;
    domain = nullptr;
    dataMaster = nullptr;
    forceMaster = nullptr;
    topology = nullptr;
    for (auto &md : modifiers) {
      if (md) delete md;
      md = nullptr;
    }
    for (auto &it : interactions) {
      if (it) delete it;
      it = nullptr;
    }
    if (boundaryConditions) delete [] boundaryConditions;
    boundaryConditions = nullptr;
  }

  bool GFlow::initialize() {
    bool non_null = true;
    // --- Initialize all the subobjects
    if (simData) simData->initialize();
    else non_null = false;

    if(integrator) integrator->initialize();
    else non_null = false;

    if (forceMaster) forceMaster->initialize(); // Should go before domain
    else non_null = false;

    if (domain) domain->initialize();
    else non_null = false;

    if (dataMaster) dataMaster->initialize();
    else non_null = false;

    // Topology does not currently have an initialization function, as it does not inherit from Base

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
    base->domain       = domain;
    base->dataMaster   = dataMaster;
    base->forceMaster  = forceMaster;
    base->topology     = topology;
    // Set vectors
    base->modifiersPtr = &modifiers;
    base->interactionsPtr = &interactions;
    // Set MPI
    base->rank = rank;
    base->numProc = numProc;
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

    // Clear timers
    integrator->clear_timer();
    domain_timer.clear();
    forceMaster->clear_timer();
    bonded_timer.clear();
    body_timer.clear();

    // Make sure we have initialized everything
    if (!initialize()) {
      // Some object was null
      throw UnexpectedNullPointer("Error: Some object was null at GFlow initialization.");
    }
    // Check that simdata has good arrays. 
    if (simData->X()==nullptr || simData->V()==nullptr || simData->F()==nullptr) 
    {
      throw UnexpectedNullPointer("Some array in simdata was null that shouldn't be.");
    }

    // Set run mode, if it is currently "idle"
    if (runMode==RunMode::IDLE) setRunMode(RunMode::SIM);

    // --> Pre-integrate
    running = true;
    elapsed_time = 0;
    iter = 0;
    // Set up all objects
    for (auto m : modifiers) m->pre_integrate();
    for (auto b : bodies) b->pre_integrate();
    simData->pre_integrate();
    integrator->pre_integrate();
    domain->pre_integrate();
    dataMaster->pre_integrate();
    forceMaster->pre_integrate();

    Timer timer;
    timer.start();

    // Do integration for the requested amount of time
    while (running && requested_time>0) {
      // --> Pre-step
      for (auto m : modifiers) m->pre_step();
      integrator->pre_step();
      dataMaster->pre_step();
      domain->pre_step();

      // --> Pre-exchange
      for (auto m : modifiers) m->pre_exchange();
      integrator->pre_exchange();
      dataMaster->pre_exchange();
      domain->pre_exchange();

      // --> Pre-force
      for (auto m : modifiers) m->pre_forces();
      integrator->pre_forces(); // -- This is where VV first half kick happens (if applicable)
      dataMaster->pre_forces();
      domain_timer.start();
      if (useForces) domain->pre_forces();   // -- This is where resectorization / verlet list creation might happen
      domain_timer.stop();
      
      // --- Do interactions
      clearForces(); // Clear force buffers

      // Reflect or repulse particles. We only need to wrap before sectorizing particles, but we need to apply forces at every timestep.
      reflectPositions(); // This only involves velocities, so it could be done before or after clear forces.
      repulsePositions(); // But this needs to be done after clear forces.
      attractPositions(); // This does too.

      // Calculate interactions and forces.
      if (useForces) {
        // Calculate interactions
        forceMaster->interact();

        // Calculate bonded interactions
        if (!bondedInteractions.empty()) {
          bonded_timer.start();
          for (auto bd : bondedInteractions) bd->interact();
          bonded_timer.stop();
        }
      }
      // Update bodies
      if (!bodies.empty()) {
        body_timer.start();
        for (auto b : bodies) b->correct();
        body_timer.stop();
      }

      // Update halo particles
      simData->updateHaloParticles();

      // Do modifier removal
      handleModifiers();

      // --> Post-forces
      for (auto m : modifiers) m->post_forces(); // -- This is where modifiers should do forces (if they need to)
      integrator->post_forces();                 // -- This is where VV second half kick happens (if applicable)
      dataMaster->post_forces();
      domain->post_forces();

      // --> Post-step
      if (requested_time<=elapsed_time) running = false;
      for (auto m : modifiers) m->post_step();
      integrator->post_step();
      dataMaster->post_step();
      domain->post_step();
      // Timer updates
      ++iter;
      RealType dt = integrator->getTimeStep();
      elapsed_time += dt;
      total_time += dt;

      // Check for bad numerical precision
      if (total_time - dt == total_time) {
        cout << "Loss of precision. Stopping simulation.\n";
	      running = false;
      }
      // Possibly print updates to the screen or to a file.
      if (print_updates && runMode==RunMode::SIM && 
        static_cast<int>((elapsed_time-dt)/update_interval) < static_cast<int>((elapsed_time)/update_interval)) 
      {
        RealType live_ratio = elapsed_time / timer.current();
        (*monitor) << "Simulation time: " << static_cast<int>(elapsed_time) << "\t";
        (*monitor) << "Ratio: " << elapsed_time / timer.current() << "\t";
        (*monitor) << "Est. time: " << (requested_time - elapsed_time)/live_ratio << "\t";
	      (*monitor) << endl;
      }
      
      // Reset simdata needs remake flag
      simData->setNeedsRemake(false);
    }

    // Possibly print a closing update
    if (print_updates && runMode==RunMode::SIM) {
      (*monitor) << " ---- End of run ----\n\n";
    }

    // --> Post-integrate
    requested_time = 0;
    simData->post_integrate();
    integrator->post_integrate();
    domain->post_integrate();
    dataMaster->post_integrate();
    forceMaster->post_integrate();
    for (auto m : modifiers) m->post_integrate();

    // Back to idle mode
    setRunMode(RunMode::IDLE);
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
    return simData->number();
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

  pair<int, char**> GFlow::getCommand() const {
    return pair<int, char**>(argc, argv);
  }

  const vector<class Interaction*>& GFlow::getInteractions() const {
    return interactions;
  }

  SimData* GFlow::getSimData() {
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
    RealType dist = 0;
    for (int d=0; d<sim_dimensions; ++d) {
      RealType ds = x[d] - y[d];
      if (boundaryConditions[d]==BCFlag::WRAP) {
        RealType dx = bounds.max[d] - bounds.min[d] - fabs(ds);
        if (dx<fabs(ds)) ds = ds>0 ? -dx : dx;
      }
      dist += sqr(ds);
    }
    return sqrt(dist);
  }

  RunMode GFlow::getRunMode() {
    return runMode;
  }

  void GFlow::addInteraction(Interaction *inter) {
    if (inter!=nullptr && !contains(interactions, inter))
      interactions.push_back(inter);
  }

  void GFlow::addBonded(Bonded *bnd) {
    if (bnd!=nullptr && !contains(bondedInteractions, bnd))
      bondedInteractions.push_back(bnd);
  }

  void GFlow::addBody(class Body *bdy) {
    if (bdy!=nullptr && !contains(bodies, bdy))
      bodies.push_back(bdy);
  }

  void GFlow::setCommand(int argc, char **argv) {
    if (argv) {
      this->argc = argc;
      this->argv = argv;
    }
  }

  void GFlow::setBounds(Bounds b) {
    bounds = b;

    topology->setSimulationBounds(b);
    // Sectorization updates its bounds in pre_integrate()
  }

  void GFlow::setAllBCs(BCFlag type) {
    for (int d=0; d<sim_dimensions; ++d) 
      boundaryConditions[d] = type;
  }

  void GFlow::setBC(const int d, const BCFlag type) {
    boundaryConditions[d] = type;
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

  void GFlow::setRunning(bool r) {
    running = r;
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
    int size = simData->size();
    for (int n=0; n<size; ++n) 
      for (int d=0; d<sim_dimensions; ++d) {
        if (boundaryConditions[d]==BCFlag::WRAP) {
          // Create a local copy
          RealType xlocal = simData->X(n, d); //x[n][d];
          // Wrap xlocal
          if (xlocal<bounds.min[d])
            xlocal = bounds.max[d]-fmod(bounds.min[d]-xlocal, bounds.wd(d));
          else if (bounds.max[d]<=xlocal)
            xlocal = fmod(xlocal-bounds.min[d], bounds.wd(d))+bounds.min[d];
          // Set
          simData->X(n, d) = xlocal; // x[n][d] = xlocal;
        }
      }

    /*
    // Wrap all particles
    for (int d=0; d<sim_dimensions; ++d) {
      if (boundaryConditions[d]==BCFlag::WRAP) { // Wrap the d-th dimension
        for (int n=0; n<size; ++n) {
          // Create a local copy
          RealType xlocal = x[n][d];
          // Wrap xlocal
          if (xlocal<bounds.min[d])
            xlocal = bounds.max[d]-fmod(bounds.min[d]-xlocal, bounds.wd(d));
          else if (bounds.max[d]<=x[n][d])
            xlocal = fmod(xlocal-bounds.min[d], bounds.wd(d))+bounds.min[d];
          // Set
          x[n][d] = xlocal;
        }
      }
    }
    */
    
  }

  void GFlow::reflectPositions() {
    // Get a pointer to position data and the number of particles in simData
    RealType **x = simData->X(), **v = simData->V();
    int size = simData->size();

    // Reflect all the particles
    for (int d=0; d<sim_dimensions; ++d)
      if (boundaryConditions[d]==BCFlag::REFL) { 
        for (int n=0; n<size; ++n) {
          // Create a local copy
          RealType xlocal = x[n][d];
          if (xlocal<bounds.min[d]) {
            xlocal = 2*bounds.min[d] - xlocal ;
            v[n][d] = -v[n][d];
          }
          else if (bounds.max[d]<xlocal) {
            xlocal = 2*bounds.max[d] - xlocal;
            v[n][d] = -v[n][d];
          }
          x[n][d] = xlocal;
        }
      }
  }

  void GFlow::repulsePositions() {
    // Get a pointer to position data and the number of particles in simData
    RealType **x = simData->X(), **v = simData->V(), **f = simData->F();
    int size = simData->size();
    // Reset boundary force
    boundaryForce = 0;
    // Reflect all the particles
    for (int d=0; d<sim_dimensions; ++d)
      if (boundaryConditions[d]==BCFlag::REPL) { 
        for (int n=0; n<size; ++n) {
          // Create a local copy
          if (x[n][d]<bounds.min[d]) {
            RealType F = (repulsion*(bounds.min[d] - x[n][d]) + dissipation*clamp(-v[n][d]));
            f[n][d] += F;
            boundaryForce += F;
          }
          else if (bounds.max[d]<x[n][d]) {
            RealType F = repulsion*(x[n][d] - bounds.max[d]) + dissipation*clamp(v[n][d]);
            f[n][d] -= F;
            boundaryForce += F;
          }
        }
      }
  }

  void GFlow::attractPositions() {
    // Only do this if center_attraction is nonzero
    if (center_attraction==0) return;
    // Get a pointer to position data and the number of particles in simData
    RealType **x = simData->X(), **f = simData->F();
    RealType *im = simData->Im();
    int size = simData->size();
    // Find the center of the simulation
    RealType *center = new RealType[sim_dimensions];
    RealType *X = new RealType[sim_dimensions], *dX = new RealType[sim_dimensions];
    bounds.center(center);

    // Attract particles towards center with constant acceleration
    for (int n=0; n<size; ++n) {
      if (simData->Type(n)<0) continue;
      copyVec(x[n], X, sim_dimensions);
      subtractVec(center, X, dX, sim_dimensions);
      normalizeVec(dX, sim_dimensions);
      scalarMultVec(center_attraction/im[n], dX, sim_dimensions);
      plusEqVec(f[n], dX, sim_dimensions);
    }
    // Clean up
    delete [] center;
    delete [] X;
  }

  void GFlow::removeOverlapping(RealType fraction) {
    if (domain) domain->removeOverlapping(fraction);
  }

  void GFlow::addDataObject(class DataObject* dob) {
    if (dob) dataMaster->addDataObject(dob);
  }

  void GFlow::addModifier(class Modifier* mod) {
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

  inline void GFlow::clearForces() {
    simData->clearF();
  }

  inline void GFlow::handleModifiers() {
    // If there are no modifiers, there is nothing to do!
    if (modifiers.empty()) return;
    // List of modifiers to remove
    vector< std::list<Modifier*>::iterator > remove;
    // Find modifiers that need to be removed
    for (auto it = modifiers.begin(); it!=modifiers.end(); ++it) {
      if ((*it)->getRemove()) remove.push_back(it);
    }
    // Remove modifiers 
    for (auto &m : remove) modifiers.erase(m);
  }

  RealType GFlow::getKB() {
    return KB;
  }

}
