#include "gflow.hpp"
// Other files
#include "allbaseobjects.hpp"

namespace GFlowSimulation {

  GFlow::GFlow() : running(false), requested_time(0), elapsed_time(0), total_time(0), 
    iter(0), argc(0), argv(nullptr), repulsion(DEFAULT_HARD_SPHERE_REPULSION) 
  {
    simData      = new SimData(this);
    // Integrator will be created by the creator
    domain       = new Sectorization(this); // Domain(this); // Sectorization(this);
    communicator = new Communicator(this);
    dataMaster   = new DataMaster(this);
    forceMaster  = new ForceMaster(this);

    // Set wrapping to true by defaule
    setAllBCs(BCFlag::WRAP);
  }

  GFlow::~GFlow() {
    if (simData)       delete simData;
    if (domain)        delete domain;
    if (integrator)    delete integrator;
    if (communicator)  delete communicator;
    if (dataMaster)    delete dataMaster;
    for (auto &md : modifiers) 
      if (md) delete md;
    for (auto &fr : forces)
      if (fr) delete fr;
  }

  bool GFlow::initialize() {
    bool non_null = true;
    // Initialize all the subobjects
    if (simData) simData->initialize();
    else non_null = false;

    if(integrator) integrator->initialize();
    else non_null = false;

    if (domain) domain->initialize();
    else non_null = false;

    if (communicator) communicator->initialize();
    else non_null = false;

    if (dataMaster) dataMaster->initialize();
    else non_null = false;

    if (forceMaster) forceMaster->initialize();
    else non_null = false;

    for (auto &md : modifiers) {
      if (md) md->initialize();
      else non_null = false; 
    }
    for (auto &fr: forces) {
      if (fr) fr->initialize();
      else non_null = false;
    }
    // Return whether pointers were non-null
    return non_null;
  }

  void GFlow::initializeBase(Base *base) {
    // Make sure we aren't handed a null pointer
    if (base==nullptr) return;
    // Give pointer to this GFlow object
    base->gflow = this;
    // Set other objects
    base->simData      = simData;
    base->integrator   = integrator;
    base->domain       = domain;
    base->communicator = communicator;
    base->dataMaster   = dataMaster;
    base->forceMaster  = forceMaster;
    // Set vectors
    base->modifiersPtr = &modifiers;
    base->forcesPtr    = &forces;
  }

  void GFlow::run(RealType rt) {
    // If a parameter was passed in, it is the requested time
    if (rt>0) requested_time = rt;

    // Only run if time has been requested
    if (requested_time<=0) return;

    // Record this request
    total_requested_time += requested_time;

    // Make sure we have initialized everything
    if (!initialize()) {
      // Some object was null
      UnexpectedNullPointer("Error: Some object was null at GFlow initialization.");
    }

    // --> Pre-integrate
    running = true;
    elapsed_time = 0;
    iter = 0;
    integrator->pre_integrate();
    dataMaster->pre_integrate();
    domain->pre_integrate();
    for (auto m : modifiers) m->pre_integrate();

    // Do integration for the requested amount of time
    while (running) {
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

      // --- Exchange particles (potentially) ---
      domain->exchange_particles();

      // --> Pre-force
      for (auto m : modifiers) m->pre_forces();
      integrator->pre_forces(); // -- This is where VV first half kick happens
      dataMaster->pre_forces();

      domain->pre_forces(); // -- This is where resectorization / verlet list creation might happen
      // Wrap or reflect particles
      wrapPositions();
      reflectPositions();
      repulsePositions();

      // --- Do forces
      clearForces(); // Clear force buffers
      // Calculate current forces
      for (auto &f : forces) f->calculateForces();

      // --> Post-forces
      for (auto m : modifiers) m->post_forces(); // -- This is where modifiers should do forces (if they need to)
      integrator->post_forces(); // -- This is where VV second half kick happens
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
    }

    // --> Post-integrate
    requested_time = 0;
    integrator->post_integrate();
    dataMaster->post_integrate();
    domain->post_integrate();
    for (auto m : modifiers) m->post_integrate();
  }

  void GFlow::writeData(string dirName) {
    // Make sure data master is non-null, print a warning if any writes failed
    if (dataMaster && !dataMaster->writeToDirectory(dirName))
      cout << "Warning: Some writes failed.\n";
  }

  RealType GFlow::getRequestedTime() const {
    return requested_time;
  }

  RealType GFlow::getTotalRequestedTime() const {
    return total_requested_time;
  }

  RealType GFlow::getElapsedTime() const{
    return elapsed_time;
  }

  RealType GFlow::getTotalTime() const {
    return total_time;
  }

  RealType GFlow::getBoundaryForce() const {
    return boundaryForce;
  }

  RealType GFlow::getDT() const {
    return integrator->getTimeStep();
  }

  int GFlow::getIter() const {
    return iter;
  }

  int GFlow::getNumForces() const { 
    return forces.size(); 
  }

  Bounds GFlow::getBounds() const {
    return bounds;
  }

  const BCFlag* GFlow::getBCs() const {
    return boundaryConditions;
  }

  BCFlag GFlow::getBC(int dim) const {
    if (dim<0 || DIMENSIONS<=dim) 
      throw BadDimension("Bad dim in get BC.");
    return boundaryConditions[dim];
  }

  int GFlow::getNTypes() const {
    return forceMaster->getNTypes();
  }

  pair<int, char**> GFlow::getCommand() const {
    return pair<int, char**>(argc, argv);
  }

  const vector<class Force*>& GFlow::getForces() const {
    return forces;
  }

  void GFlow::setCommand(int argc, char **argv) {
    if (argv) {
      this->argc = argc;
      this->argv = argv;
    }
  }

  void GFlow::setBounds(Bounds b) {
    bounds = b;
    // Sectorization updates its bounds in pre_integrate()
  }

  void GFlow::setAllBCs(BCFlag type) {
    for (int d=0; d<DIMENSIONS; ++d) 
      boundaryConditions[d] = type;
  }

  void GFlow::setRepulsion(RealType r) {
    if (r<0) return;
    repulsion = r;
  }

  void GFlow::requestTime(RealType t) {
    if (t<0) t = 0;
    requested_time = t;
  }

  void GFlow::wrapPositions() {
    // Get a pointer to position data and the number of particles in simData
    RealType **x = simData->x;
    int number = simData->number;

    // Wrap all particles
    for (int d=0; d<DIMENSIONS; ++d) {
      if (boundaryConditions[d]==BCFlag::WRAP) { // Wrap the d-th dimension
        for (int n=0; n<number; ++n) {
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
  }

  void GFlow::reflectPositions() {
    // Get a pointer to position data and the number of particles in simData
    RealType **x = simData->x, **v = simData->v;
    int number = simData->number;

    // Reflect all the particles
    for (int d=0; d<DIMENSIONS; ++d)
      if (boundaryConditions[d]==BCFlag::REFL) { 
        for (int n=0; n<number; ++n) {
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
    RealType **x = simData->x, **f = simData->v;
    int number = simData->number;
    // Reset boundary force
    boundaryForce = 0;
    // Reflect all the particles
    for (int d=0; d<DIMENSIONS; ++d)
      if (boundaryConditions[d]==BCFlag::REPL) { 
        for (int n=0; n<number; ++n) {
          // Create a local copy
          if (x[n][d]<bounds.min[d]) {
            f[n][d] += repulsion*(bounds.min[d] - x[n][d]);
            boundaryForce += repulsion*(bounds.min[d] - x[n][d]);
          }
          else if (bounds.max[d]<x[n][d]) {
            f[n][d] -= repulsion*(x[n][d] - bounds.max[d]);
            boundaryForce += repulsion*(x[n][d] - bounds.max[d]);
          }
        }
      }
  }

  void GFlow::addDataObject(class DataObject* dob) {
    dataMaster->addDataObject(dob);
  }

  void GFlow::addModifier(class Modifier* mod) {
    modifiers.push_back(mod);
  }

  void GFlow::resetAllTimes() {
    requested_time       = 0.;
    total_requested_time = 0.;
    elapsed_time         = 0.;
    total_time           = 0.;
    iter                 = 0 ;
  }

  void GFlow::setStartRecTime(RealType t) {
    dataMaster->setStartRecTime(t);
  }

  void GFlow::setFPS(RealType fps) {
    dataMaster->setFPS(fps);
  }

  void GFlow::setFPS(int dob_id, RealType fps) {
    dataMaster->setFPS(dob_id, fps);
  }

  void GFlow::setDT(RealType dt) {
    integrator->setDT(dt);
  }

  void GFlow::setDMCmd(int argc, char** argv) {
    dataMaster->setCommand(argc, argv);
  }

  inline void GFlow::clearForces() {
    simData->clearF();
  }

}
