#include "gflow.hpp"
#include "simdata.hpp"
#include "integrator.hpp"
#include "velocityverlet.hpp"
#include "force.hpp"
#include "sectorization.hpp"
#include "verletlist.hpp"
#include "communicator.hpp"
#include "datamaster.hpp"
#include "forcemaster.hpp"
#include "modifier.hpp"

namespace GFlowSimulation {

  GFlow::GFlow() : running(false), requested_time(0), elapsed_time(0), total_time(0), iter(0), argc(0), argv(nullptr) {
    simData = new SimData(this);
    // Integrator will be created by the creator
    sectorization = new Sectorization(this);
    communicator = new Communicator(this);
    dataMaster = new DataMaster(this);
    forceMaster = new ForceMaster(this);

    // Set wrapping to true by defaule
    setAllWrap(true);
  }

  GFlow::~GFlow() {
    if (simData)       delete simData;
    if (sectorization) delete sectorization;
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

    if (sectorization) sectorization->initialize();
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
    base->simData       = simData;
    base->integrator    = integrator;
    base->sectorization = sectorization;
    base->communicator  = communicator;
    base->dataMaster    = dataMaster;
    base->forceMaster   = forceMaster;
    // Set vectors
    base->modifiersPtr  = &modifiers;
    base->forcesPtr     = &forces;
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
      cout << "Some object was null. Exiting.\n";
      return;
    }

    // --> Pre-integrate
    running = true;
    elapsed_time = 0;
    iter = 0;
    integrator->pre_integrate();
    dataMaster->pre_integrate();
    sectorization->pre_integrate();
    for (auto m : modifiers) m->pre_integrate();

    // Do integration for the requested amount of time
    while (running) {
      // --> Pre-step
      integrator->pre_step();
      dataMaster->pre_step();
      sectorization->pre_step();
      for (auto m : modifiers) m->pre_step();

      // --> Pre-exchange
      integrator->pre_exchange();
      dataMaster->pre_exchange();
      sectorization->pre_exchange();
      for (auto m : modifiers) m->pre_exchange();

      // --- Exchange particles (potentially) ---

      // --> Pre-force
      integrator->pre_forces(); // -- This is where VV first half kick happens

      dataMaster->pre_forces();
      sectorization->pre_forces(); // -- This is where resectorization / verlet list creation might happen
      for (auto m : modifiers) m->pre_forces(); // -- This is where modifiers should do forces (if they need to)

      // --- Do forces
      for (auto &f : forces) f->calculateForces();

      // --> Post-forces
      integrator->post_forces(); // -- This is where VV second half kick happens
      dataMaster->post_forces();
      sectorization->post_forces();
      for (auto m : modifiers) m->post_forces();

      // --- Clear force bufferes
      clearForces();

      // --> Post-step
      if (requested_time<=elapsed_time) running = false;
      integrator->post_step();
      dataMaster->post_step();
      sectorization->post_step();
      for (auto m : modifiers) m->post_step();
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
    sectorization->post_integrate();
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

  int GFlow::getIter() const {
    return iter;
  }

  Bounds GFlow::getBounds() const {
    return bounds;
  }

  const bool* GFlow::getWrap() const {
    return wrap;
  }

  pair<int, const char**> GFlow::getCommand() const {
    return pair<int, const char**>(argc, argv);
  }

  void GFlow::setCommand(int argc, char **argv) {
    if (argv) {
      this->argc = argc;
      this->argv = argv;
    }
  }

  void GFlow::setAllWrap(bool w) {
    for (int d=0; d<DIMENSIONS; ++d) wrap[d] = w;
  }

  void GFlow::requestTime(RealType t) {
    requested_time = t;
  }

  void GFlow::wrapPositions() {
    // Get a pointer to position data and the number of particles in simData
    RealType **x = simData->x;
    int number = simData->number;

    // Wrap all particles
    for (int d=0; d<DIMENSIONS; ++d) {
      if (wrap[d]) { // Wrap the d-th dimension
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

  void GFlow::addDataObject(class DataObject* dob) {
    dataMaster->addDataObject(dob);
  }

  inline void GFlow::clearForces() {
    simData->clearF();
  }

}
