#include "gflow.hpp"
#include "simdata.hpp"
#include "integrator.hpp"
#include "velocityverlet.hpp"
#include "force.hpp"
#include "sectorization.hpp"
#include "neighbors.hpp"
#include "communicator.hpp"
#include "datamaster.hpp"
#include "modifier.hpp"

namespace GFlowSimulation {

  GFlow::GFlow() : running(false), elapsed_time(0), iter(0) {
    simData = new SimData(this);
    // Integrator will be created by the creator
    sectorization = new Sectorization(this);
    communicator = new Communicator(this);
    dataMaster = new DataMaster(this);

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
    base->simData = simData;
    base->integrator = integrator;
    base->sectorization = sectorization;
    base->communicator = communicator;
  }

  void GFlow::run(RealType requested_time) {
    // Make sure we have initialized everything
    if (!initialize()) {
      // Some object was null
      cout << "Some object was null. Exiting.\n";
      return;
    }

    // Pre-integrate
    running = true;
    elapsed_time = 0;
    iter = 0;

    // Do integration for the requested amount of time
    while (running) {
      // --> Pre-step
      integrator->pre_step();
      dataMaster->pre_step();
      for (auto m : modifiers) m->pre_step();

      // --> Pre-exchange
      // Wrap positions
      wrapPositions();
      integrator->pre_exchange();
      dataMaster->pre_exchange();
      for (auto m : modifiers) m->pre_exchange();

      // --- Exchange particles (potentially) ---

      // --> Pre-neighbors
      integrator->pre_neighbors();
      dataMaster->pre_neighbors();
      for (auto m : modifiers) m->pre_neighbors();

      // --- Do neighbor updates (potentially)

      // --> Pre-force
      integrator->pre_forces(); // -- This is where VV first half kick happens
      dataMaster->pre_forces();
      for (auto m : modifiers) m->pre_forces();
      // Clip positions (if necessary)

      // --- Do forces

      // --> Post-forces
      integrator->post_forces(); // -- This is where VV second half kick happens
      dataMaster->post_forces();
      for (auto m : modifiers) m->post_forces();

      // --> Post-step
      if (requested_time<=elapsed_time) running = false;
      integrator->post_step();
      dataMaster->post_step();
      for (auto m : modifiers) m->post_step();
      // Timer updates
      ++iter;
      elapsed_time += integrator->getTimeStep();
    }

    // Post-integrate

  }

  void GFlow::writeData(string dirName) {
    // Make sure data master is non-null, print a warning if any writes failed
    if (dataMaster && !dataMaster->writeToDirectory(dirName))
      cout << "Warning: Some writes failed.\n";
  }

  RealType GFlow::getElapsedTime() {
    return elapsed_time;
  }

  void GFlow::setAllWrap(bool w) {
    for (int d=0; d<DIMENSIONS; ++d) wrap[d] = w;
  }

  inline void GFlow::wrapPositions() {
    // Get a pointer to position data and the number of particles in simData
    RealType **x = simData->x;
    int number = simData->number;
    // Wrap, if applicable
    for (int d=0; d<DIMENSIONS; ++d)
      if (wrap[d]) { // Wrap the d-th dimension
        // The width of the simulation in the d-th dimension
        RealType dim = bounds.max[d] - bounds.min[d];
        // Wrap each particle
        for (int n=0; n<number; ++n) {
          if (x[n][d]>bounds.max[d])       x[n][d] -= dim;
          else if (x[n][d]<=bounds.min[d]) x[n][d] += dim;
        }
      }
  }

}