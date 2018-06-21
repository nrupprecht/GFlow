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
    integrator = new VelocityVerlet(this);
    sectorization = new Sectorization(this);
    communicator = new Communicator(this);
    dataMaster = new DataMaster(this);
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
      // --- Pre-step
      integrator->pre_step();
      for (auto m : modifiers) m->pre_step();

      // --- Pre-exchange
      // Wrap positions
      wrapPositions();

      // --- Exchange particles (potentially) ---

      // --- Pre-neighbors

      // --- Do neighbor updates (potentially)

      // --- Pre-force
      integrator->pre_forces();
      for (auto m : modifiers) m->pre_forces();
      // Clip positions (if necessary)

      // --- Do forces

      // --- Post-forces
      integrator->post_forces();
      for (auto m : modifiers) m->post_forces();

      // --- Post-step
      if (requested_time<=elapsed_time) running = false;
      integrator->post_step();
      for (auto m : modifiers) m->post_step();
      // Timer updates
      ++iter;
      elapsed_time += integrator->getTimeStep();
    }

    // Post-integrate

  }

  RealType GFlow::getElapsedTime() {
    return elapsed_time;
  }

  inline void GFlow::wrapPositions() {
    
  }

}