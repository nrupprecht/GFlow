#include "Integrator.hpp"

namespace GFlow {
  
  Integrator::Integrator() : running(false), dt(default_epsilon), time(0), runTime(0), iter(0), simData(nullptr), sectors(nullptr), dataRecord(nullptr), forceHandler(nullptr) {};

  Integrator::Integrator(SimData *sim) : running(false), dt(default_epsilon), time(0), runTime(0), iter(0), simData(sim), sectors(nullptr), dataRecord(nullptr), forceHandler(nullptr) {};

  Integrator::Integrator(SimData *sim, DataRecord *dat) : running(false), dt(default_epsilon), time(0), runTime(0), iter(0), simData(sim), sectors(nullptr), dataRecord(dat), forceHandler(nullptr) {};

  Integrator::~Integrator() {
    delete sectors;
    delete forceHandler;
    // We are responsible for deleting these
    for (auto& t : termination) delete t;
  }

  void Integrator::initialize(RealType rt) {
    // Set run time
    if (simData) runTime = rt;
    else runTime = 0;

    // Create and set up a sectorization
    initializeSectors();

    // Make sure we have initial verlet lists
    sectors->createVerletLists();
    sectors->createWallLists();

    // Create and set up a force handler
    initializeForceHandler();
  }

  void Integrator::setDataRecord(DataRecord* dr) {
    dataRecord = dr;
  }

  void Integrator::integrate(RealType time) {
    // Initialize
    initialize(time);

    // Pre
    preIntegrate();

    // Do the integration
    while (running) step();

    // Post
    postIntegrate();
  }

  void Integrator::step() {
    preStep();    // Before doing an integration update
    _integrate(); // Integration update
    postStep();   // After doing an integration update
  }

  void Integrator::initializeSectors() {
    // Create a sectorization
    if (sectors) delete sectors;
    sectors = new StandardSectorization;
    // Set up the sectorization
    if (simData) sectors->initialize(simData);
  }

  void Integrator::initializeForceHandler() {
    // Create force handler
    forceHandler = new ForceHandler;

    // Give the force handler to the simulation data
    if (simData) simData->setForceHandler(forceHandler);
  }

  void Integrator::preIntegrate() {
    // Reset iter
    iter = 0;
    
    // Reset time
    time = 0;

    // Set running to true
    running = true;

    // Beginning data record stuff
    if (dataRecord) {
      dataRecord->initialize();
      dataRecord->markTime();
      dataRecord->setRunTime(runTime);
      dataRecord->startTiming();
      // Initial record
      dataRecord->record(simData, time);
    }
  }

  void Integrator::preStep() {};

  void Integrator::postStep() {
    // Increment iter 
    ++iter;

    // Update time
    time += dt;

    // End if we have simulated for enough time
    if (runTime < time) running = false;

    // End if the terminate flag has been set in our simulation data
    if (simData->getTerminate()) running = false;

    // Update data recorder
    if (dataRecord) dataRecord->record(simData, time);
  }

  void Integrator::postIntegrate() {
    // Ending data record stuff
    if (dataRecord) {
      dataRecord->endTiming();
      dataRecord->getSectorizationData(sectors);
      dataRecord->setActualTime(time);
    }
  }
  
}
