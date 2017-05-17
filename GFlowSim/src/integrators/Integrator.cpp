#include "Integrator.hpp"

namespace GFlow {
  
  Integrator::Integrator() : running(false), dt(default_epsilon), time(0), runTime(0), iter(0), simData(nullptr), sectors(nullptr), dataRecord(nullptr) {};

  Integrator::Integrator(SimData *sim) : running(false), dt(default_epsilon), time(0), runTime(0), iter(0), simData(sim), sectors(nullptr), dataRecord(nullptr) {};

  Integrator::Integrator(SimData *sim, DataRecord *dat) : running(false), dt(default_epsilon), time(0), runTime(0), iter(0), simData(sim), sectors(nullptr), dataRecord(dat) {};

  Integrator::~Integrator() {
    delete sectors;
    delete forceHandler;
  }

  void Integrator::initialize(RealType rt) {
    // Set run time
    if (simData) runTime = rt;
    else runTime = 0;

    // Create a sectorization
    if (sectors) delete sectors;
    sectors = new Sectorization;

    // Set up the sectorization
    if (simData) sectors->initialize(simData);

    // Create force handler
    forceHandler = new ForceHandler;

  }

  void Integrator::integrate() {
    // Pre
    preIntegrate();

    // Do the integration
    while (running) {
      preStep();    // Before doing an integration update
      _integrate(); // Integration update
      postStep();   // After doing an integration update
    }

    // Post
    postIntegrate();
  }

  void Integrator::setDataRecord(DataRecord* dr) {
    dataRecord = dr;
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
      dataRecord->markTime();
      dataRecord->setRunTime(runTime);
      dataRecord->startTiming();
    }
  }

  void Integrator::preStep() {
    // Increment times - do this here so if we change the timestep later on in integrators, we still have the correct time (without having to store a 'last dt' variable)
    time += dt;
  };

  void Integrator::postStep() {
    // Increment iter 
    ++iter;

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
