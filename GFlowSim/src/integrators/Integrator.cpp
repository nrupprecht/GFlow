#include "Integrator.hpp"

namespace GFlow {
  Integrator::Integrator() : dt(default_epsilon), time(0), runTime(0), iter(0), maxIter(0), simData(nullptr), sectors(nullptr), dataRecord(nullptr) {};

  Integrator::Integrator(SimData *sim) : dt(default_epsilon), time(0), runTime(0), iter(0), maxIter(0), simData(sim), sectors(nullptr), dataRecord(nullptr) {};

  Integrator::Integrator(SimData *sim, DataRecord *dat) : dt(default_epsilon), time(0), runTime(0), iter(0), maxIter(0), simData(sim), sectors(nullptr), dataRecord(dat) {};

  Integrator::~Integrator() {
    delete sectors;
    delete forceHandler;
  }

  void Integrator::initialize(RealType rt) {
    // Set run time
    if (simData) runTime = rt;
    else runTime = 0;

    // Create a sectorization
    sectors = new Sectorization;

    // Set up the sectorization
    if (simData) sectors->initialize(simData);

    // Create force handler
    forceHandler = new ForceHandler;

    // Set max iter
    if (dt>0) maxIter = runTime/dt;
    else maxIter = 0;
  }

  void Integrator::integrate() {
    // Pre
    preIntegrate();
    // Do the integration
    _integrate();
    // Post
    postIntegrate();
  }

  void Integrator::setDataRecord(DataRecord* dr) {
    dataRecord = dr;
  }

  void Integrator::preIntegrate() {
    // Beginning data record stuff
    if (dataRecord) {
      dataRecord->setRunTime(runTime);
      dataRecord->startTiming();
    }
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
