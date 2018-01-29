#include "VelocityVerletIntegrator.hpp"

namespace GFlow {
  VelocityVerletIntegrator::VelocityVerletIntegrator() : Integrator() {
    _init();
  }

  VelocityVerletIntegrator::VelocityVerletIntegrator(SimData* sim) : Integrator(sim) {
    _init();
  }

  VelocityVerletIntegrator::VelocityVerletIntegrator(SimData* sim, DataRecord* dat) : Integrator(sim, dat) {
    _init();
  }

  void VelocityVerletIntegrator::addExternalForce(ExternalForce* force) {
    if (simData==nullptr) return;
    simData->addExternalForce(force);
  }

  RealType VelocityVerletIntegrator::getAveUpdateDelay() { 
    return time/sectorUpdates;
  }

  // Get the average timestep
  RealType VelocityVerletIntegrator::getAveTimeStep() {
    return time/iter;
  }

  RealType VelocityVerletIntegrator::getMvRatio() {
    if (mvRatio==0) // If we are not updating delay
      return sectors->checkNeedRemake();
    return mvRatio;
  }

  void VelocityVerletIntegrator::preIntegrate() {
    // Make sure we are set up to run
    if (simData==nullptr || sectors==nullptr) {
      running = false;
      return;
    }

    // Do the normal pre-integration
    Integrator::preIntegrate();
  }
  
  void VelocityVerletIntegrator::preStep() {
    // Update timers
    updateTimer += dt;
    verletListTimer += dt;
    adjustTimer += dt;

    // Do the normal pre-step
    Integrator::preStep();
  }
  
  inline void VelocityVerletIntegrator::_integrate() {
    // First VV half kick
    firstHalfKick();
    
    // Update sectors, do MPI
    updates();

    // Calculate forces
    forces();
    
    // Second VV half kick
    secondHalfKick();
  }
  
  void VelocityVerletIntegrator::postStep() {
    // Record data from the integrator (before time is updated)
    if (dataRecord) dataRecord->record(this, time);

    // Do the normal post-step
    Integrator::postStep();

    // Update adjustments 
    if (updateDelay < adjustTimer) {
      // Rest timer
      adjustTimer = 0;
      // Update time step
      if (adjustTimeStep) doAdjustTimeStep();
    }
  }

  inline void VelocityVerletIntegrator::_init() {
    updateDelay = default_update_delay;
    delayFactor = default_delay_factor;
    verletListTimer = 0;
    updateTimer = 0;
    adjustTimer = 0;
    adjustUpdateDelay = true;
    adjustTimeStep = true;
    maxTimestep = default_max_timestep;
    minTimestep = default_min_timestep;
    periodIterations = default_period_iterations;
    aveUpdateDelay = 0;
    sectorUpdates = 0;
    mvRatio = 0;

    // Add a termination condition
    termination.push_back(new OutsideRegion);
  }

  inline void VelocityVerletIntegrator::firstHalfKick() {
    // Get data
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    RealType *fx = simData->getFxPtr();
    RealType *fy = simData->getFyPtr();
    RealType *th = simData->getThPtr();
    RealType *om = simData->getOmPtr();
    RealType *tq = simData->getTqPtr();
    RealType *im = simData->getImPtr();
    RealType *iI = simData->getIiPtr();

    // Get the number of particles we need to update
    int domain_end = simData->getDomainEnd();
    double hdt = 0.5*dt;
    if (simData->getWrapX() || simData->getWrapY()) {
#if _INTEL_ == 1
#pragma vector aligned
#pragma simd
#endif
#if _CLANG_ == 1
#pragma clang loop vectorize(enable)
#pragma clang loop interleave(enable)
#endif
      for (int i=0; i<domain_end; ++i) {
        // Update linear variables
        vx[i] += hdt*im[i]*fx[i];
	vy[i] += hdt*im[i]*fy[i];
	px[i] += dt*vx[i];
        py[i] += dt*vy[i];
        // Wrap position
	simData->wrap(px[i], py[i]);
        // Update angular variables
	om[i] += hdt*iI[i]*tq[i];
	th[i] += dt*om[i];
        // Wrap theta
	simData->wrap(th[i]);
      }
    }
    else {
#if _INTEL_ == 1
#pragma vector aligned
#pragma simd
#endif
#if _CLANG_ == 1
#pragma clang loop vectorize(enable)
#pragma clang loop interleave(enable)
#endif
      for (int i=0; i<domain_end; ++i) {
	// Update linear variables
	vx[i] += hdt*im[i]*fx[i];
	vy[i] += hdt*im[i]*fy[i];
	px[i] += dt*vx[i];
	py[i] += dt*vy[i];
	// Update angular variables
	om[i] += hdt*iI[i]*tq[i];
	th[i] += dt*om[i];
      }
    }
    // Clear force and torque
    simData->clearForceTorque();
  }

  inline void VelocityVerletIntegrator::updates() {
    // Do characteristics first
    for (auto &c : simData->getCharacteristics())
      c.second->modify(simData, c.first, dt);

    // Update sectors
    if (updateDelay < updateTimer) {
      // Check termination conditions
      for (auto & t : termination)
	if (t->check(simData)) {
	  running = false;
	  // Say why we ended
	  std::cerr << t->getMessage() << endl;
	  return;
	}

#if USE_MPI == 1 // If using mpi, exchange data with other processors
      simData->atomMove();
      simData->atomCopy();
#endif

      // Adjust update delay
      if (adjustUpdateDelay) doAdjustDelay();
      
      // Update lists and reset list timer
      sectors->createVerletLists();
      sectors->createWallLists();
      verletListTimer = 0;

      // Reset update timer
      updateTimer = 0;
    }
  }
  
  inline void VelocityVerletIntegrator::forces() {
    // Apply external forces
    for (const auto& force : simData->getExternalForces())
      force->applyForce(simData);

    // Do particle and wall forces
    const auto& verletList = sectors->getVerletList();
    const auto& wallList   = sectors->getWallList();
    forceHandler->pForces(verletList, simData);
    forceHandler->wForces(wallList, simData);
  }
  
  inline void VelocityVerletIntegrator::secondHalfKick() {
    // Get the neccessary data
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    RealType *fx = simData->getFxPtr();
    RealType *fy = simData->getFyPtr();
    RealType *om = simData->getOmPtr();
    RealType *tq = simData->getTqPtr();
    RealType *im = simData->getImPtr();
    RealType *iI = simData->getIiPtr();
    
    // Get the number of particles we need to update
    int domain_end = simData->getDomainEnd();
    // Do second half-kick
    RealType hdt = 0.5*dt;
    
#if _INTEL_ == 1
#pragma vector aligned
#pragma simd
#endif
#if _CLANG_ == 1
#pragma clang loop vectorize(enable)
#pragma clang loop interleave(enable)
#endif
    for (int i=0; i<domain_end; ++i) {
      vx[i] += hdt*im[i]*fx[i];
      vy[i] += hdt*im[i]*fy[i];
      om[i] += hdt*iI[i]*tq[i];
    }
  }
  
  inline void VelocityVerletIntegrator::doAdjustDelay() {
    // Calculate new update delay
    mvRatio = sectors->checkNeedRemake(); // Want movement to be slightly greater then skin depth after every update delay
    if (mvRatio==0) return;

    // Set the update delay
    updateDelay = delayFactor*verletListTimer/mvRatio;
    
    // If something suddenly starts moving, and the update delay is to long, there could be problems. For now, we deal with that by capping the update delay
    updateDelay = updateDelay>default_max_update_delay ? default_max_update_delay : updateDelay;
    
    // Record data for average update delay
    ++sectorUpdates;
    aveUpdateDelay += updateDelay;
  }

  inline void VelocityVerletIntegrator::doAdjustTimeStep() {
    // Find the minimum amount of time it takes for one particle to traverse its own radius
    RealType minPeriod = 1., dt1=1, dt2=1;
    int domain_end = simData->getDomainEnd();

    // Linear period finding
    for (int i=0; i<domain_end; ++i) {
      if (-1<simData->getIt(i)) {
	RealType period = sqr(simData->getSg(i))/(sqr(simData->getVx(i))+sqr(simData->getVy(i)));
	if (period<minPeriod) minPeriod = period;
      }
    }
    dt2 = sqrt(minPeriod)/periodIterations;

    // Set the time step
    //dt = minPeriod/periodIterations;
    dt = min(dt1, dt2);
    
    // If something suddenly starts moving, and the update delay is to long, there could be problems. For now, we deal with that by capping the time step
    dt = dt>maxTimestep ? maxTimestep : dt;
    dt = dt<minTimestep ? minTimestep : dt;
  }
  
}
