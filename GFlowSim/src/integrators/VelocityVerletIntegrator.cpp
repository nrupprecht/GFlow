#include "VelocityVerletIntegrator.hpp"

namespace GFlow {
  VelocityVerletIntegrator::VelocityVerletIntegrator() : Integrator(), updateDelay(default_update_delay), delayFactor(default_delay_factor), verletListTimer(0), updateTimer(0), adjustUpdateDelay(true), adjustTimeStep(true), periodIterations(default_period_iterations), aveUpdateDelay(0), sectorUpdates(0) {};

  VelocityVerletIntegrator::VelocityVerletIntegrator(SimData* sim) : Integrator(sim), updateDelay(default_update_delay), delayFactor(default_delay_factor), verletListTimer(0), updateTimer(0), adjustUpdateDelay(true), adjustTimeStep(true), periodIterations(default_period_iterations), aveUpdateDelay(0), sectorUpdates(0) {};

  VelocityVerletIntegrator::VelocityVerletIntegrator(SimData* sim, DataRecord* dat) : Integrator(sim, dat), updateDelay(default_update_delay), delayFactor(default_delay_factor), verletListTimer(0), updateTimer(0), adjustUpdateDelay(true), adjustTimeStep(true), periodIterations(default_period_iterations), aveUpdateDelay(0), sectorUpdates(0) {};

  void VelocityVerletIntegrator::addExternalForce(ExternalForce* force) {
    if (simData==nullptr) return;
    simData->addExternalForce(force);
  }

  RealType VelocityVerletIntegrator::getAveUpdateDelay() { 
    return sectorUpdates>0 && adjustUpdateDelay ? aveUpdateDelay/sectorUpdates :updateDelay; 
  }

  // Get the average timestep
  RealType VelocityVerletIntegrator::getAveTimeStep() {
    return sectorUpdates>0 && adjustTimeStep ? aveDt/sectorUpdates : dt; 
  }

  void VelocityVerletIntegrator::preIntegrate() {
    // Make sure we are set up to run
    if (simData==nullptr || sectors==nullptr) {
      running = false;
      return;
    }

    // Make sure we have initial verlet lists
    sectors->createVerletLists(true);
    sectors->createWallLists(true);

    // Do the normal pre-integration
    Integrator::preIntegrate();
  }
  
  void VelocityVerletIntegrator::preStep() {
    // Update timers
    updateTimer += dt;
    verletListTimer += dt;

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
    // End if we have simulated for enough time
    if (runTime < time) running = false;

    // Do the normal post-step
    Integrator::postStep();
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
    int domain_size = simData->getDomainSize();
    double hdt = 0.5*dt;
#pragma vector aligned
#pragma simd
    for (int i=0; i<domain_size; ++i) {
      // Update linear variables
      vx[i] += hdt*im[i]*fx[i];
      vy[i] += hdt*im[i]*fy[i];
      px[i] += dt*vx[i];
      py[i] += dt*vy[i];
      // Wrap position
      simData->wrap(px[i], py[i]);
      // Zero force
      fx[i] = 0;
      fy[i] = 0;
      // Update angular variables
      om[i] += hdt*iI[i]*tq[i];
      th[i] += dt*om[i];
      // Wrap theta
      simData->wrap(th[i]);
      // Zero torque
      tq[i] = 0;
    }
  }

  inline void VelocityVerletIntegrator::updates() {
    // Update sectors
    if (updateDelay<updateTimer) {
      
#ifdef USE_MPI // If using mpi, exchange data with other processors
      simData->atomMove();
      simData->atomCopy();
#endif

      // Adjust update delay
      if (adjustUpdateDelay) doAdjustDelay();
      // Update time step
      if (adjustTimeStep) doAdjustTimeStep();
      
      // Update lists and reset list timer
      sectors->createVerletLists();
      sectors->createWallLists(true);
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
    int domain_size = simData->getDomainSize();
    // Do second half-kick
    RealType hdt = 0.5*dt;
#pragma vector aligned
#pragma simd
    for (int i=0; i<domain_size; ++i) {
      vx[i] += hdt*im[i]*fx[i];
      vy[i] += hdt*im[i]*fy[i];
      om[i] += hdt*iI[i]*tq[i];
    }
  }
  
  inline void VelocityVerletIntegrator::doAdjustDelay() {
    // Calculate new update delay
    RealType ratio = sectors->checkNeedRemake(); // Want movement to be slightly greater then skin depth after every update delay
    if (ratio==0) return;

    // Set the update delay
    updateDelay = delayFactor*verletListTimer/ratio;
    
    // If something suddenly starts moving, and the update delay is to long, there could be problems. For now, we deal with that by capping the update delay
    updateDelay = updateDelay>default_max_update_delay ? default_max_update_delay : updateDelay;
    
    // Record data for average update delay
    ++sectorUpdates;
    aveUpdateDelay += updateDelay;

    // Tell data recorder about this
    if (dataRecord) dataRecord->push_mvRatio(ratio);
  }

  inline void VelocityVerletIntegrator::doAdjustTimeStep() {
    // Find the minimum amount of time it takes for one particle to traverse its own radius
    RealType minPeriod = 1.;
    int domain_size = simData->getDomainSize();
    for (int i=0; i<domain_size; ++i) {
      if (-1<simData->getIt(i)) {
	RealType period = sqr(simData->getSg(i))/(sqr(simData->getVx(i))+sqr(simData->getVy(i)));
	if (period<minPeriod) minPeriod = period;
      }
    }
    minPeriod = sqrt(minPeriod);
    
    // Set the time step
    dt = minPeriod/periodIterations;
    // If something suddenly starts moving, and the update delay is to long, there could be problems. For now, we deal with that by capping the time step
    dt = dt>default_max_timestep ? default_max_timestep : dt;
    
    // Record data for average timestep data
    aveDt += dt;
  }
  
}
