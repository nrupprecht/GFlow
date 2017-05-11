#include "VelocityVerletIntegrator.hpp"

namespace GFlow {
  VelocityVerletIntegrator::VelocityVerletIntegrator() : Integrator() {};
  
  void VelocityVerletIntegrator::_integrate() {
    // Reset iter
    iter = 0;
    // Initialize 
    initialize();
    // Do the integration loop
    while (iter<maxIter) {
      preStep();
      integrateStep();
      postStep();
    }
  }
  
  void VelocityVerletIntegrator::preStep() {
    
  }
  
  void VelocityVerletIntegrator::integrateStep() {
    // First VV half kick
    
    // Calculate forces
    simData->doForces();
    // Exchange data with other processors
#ifdef USE_MPI
    simData->doMPI();
#endif
    // Second VV half kick
    
  }
  
  void VelocityVerletIntegrator::postStep() {
    // Update iteration
    ++iter;
    // Update time
    time += dt;
    // Update data recorder
    dataRecord->update(dt);
    dataRecord->record(simData);
  }
  
}
