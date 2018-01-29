/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#ifndef __VELOCITY_VERLET_INTEGRATOR_HPP__
#define __VELOCITY_VERLET_INTEGRATOR_HPP__

// Includes
#include "Integrator.hpp"
#include "../forces/ExternalForce.hpp"

namespace GFlow {

  /*
   * @class VelocityVerletIntegrator
   * This integrator implements the Velocity Verlet integration scheme
   *  
   */
  class VelocityVerletIntegrator : public Integrator {
  public:
    // Default constructor
    VelocityVerletIntegrator();

    // SimData initializing constructor
    VelocityVerletIntegrator(SimData*);

    // SimData and DataRecord initializing constructor
    VelocityVerletIntegrator(SimData*, DataRecord*);

    // Add a drag force
    void addExternalForce(ExternalForce*);

    // Set whether to adjust update delay or not
    void setAdjustUpdateDelay(bool b) { adjustUpdateDelay = b; }

    // Set whether to adjust the time step or not
    void setAdjustTimeStep(bool b) { adjustTimeStep = b; }

    // Set minimum timestep
    void setMinTimeStep(RealType m) { minTimestep = m; }

    // Set maximum timestep
    void setMaxTimeStep(RealType m) { maxTimestep = m; }

    // Get the update delay
    RealType getUpdateDelay() { return updateDelay; }

    // Get the average update delay
    RealType getAveUpdateDelay();

    // Get the average timestep
    RealType getAveTimeStep();

    // Get the move ratio
    RealType getMvRatio();

  private:
    // Inherited private virtual functions (from Integrator)
    virtual void preIntegrate();
    virtual inline void _integrate();
    virtual void preStep();
    virtual void postStep();

    // Private virtual functions
    virtual inline void _init();
    virtual inline void firstHalfKick();
    virtual inline void updates();
    virtual inline void forces();
    virtual inline void secondHalfKick();

    // Adjustment functions
    inline void doAdjustDelay();
    inline void doAdjustTimeStep();
    
    // Delay between updating sectors
    RealType updateDelay;
        
    // Delay multiplying factor
    RealType delayFactor;
    
    // Time since last verlet list creation
    RealType verletListTimer;
    
    // How long it's been since we updated the sectors
    RealType updateTimer;
    
    // How long it's been since we updated the time step
    RealType adjustTimer;
    
    // Whether to auto adjust update delay
    bool adjustUpdateDelay;
    
    // Whether to auto adjust the time step
    bool adjustTimeStep;
    RealType maxTimestep;
    RealType minTimestep;
    
    // How many time steps we want it to take for a particle to move a distance equal to its radius
    int periodIterations;
    
    // Number of updates and average delay and timestep
    RealType aveUpdateDelay;
    int sectorUpdates;
    
    // Movement ratio
    RealType mvRatio;
    
    // External force objects - we are not responsible for cleaning these up, simData is
    vector<ExternalForce*> externalForces;
  };

}
#endif // __VELOCITY_VERLET_INTEGRATOR_HPP__
