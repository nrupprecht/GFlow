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

  private:
    // Inherited private virtual functions
    virtual void _integrate(); // Inherits from Integrator

    // Private virtual functions
    virtual void preStep();
    virtual void integrateStep();
    virtual void postStep();
    virtual void firstHalfKick();
    virtual void secondHalfKick();

    // Delay between updating sectors
    RealType updateDelay;

    // Delay multiplying factor
    RealType delayFactor;

    // Time since last verlet list creation
    RealType verletListTimer;

    // How long it's been since we updated the sectors
    RealType updateTimer;

    // Drag force objects
    vector<ExternalForce*> externalForces;

  };

}
#endif // __VELOCITY_VERLET_INTEGRATOR_HPP__
