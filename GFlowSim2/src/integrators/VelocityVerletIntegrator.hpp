/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#ifndef __VELOCITY_VERLET_INTEGRATOR_HPP__
#define __VELOCITY_VERLET_INTEGRATOR_HPP__

// Includes
#include "Integrator.hpp"

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

  private:
    // Private virtual functions
    virtual void _integrate();
    virtual void preStep();
    virtual void integrateStep();
    virtual void postStep();
    virtual void firstHalfKick();
    virtual void secondHalfKick();

  };

}
#endif // __VELOCITY_VERLET_INTEGRATOR_HPP__
