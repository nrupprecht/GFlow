/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#ifndef __GRAVITY_VELOCITY_VERLET_INTEGRATOR_HPP__
#define __GRAVITY_VELOCITY_VERLET_INTEGRATOR_HPP__

// Includes
#include "VelocityVerletIntegrator.hpp"

namespace GFlow {

  /*
   * @class VelocityVerletIntegrator
   * This integrator implements the Velocity Verlet integration scheme with gravity
   *
   */
  class GravityVelocityVerletIntegrator : public VelocityVerletIntegrator {
  public:

  private:
    virtual void firstHalfKick(); // Gravity is incorporated in the first half kick
  };

}
#endif // __GRAVITY_VELOCITY_VERLET_INTEGRATOR_HPP__
