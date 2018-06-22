/*
 * Author: Nathaniel Rupprecht
 * Start Data: January 29, 2018
 *
 */

#ifndef __TEMPERATURE_FORCE_HPP__
#define __TEMPERATURE_FORCE_HPP__

#include "ExternalForce.hpp"

namespace GFlow {
  
  /*
   * @class ConstantAcceleration
   * A random force that simulates temperature
   *
   */
  class TemperatureForce : public ExternalForce {
  public:
    // Constructor
    TemperatureForce(RealType); // Temperature
    TemperatureForce(vec2, RealType, RealType); // Position, radius, temperature
    TemperatureForce(RealType, RealType); // Temperature and viscosity

  protected:
    // Private virtual functions, must overload
    virtual void _applyForce(SimData*) const;
    virtual string _summary() const;

    // The temperature and viscosity
    RealType temperature; // Temperature
    RealType viscosity;   // Viscosity
    RealType radius;      // Radius around the position that the temperature applies
    RealType tempDelay;   // How long we should wait between perturbations
    vec2 position;
    mutable RealType lastUpdate;  // The last time a perturbation was applied
  };
};

#endif // __TEMPERATURE_FORCE_HPP__