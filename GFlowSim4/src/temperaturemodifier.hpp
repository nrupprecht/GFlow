#ifndef __TEMPERATURE_MODIFIER_HPP__GFLOW__
#define __TEMPERATURE_MODIFIER_HPP__GFLOW__

#include "modifier.hpp"

namespace GFlowSimulation {

  class TemperatureModifier : public Modifier {
  public:
    // Constructor
    TemperatureModifier(GFlow*, RealType);

    // Constructor
    TemperatureModifier(GFlow*, RealType, RealType);

    // Apply forces
    virtual void post_forces();

  private:
    // Parameters
    RealType temperature, viscosity;

    // So we don't apply the force all the time - this is taxing because of all the random numbers
    // we need to generate
    RealType lastUpdate, updateDelay;
  };
}
#endif // __TEMPERATURE_MODIFIER_HPP__GFLOW__