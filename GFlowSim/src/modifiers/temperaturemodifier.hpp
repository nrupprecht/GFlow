#ifndef __TEMPERATURE_MODIFIER_HPP__GFLOW__
#define __TEMPERATURE_MODIFIER_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  class TemperatureModifier : public Modifier {
  public:
    //! @brief Constructor
    TemperatureModifier(GFlow*, RealType);

    //! @brief Constructor
    TemperatureModifier(GFlow*, RealType, RealType);

    //! @brief Apply forces
    virtual void post_forces() override;

  private:
    //! @brief Parameters
    RealType temperature, viscosity;

    // So we don't apply the force all the time - this is taxing because of all the random numbers
    // we need to generate
    RealType lastUpdate, updateDelay;
  };
}
#endif // __TEMPERATURE_MODIFIER_HPP__GFLOW__