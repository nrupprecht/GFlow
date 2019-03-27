#ifndef __KINETIC_ENERGY_DATA_HPP__GFLOW__
#define __KINETIC_ENERGY_DATA_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"

namespace GFlowSimulation {

  class KineticEnergyData : public GraphObject {
  public:
    //! \brief Constructor
    KineticEnergyData(GFlow*, bool=true);

    //! \brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step() override;

  private:
    //! \brief Whether to use the average kinetic energy per particle (as opposed to the total kinetic energy).
    bool useAve; 
  };

}
#endif // __KINETIC_ENERGY_DATA_HPP__GFLOW__