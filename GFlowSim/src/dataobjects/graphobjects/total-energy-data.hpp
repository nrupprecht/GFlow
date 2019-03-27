#ifndef __TOTAL_ENERGY_DATA_HPP__GFLOW__
#define __TOTAL_ENERGY_DATA_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"

namespace GFlowSimulation {

  class TotalEnergyData : public GraphObject {
  public:
    //! \brief Default constructor.
    TotalEnergyData(GFlow*, bool=true);

    //! \brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step() override;

  private:
    //! \brief Whether to use the average energy per particle (as opposed to the total energy).
    bool useAve;
  };

}
#endif // __TOTAL_ENERGY_DATA_HPP__GFLOW__