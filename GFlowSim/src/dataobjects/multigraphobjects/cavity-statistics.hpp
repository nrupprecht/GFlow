#ifndef __CAVITY_STATISTICS_HPP__GFLOW__
#define __CAVITY_STATISTICS_HPP__GFLOW__

#include "../dataobjecttypes/multigraphobject.hpp"
#include "../../compute/store_data.hpp"

namespace GFlowSimulation {

  class CavityStatistics : public MultiGraphObject {
  public:
    CavityStatistics(GFlow*, real);

    virtual void pre_integrate() override;
    virtual void post_step() override;

  private:
    //! \brief A store data object.
    StoreData storeData;
    //! \brief Only collect data from particles with x velocity components less than this value.
    real limit_velocity;
    //! \brief The select function that is passed into store data when it collects data.
    std::function<bool(std::shared_ptr<SimData>, int)> select_function;
  };

}
#endif // __CAVITY_STATISTICS_HPP__GFLOW__