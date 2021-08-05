#ifndef __CAVITY_POSITIONS_HPP__GFLOW__
#define __CAVITY_POSITIONS_HPP__GFLOW__

#include "positiondata.hpp"

#include "dataobjecttypes/data-base/multigraphdata.hpp"

namespace GFlowSimulation {

  /**
  *  \brief A specialization of PositionData, used for looking at the cavity.
  *  
  **/
  class CavityPositions : public PositionData, public MultiGraphData {
  public:
    CavityPositions(GFlow*, real);

  private:
    //! \brief Only particles with v_x < limit_velocity will be recorded.
    real limit_velocity;
  };

}
#endif // __CAVITY_POSITIONS_HPP__GFLOW__