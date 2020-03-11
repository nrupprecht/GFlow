#include "cavity-positions.hpp"

namespace GFlowSimulation {

  CavityPositions::CavityPositions(GFlow *gflow, real v) : PositionData(gflow), limit_velocity(v) {
    dataName = "Cav";
    select_function = [=](shared_ptr<SimData> simdata, int n)->bool {
      return simdata->V(n, 0)<limit_velocity && 0<simdata->X(n, 0) && 0<simdata->Im(n);
    };
  };

}