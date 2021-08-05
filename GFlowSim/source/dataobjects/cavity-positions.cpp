#include <dataobjects/cavity-positions.hpp>

using namespace GFlowSimulation;

CavityPositions::CavityPositions(GFlow *gflow, real v)
    : PositionData(gflow), MultiGraphData("time", "", 3), limit_velocity(v) {
  dataName = "Cav";
  select_function = [=](const auto& simdata, int n) -> bool {
    return simdata->V(n, 0) < 0.5f * limit_velocity && 0 < simdata->X(n, 0) && 0 < simdata->Im(n);
  };
};
