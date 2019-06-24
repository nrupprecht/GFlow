#include "volumeplotobject2d.hpp"

namespace GFlowSimulation {

  VolumePlotObject2D::VolumePlotObject2D(GFlow *gflow, const string& name) : DataObject(gflow, name, DataObjectType::VOLUMEPLOT) {};

}