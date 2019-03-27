#include "avevelocitydata.hpp"
// Other files
#include "../../base/simdata.hpp"
#include "../../utility/vectormath.hpp"
#include "../../visualization/palette.hpp"

namespace GFlowSimulation {
  // Constructor
  AveVelocityData::AveVelocityData(GFlow *gflow) : GraphObject(gflow, "AveV", "time", "average velocity") {};

  void AveVelocityData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Get and store data
    RealType av = 0;
    RealType **v = Base::simData->V();
    RealType *im = Base::simData->Im();
    int size = Base::simData->size(), *type = Base::simData->Type();
    int count = 0;
    for (int n=0; n<size; ++n)
      if (im[n]>0 && type[n]>-1) {
        av += magnitudeVec(v[n], sim_dimensions);
        ++count;
      }
    // If we want the average
    av /= count;
    // Store data
    RealType time = Base::gflow->getElapsedTime();
    data.push_back(RPair(time, av));
  }

}