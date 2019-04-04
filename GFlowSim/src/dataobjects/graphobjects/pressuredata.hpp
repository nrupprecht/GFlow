#ifndef __PRESSURE_DATA_HPP__GFLOW__
#define __PRESSURE_DATA_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"

namespace GFlowSimulation {

  class PressureData : public GraphObject {
  public:
    //! Constructor
    PressureData(GFlow*);

    //! Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    //! \brief Calculate the kinetic energy of the particles in the simdata.
    static RealType calculate_pressure(GFlow*);

  };

}
#endif // __PRESSURE_DATA_HPP__GFLOW__