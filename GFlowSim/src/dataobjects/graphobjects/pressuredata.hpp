#ifndef __PRESSURE_DATA_HPP__GFLOW__
#define __PRESSURE_DATA_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"

namespace GFlowSimulation {

  /** 
   *  \brief Graphs the instantaneous pressure of all particles
   */
  class PressureData : public GraphObject {
  public:
    //! \brief Constructor
    PressureData(GFlow*);

    //! \brief Turn on virial calculations for all interactions.
    virtual void pre_integrate() override;

    //! \brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step() override;

    //! \brief Calculate the kinetic energy of the particles in the simdata.
    static RealType calculate_pressure(GFlow*);
  };

}
#endif // __PRESSURE_DATA_HPP__GFLOW__