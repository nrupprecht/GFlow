#ifndef __KINETIC_ENERGY_BIN_HPP__GFLOW__
#define __KINETIC_ENERGY_BIN_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"

namespace GFlowSimulation {

  class KineticEnergyBin : public GraphObject {
  public:
    //! \brief Default constructor
    KineticEnergyBin(GFlow*);

    //! \brief Default constructor
    KineticEnergyBin(GFlow*, RealType);

    void post_step() override {};

    //! \brief Look at the data once the simulation ends.
    void post_integrate() override;

  private:
    //! \brief Number of bins to use.
    int bins = 50;

    //! \brief Cutoff radius.
    RealType radius;

    //! \brief Reference point.
    Vec center;
  };

}

#endif // __KINETIC_ENERGY_BIN_HPP__GFLOW__