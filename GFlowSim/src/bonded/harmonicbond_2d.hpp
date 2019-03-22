#ifndef __HARMONIC_BOND__2D_HPP__GFLOW__
#define __HARMONIC_BOND__2D_HPP__GFLOW__

#include "harmonicbond.hpp"

namespace GFlowSimulation {

  class HarmonicBond_2d : public HarmonicBond {
  public:
    //! \brief Default constructor.
    HarmonicBond_2d(GFlow*);

    //! \brief Constructor that sets spring constant.
    HarmonicBond_2d(GFlow*, RealType);

    //! \brief Calculate the interparticle forces.
    virtual void post_forces() override;
  };

}
#endif // __HARMONIC_BOND__2D_HPP__GFLOW__