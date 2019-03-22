#ifndef __HARMONIC_BOND__3D_HPP__GFLOW__
#define __HARMONIC_BOND__3D_HPP__GFLOW__

#include "harmonicbond.hpp"

namespace GFlowSimulation {

  class HarmonicBond_3d : public HarmonicBond {
  public:
    //! \brief Default constructor.
    HarmonicBond_3d(GFlow*);

    //! \brief Constructor that sets spring constant.
    HarmonicBond_3d(GFlow*, RealType);

    //! \brief Calculate the interparticle forces.
    virtual void post_forces() override;
  };
}
#endif // __HARMONIC_BOND__3D_HPP__GFLOW__