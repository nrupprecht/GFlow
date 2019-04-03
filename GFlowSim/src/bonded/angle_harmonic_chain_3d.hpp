#ifndef __ANGLE_HARMONIC_CHAIN_3D_HPP__GFLOW__
#define __ANGLE_HARMONIC_CHAIN_3D_HPP__GFLOW__

#include "angle-harmonic-chain.hpp"

namespace GFlowSimulation {

  class AngleHarmonicChain_3d : public AngleHarmonicChain {
  public:
    AngleHarmonicChain_3d(GFlow*);
    AngleHarmonicChain_3d(GFlow*, RealType);

    virtual void interact() const override;
  };

}
#endif // __ANGLE_HARMONIC_CHAIN_3D_HPP__GFLOW__