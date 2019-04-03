#ifndef __ANGLE_HARMONIC_CHAIN_2D_HPP__GFLOW__
#define __ANGLE_HARMONIC_CHAIN_2D_HPP__GFLOW__

#include "angle-harmonic-chain.hpp"

namespace GFlowSimulation {

  class AngleHarmonicChain_2d : public AngleHarmonicChain {
  public:
    AngleHarmonicChain_2d(GFlow*);
    AngleHarmonicChain_2d(GFlow*, RealType);

    virtual void interact() const override;
  };

}
#endif // __ANGLE_HARMONIC_CHAIN_2D_HPP__GFLOW__