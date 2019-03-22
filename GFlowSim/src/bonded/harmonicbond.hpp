#ifndef __HARMONIC_BOND_HPP__GFLOW__
#define __HARMONIC_BOND_HPP__GFLOW__

#include "bond.hpp"

namespace GFlowSimulation {

  /**
  *  \brief The parent class for harmonic bonds.
  *
  *  The child classes are specific to the dimensionality of the system.
  */
  class HarmonicBond : public Bond {
  public:
    //! \brief Default constructor.
    HarmonicBond(GFlow*);

    //! \brief Constructor that sets spring constant.
    HarmonicBond(GFlow*, RealType);

    //! \brief Add a pair of bonded particles.
    virtual void addBond(int, int) override;

    //! \brief Calculate the interparticle forces.
    virtual void post_forces() override;

  protected:
    //! \brief The relaxation distances for the bonds.
    vector<RealType> distance;

    //! \brief The strength of the harmonic bond, K.
    RealType springConstant;
  };

}
#endif // __HARMONIC_BOND_HPP__GFLOW__