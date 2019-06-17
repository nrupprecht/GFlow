#ifndef __ANGLE_HARMONIC_CHAIN_HPP__GFLOW__
#define __ANGLE_HARMONIC_CHAIN_HPP__GFLOW__

// Useful: "Determination of Forces from a Potential in Molecular Dynamics"

#include "../base/bonded.hpp"
#include "../other/group.hpp"

namespace GFlowSimulation {

  class AngleHarmonicChain : public Bonded, public Group {
  public:
    //! \brief Default constructor.
    AngleHarmonicChain(GFlow *);

    //! \brief Constructor that sets spring constant.
    AngleHarmonicChain(GFlow*, RealType);

    //! \brief Constructor that sets spring and angle constants.
    AngleHarmonicChain(GFlow*, RealType, RealType);

    //! \brief Add another atom to the chain.
    virtual void addAtom(int);

    //! \brief Get the number of bonds.
    virtual int size() const override;

    //! \brief Set the spring constant.
    void setSpringConstant(RealType);

    //! \brief Set the spring constant
    void setAngleConstant(RealType);

  protected:

    //! \brief The relaxation distances for the bonds.
    vector<RealType> distance;

    //! \brief The strength of the harmonic bond.
    RealType springConstant;

    //! \brief The strength of the angle.
    RealType angleConstant;
  };

}
#endif // __ANGLE_HARMONIC_CHAIN_HPP__GFLOW__