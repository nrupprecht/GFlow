#ifndef __ANGLE_HARMONIC_CHAIN_HPP__GFLOW__
#define __ANGLE_HARMONIC_CHAIN_HPP__GFLOW__

#include "../base/bonded.hpp"

namespace GFlowSimulation {

  class AngleHarmonicChain : public Bonded {
  public:
    //! \brief Default constructor.
    AngleHarmonicChain(GFlow *);

    //! \brief Constructor that sets spring constant.
    AngleHarmonicChain(GFlow*, RealType);

    //! \brief Add another atom to the chain.
    virtual void addAtom(int);

    //! \brief Get the number of bonds.
    virtual int size() const override;

  protected:
    //! \brief Update the local id list from the global ids.
    void updateLocalIDs() const;

    //! \brief The global ids of the particles in the chain.
    vector<int> global_ids;

    //! \brief The local ids of the particles in the chain.
    mutable vector<int> local_ids;

    //! \brief The relaxation distances for the bonds.
    vector<RealType> distance;

    //! \brief The strength of the harmonic bond.
    RealType springConstant;

    //! \brief The strength of the angle.
    RealType angleConstant;
  };

}
#endif // __ANGLE_HARMONIC_CHAIN_HPP__GFLOW__