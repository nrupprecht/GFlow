#ifndef __HARMONIC_CHAIN_HPP__GFLOW__
#define __HARMONIC_CHAIN_HPP__GFLOW__

#include "../base/bonded.hpp"

namespace GFlowSimulation {

  class HarmonicChain : public Bonded {
  public:
    //! \brief Default constructor.
    HarmonicChain(GFlow*);

    //! \brief Constructor that sets spring constant.
    HarmonicChain(GFlow*, RealType);

    //! \brief Add another atom to the chain.
    virtual void addAtom(int);

    //! \brief Get the number of bonds.
    virtual int size() const override;

    //! \brief Bonded interactions.
    virtual void interact() const override;

  private:
    //! \brief Update the local id list from the global ids.
    void updateLocalIDs() const;

    //! \brief The global ids of the particles in the chain.
    vector<int> global_ids;

    //! \brief The local ids of the particles in the chain.
    mutable vector<int> local_ids;

    //! \brief The relaxation distances for the bonds.
    vector<RealType> distance;

    //! \brief The strength of the harmonic bond, K.
    RealType springConstant;
  };

}
#endif // __HARMONIC_CHAIN_HPP__GFLOW__