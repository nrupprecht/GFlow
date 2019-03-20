#ifndef __HARMONIC_CHAIN_HPP__GFLOW__
#define __HARMONIC_CHAIN_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  class HarmonicChain : public Modifier {
  public:
    //! \brief Default constructor.
    HarmonicChain(GFlow*);

    //! \brief Constructor that sets spring constant.
    HarmonicChain(GFlow*, RealType);

    //! \brief Add another atom to the chain.
    virtual void addAtom(int);

    //! \brief Applies the forces to the atoms.
    virtual void post_forces();

  private:
    //! \brief Update the local id list from the global ids.
    void updateLocalIDs();

    //! \brief The global ids of the particles in the chain.
    vector<int> global_ids;

    //! \brief The local ids of the particles in the chain.
    vector<int> local_ids;

    //! \brief The relaxation distances for the bonds.
    vector<RealType> distance;

    //! \brief The strength of the harmonic bond, K.
    RealType springConstant;
  };

}
#endif // __HARMONIC_CHAIN_HPP__GFLOW__