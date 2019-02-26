#ifndef __HARMONIC_BOND_HPP__GFLOW__
#define __HARMONIC_BOND_HPP__GFLOW__

#include "bond.hpp"

namespace GFlowSimulation {

  class HarmonicBond : public Bond {
  public:
    //! \brief Default constructor.
    HarmonicBond(GFlow*);

    //! \brief Constructor that sets spring constant.
    HarmonicBond(GFlow*, RealType);

    virtual void addBond(int, int);

    virtual void post_forces();

  private:
    // Helper functions
    
    //! \brief Update the local id list from the global ids.
    void updateLocalIDs();

    //! \brief The relaxation distances for the bonds.
    vector<RealType> distance;

    //! \brief The strength of the harmonic bond, K.
    RealType springConstant;
  };

}
#endif // __HARMONIC_BOND_HPP__GFLOW__