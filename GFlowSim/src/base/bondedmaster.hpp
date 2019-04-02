#ifndef __BONDED_MASTER_HPP__GFLOW__
#define __BONDED_MASTER_HPP__GFLOW__

#include "../gflow.hpp"

namespace GFlowSimulation {

  /**
  *  \brief Just as ForceMaster oversees short range interactions, BondedMaster oversees bonded interactions. This would
  *  include bonds, angles, dihedrals, chains, etc.
  *
  */
  class BondedMaster : public Base{
  public:
    //! \brief Default constructor.
    BondedMaster(GFlow*);

  private:
    //! \brief A vector of bonded objects.
    vector<Bonded*> bondedInteractions;
  };

}

#endif // __BONDED_MASTER_HPP__GFLOW__