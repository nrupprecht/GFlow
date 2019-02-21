#ifndef __BOND_HPP__GFLOW__
#define __BOND_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  class Bond : public Modifier {
  public:
    //! \brief Default constructor.
    Bond(GFlow*);

    //! \brief Add a bond - basic version.
    virtual void addBond(int, int);

  private:
    //! \brief The left and right particles of the bond.
    vector<int> left, center, right;
  };

}
#endif // __BOND_HPP__GFLOW__