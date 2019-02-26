#ifndef __BOND_HPP__GFLOW__
#define __BOND_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  class Bond : public Modifier {
  public:
    //! \brief Default constructor.
    Bond(GFlow*);

    //! \brief Add a bond - basic version.
    virtual void addBond(int, int)=0;

    //! \brief Where the bonds execute forces.
    virtual void post_forces()=0;

    //! \brief Exception class for the vectors of data are not the same size.
    class UnequalBondVectors : public Exception {};

  protected:

    //! \brief Check whether lists of bonds have the same size.
    bool checkBondVectors();

    //! \brief The left and right particles of the bond.
    vector<int> left, right;

    //! \brief The global ids of the left and right particles of the bond.
    vector<int> gleft, gright;
    
  };

}
#endif // __BOND_HPP__GFLOW__