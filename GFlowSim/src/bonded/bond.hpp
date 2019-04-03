#ifndef __BOND_HPP__GFLOW__
#define __BOND_HPP__GFLOW__

#include "../base/bonded.hpp"

namespace GFlowSimulation {

  class Bond : public Bonded {
  public:
    //! \brief Default constructor.
    Bond(GFlow*);

    //! \brief Make sure local ids are set.
    virtual void pre_integrate() override;

    //! \brief Add a bond - basic version.
    virtual void addBond(int, int)=0;

    //! \brief Get the number of bonds.
    virtual int size() const override;

    //! \brief Exception class for the vectors of data are not the same size.
    class UnequalBondVectors : public Exception {};

  protected:

    //! \brief Check whether lists of bonds have the same size.
    virtual bool checkBondVectors() const;

    //! \brief Update the local ids.
    void updateLocalIDs() const;

    //! \brief The global ids of the left and right particles in a bond.
    vector<int> gleft, gright;

    //! \brief The local ids of the left and right particles in a bond.
    mutable vector<int> left, right;
    
  };

}
#endif // __BOND_HPP__GFLOW__