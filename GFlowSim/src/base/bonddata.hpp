#ifndef __BOND_DATA_HPP__GFLOW__
#define __BOND_DATA_HPP__GFLOW__

#include "../gflow.hpp"

namespace GFlowSimulation {

  class BondData : public Base {
  public:
    //! Constructor
    BondData(GFlow*);

    //! @brief Add a pair of particles for a bond, specifying the relaxed length of the bond.
    void addBond(int, int, RealType, RealType);

    //! @brief Add a pair of particles for a bond.
    void addBond(int, int, RealType);

    virtual void post_forces() override;

  private:
    //! @brief Indices, in pair, that represent bonded particles
    vector<int> pairs;
    //! @brief The strength of bonds.
    vector<int> strength;
    //! @brief The relaxed lengths.
    vector<RealType> length;
  };

}
#endif // __BOND_DATA_HPP__GFLOW__