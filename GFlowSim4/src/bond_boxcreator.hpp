#ifndef __BOND_BOX_CREATOR_HPP__GFLOW__
#define __BOND_BOX_CREATOR_HPP__GFLOW__

#include "creator.hpp"

namespace GFlowSimulation {

  class BondBoxCreator : public Creator {
  public:
    // Constructor
    BondBoxCreator(int, char**);

    // Constructor
    BondBoxCreator(ArgParse*);

    // Seed generators
    virtual void seedGenerator(uint);

    // Create simulation
    virtual GFlow* createSimulation();

  private:
    // Normal distribution
    std::mt19937 generator;
    std::normal_distribution<RealType> normal_dist;
  };

}
#endif // __BOND_BOX_CREATOR_HPP__GFLOW__

