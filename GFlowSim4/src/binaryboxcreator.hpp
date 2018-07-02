#ifndef __BINARY_BOX_CREATOR_HPP__GFLOW__
#define __BINARY_BOX_CREATOR_HPP__GFLOW__

#include "creator.hpp"

namespace GFlowSimulation {

  class BinaryBoxCreator : public Creator {
  public: 
    // Constructor
    BinaryBoxCreator(int, char**);

    // Constructor
    BinaryBoxCreator(ArgParse*);

    // Seed generators
    virtual void seedGenerator(uint);

    // Create simulation
    virtual GFlow* createSimulation();

  private:
    // Normal distribution
    std::mt19937 generator;
    std::normal_distribution<RealType> normal_dist;
    std::uniform_real_distribution<RealType> real_dist;
  };

}
#endif // __BINARY_BOX_CREATOR_HPP__GFLOW__