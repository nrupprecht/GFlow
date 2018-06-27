#ifndef __BOX_CREATOR_HPP__GFLOW__
#define __BOX_CREATOR_HPP__GFLOW__

#include "creator.hpp"

namespace GFlowSimulation {

  class BoxCreator : public Creator {
  public:
    // Constructor
    BoxCreator(int, char**);

    // Constructor
    BoxCreator(ArgParse*);

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
#endif