#ifndef __BIPARTITE_BOX_CREATOR_HPP__GFLOW__
#define __BIPARTITE_BOX_CREATOR_HPP__GFLOW__

#include "creator.hpp"

namespace GFlowSimulation {

  class BipartiteBoxCreator : public Creator {
  public: 
    // Constructor
    BipartiteBoxCreator(int, char**);

    // Constructor
    BipartiteBoxCreator(ArgParse*);

    // Seed generators
    virtual void seedGenerator(uint);

    // Create simulation
    virtual GFlow* createSimulation();

    void setPhi(RealType p) { phi = p; }
    void setWidth(RealType w) { width = w; }
    void setRadius(RealType r) { radius = r; }

  private:
    // Data
    RealType phi, width, radius;
    
    // Normal distribution
    std::mt19937 generator;
    std::normal_distribution<RealType> normal_dist;
    std::uniform_real_distribution<RealType> real_dist;
  };

}
#endif // __BIPARTITE_BOX_CREATOR_HPP__GFLOW__