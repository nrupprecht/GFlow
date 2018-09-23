#ifndef __FLOW_CREATOR_HPP__GFLOW__
#define __FLOW_CREATOR_HPP__GFLOW__

#include "../base/creator.hpp"

namespace GFlowSimulation {

  //! @todo Commentary.
  class FlowCreator : public Creator {
  public:
    //! Constructor.
    FlowCreator(int, char**);

    //! Constructor -- pass in a pointer to an ArgParse object.
    FlowCreator(ArgParse*);

    //! Seed random number generators.
    virtual void seedGenerator(uint);

    //! Create a simulation.
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
  };

}
#endif // __FLOW_CREATOR_HPP__GFLOW__