#ifndef __POLYMER_CREATOR__HPP__GFLOW__
#define __POLYMER_CREATOR__HPP__GFLOW__

#include "area-creator.hpp"

namespace GFlowSimulation {

  class PolymerCreator : public AreaCreator {
  public:
    //! \brief Create a chain of small, noninteracting particles with random large, interacting particles.
    //!
    //! The chain maintains cohesion via harmonic bonds.
    virtual void createArea(HeadNode*, GFlow*, std::map<string, string>&);
  };

}
#endif // __POLYMER_CREATOR__HPP__GFLOW__