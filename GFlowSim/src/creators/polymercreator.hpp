#ifndef __POLYMER_CREATOR__HPP__GFLOW__
#define __POLYMER_CREATOR__HPP__GFLOW__

#include "area-creator.hpp"

namespace GFlowSimulation {

  class PolymerCreator : public AreaCreator {
  public:
    //! \brief Create a chain of small, noninteracting particles with random large, interacting particles.
    //!
    //! The chain maintains cohesion via harmonic bonds.
    virtual void createArea(HeadNode*, GFlow*, std::map<string, string>&) override;

  private:
    //! \brief Create a polymer chain.
    void createPolymer(GFlow*, int, RealType, RealType, RealType, int, int, bool);

    void createLine(HeadNode*, GFlow*, std::map<string, string>&);

    //! \brief A unified group correlation object. Note - this is handed to gflow, so this object should not attempt to delete it.
    class GroupCorrelation *correlation = nullptr;

    class HarmonicBond *harmonicbonds = nullptr;

    class AngleHarmonicChain *harmonicchain = nullptr;
  };

}
#endif // __POLYMER_CREATOR__HPP__GFLOW__