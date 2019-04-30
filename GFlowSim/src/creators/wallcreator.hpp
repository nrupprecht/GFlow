#ifndef __WALL_CREATOR_HPP__GFLOW__
#define __WALL_CREATOR_HPP__GFLOW__

#include "area-creator.hpp"

namespace GFlowSimulation {

  class WallCreator : public AreaCreator  {
  public:
    //! \brief Create some particles and objects from part of a parse tree.
    //!
    //! \param head The head node of the parse tree.
    //! \param gflow The gflow object.
    //! \param variable The map of variables.
    virtual void createArea(HeadNode*, GFlow*, const std::map<string, string>&, vector<ParticleFixer>&) override;
  };

}
#endif // __WALL_CREATOR_HPP__GFLOW__