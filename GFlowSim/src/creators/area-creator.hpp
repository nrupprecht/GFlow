#ifndef __AREA_CREATOR__HPP__GFLOW__
#define __AREA_CREATOR__HPP__GFLOW__

#include "../base/creator.hpp"
#include "../utility/parsehelper.hpp"
#include "particletemplate.hpp"
#include "fillbounds.hpp"

namespace GFlowSimulation {

  class AreaCreator {
  public:
    //! \brief Create some particles and objects from part of a parse tree.
    //! \param head The head node of the parse tree.
    //! \param gflow The gflow object.
    //! \param variable The map of variables.
    virtual void createArea(HeadNode*, GFlow*, std::map<string, string>&)=0;
  };

}
#endif // __AREA_CREATOR__HPP__GFLOW__