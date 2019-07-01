#ifndef __PROJECTION_SORTER_HPP__GFLOW__
#define __PROJECTION_SORTER_HPP__GFLOW__

#include "../base/interactionhandler.hpp"

namespace GFlowSimulation {

  class ProjectionSorter : public InteractionHandler {
  public:
    ProjectionSorter(GFlow*);

    //! \brief Construct verlet lists.
    virtual void construct() override;
    
  private:
    //! \brief The axis that particles are projected onto.
    Vec axis;
  };

}
#endif // __PROJECTION_SORTER_HPP__GFLOW__