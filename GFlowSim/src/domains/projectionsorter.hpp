#ifndef __PROJECTION_SORTER_HPP__GFLOW__
#define __PROJECTION_SORTER_HPP__GFLOW__

#include "../base/interactionhandler.hpp"

namespace GFlowSimulation {

  class ProjectionSorter : public InteractionHandler {
  public:
    ProjectionSorter(GFlow*);

    //! \brief Construct verlet lists.
    virtual void construct() override;

    virtual void getAllWithin(int, vector<int>&, RealType=-1.) override;

    //! \brief Get all the particles withing a radius of some position.
    //!
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(Vec, vector<int>&, RealType=-1.) override;

    //! \brief Remove particles that are overlapping by more than some fraction.
    virtual void removeOverlapping(RealType) override;

  private:
    //! \brief The axis that particles are projected onto.
    Vec axis;
  };

}
#endif // __PROJECTION_SORTER_HPP__GFLOW__