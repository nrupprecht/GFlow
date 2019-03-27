#ifndef __BOUNDARY_FORCE_DATA_HPP__GFLOW__
#define __BOUNDARY_FORCE_DATA_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"

namespace GFlowSimulation {

  class BoundaryForceData : public GraphObject {
  public:
    //! \brief Default constructor.
    BoundaryForceData(GFlow*);

    //! \brief Collect boundary force data.
    virtual void post_step() override;

    //! \brief Get the average of the boundary forces, over time.
    RealType getAverage() const;
  };

}
#endif // __BOUNDARY_FORCE_DATA_HPP__GFLOW__