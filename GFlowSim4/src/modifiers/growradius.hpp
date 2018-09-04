#ifndef __GROW_RADIUS_HPP__GFLOW__
#define __GROW_RADIUS_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  /**
  *  @brief A modifier that grows the radius (sigma - cutoff) of a particle at a constant rate over time.
  * 
  *  Note this class acts on a single particle.
  */
  class GrowRadius : public Modifier {
  public:
    GrowRadius(GFlow*, int, RealType, RealType, RealType);

    virtual void post_forces() override;

    //! @brief An exception class indicating that a negative radius has been indicated.
    class BadRadius {};

  private:
    //! @brief The global id of the particle we are tracking
    int global_id;
    //! @brief When the grow radius started.
    RealType time0;
    //! @brief The rate of change of the radius.
    RealType rdot;
    //! @brief The initial radius.
    RealType sigma0;
    //! @brief The final radius.
    RealType sigmaf;
  };

}
#endif // __GROW_RADIUS_HPP__GFLOW__