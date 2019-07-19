#ifndef __HARD_SPHERE_HPP__GFLOW__
#define __HARD_SPHERE_HPP__GFLOW__

#include "../base/interaction.hpp"

namespace GFlowSimulation {

  /**
  *  \brief Parent class for hard sphere interactions.
  */
  class HardSphere : public Interaction {
  public:
    //! \brief Default constructor.
    HardSphere(GFlow *gflow) : Interaction(gflow), repulsion(DEFAULT_HARD_SPHERE_REPULSION) {};
  
    //! \brief Set the repulsion parameter.
    void setRepulsion(RealType r) { repulsion = r; }

    //! \brief Suggests a safe timescale given the minimum mass of a particle that has this interaction.
    RealType suggest_timescale(RealType mass) const override {
      return 2*PI/sqrt(2*repulsion/mass);
    }

  protected:
    //! \brief The hard sphere repulsion parameter.
    RealType repulsion;
  };

}

#endif // __HARD_SPHERE_HPP__GFLOW__