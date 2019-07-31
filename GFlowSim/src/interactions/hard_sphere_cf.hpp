#ifndef __HARD_SPHERE_CF_HPP__GFLOW__
#define __HARD_SPHERE_CF_HPP__GFLOW__

#include "hard_sphere.hpp"

namespace GFlowSimulation {

  /** 
  *  \brief Parent class for dissipative hard sphere interactions.
  */
  class HardSphereCf : public HardSphere {
  public:
    //! \brief Default constructor.
    HardSphereCf(GFlow *gflow) 
      : HardSphere(gflow), dissipation(DEFAULT_HARD_SPHERE_DISSIPATION), angular_dissipation(DEFAULT_HARD_SPHERE_DISSIPATION), mu(0.25) {};

    //! \brief Set the dissipation constant.
    void setDissipation(RealType d) { dissipation = d>=0 ? d : dissipation; }

    //! \brief Set the dissipation constant.
    void setAngularDissipation(RealType d) { angular_dissipation = d>=0 ? d : angular_dissipation; }

    //! \brief Set the coefficient of friction.
    void setMu(RealType c) { mu = c>=0 ? c : mu; }

    //! \brief Suggests a safe timescale given the minimum mass of a particle that has this interaction.
    RealType suggest_timescale(RealType mass) const override {
      return 2*PI/sqrt(2*repulsion/mass);
    }

  protected:
    //! \brief The dissipation constant.
    RealType dissipation;

    //! \brief The angular dissipation constant.
    RealType angular_dissipation;

    //! \brief The coefficient of friction.
    RealType mu;
  };

}

#endif // __HARD_SPHERE_DS_HPP__GFLOW__