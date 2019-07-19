#ifndef __HARD_SPHERE_DS_HPP__GFLOW__
#define __HARD_SPHERE_DS_HPP__GFLOW__

#include "../base/interaction.hpp"

namespace GFlowSimulation {

  /** 
  *  \brief Parent class for dissipative hard sphere interactions.
  */
  class HardSphereDs : public Interaction {
  public:
    //! \brief Default constructor.
    HardSphereDs(GFlow *gflow) 
      : Interaction(gflow), repulsion(DEFAULT_HARD_SPHERE_REPULSION), dissipation(DEFAULT_HARD_SPHERE_DISSIPATION) {};

    //! \brief Set the repulsion constant.
    void setRepulsion(RealType r) { repulsion = r; }

    //! \brief Set the dissipation constant.
    void setDissipation(RealType d) { dissipation = d>0 ? d : dissipation; }

    //! \brief Suggests a safe timescale given the minimum mass of a particle that has this interaction.
    RealType suggest_timescale(RealType mass) const override {
      return 2*PI/sqrt(2*repulsion/mass);
    }

  protected:
    //! \brief The repulsion constant for the spheres.
    RealType repulsion;

    //! \brief The dissipation constant.
    RealType dissipation;
  };

}

#endif // __HARD_SPHERE_DS_HPP__GFLOW__