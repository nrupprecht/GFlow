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
    HardSphere(GFlow *gflow, InteractionHandler *hndlr) 
      : Interaction(gflow, hndlr), repulsion(DEFAULT_HARD_SPHERE_REPULSION) {};
  
    //! \brief Set the repulsion parameter.
    void setRepulsion(RealType r) { repulsion = r; }

  protected:
    //! \brief The hard sphere repulsion parameter.
    RealType repulsion;
  };

}

#endif // __HARD_SPHERE_HPP__GFLOW__