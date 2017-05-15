#ifndef __VISCOUS_DRAG_HPP__
#define __VISCOUS_DRAG_HPP__

// Includes
#include "DragForce.hpp"

namespace GFlow {
  
  /*
   * @class ViscousForce
   * A viscous force
   *
   */
  class ViscousDrag : public DragForce {
  public:
    // Default constructor - sets the viscosity to that of water
    ViscousDrag();

    // Viscosity setting constructor
    ViscousDrag(RealType);

  protected:
    // Inherited private virtual functions
    virtual void _applyForce(SimData*, int) const;

    RealType viscosity;
  };

}
#endif // __VISCOUS_DRAG_HPP__
