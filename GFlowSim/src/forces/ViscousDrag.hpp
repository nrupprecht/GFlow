#ifndef __VISCOUS_DRAG_HPP__
#define __VISCOUS_DRAG_HPP__

// Includes
#include "ExternalForce.hpp"

namespace GFlow {
  
  /*
   * @class ViscousForce
   * A viscous force
   *
   */
  class ViscousDrag : public ExternalForce {
  public:
    // Default constructor - sets the viscosity to that of water
    ViscousDrag();

    // Viscosity setting constructor
    ViscousDrag(RealType);

  protected:
    // Inherited private virtual functions
    virtual void _applyForce(SimData*) const;

    // The viscosity
    RealType viscosity;
  };

}
#endif // __VISCOUS_DRAG_HPP__
