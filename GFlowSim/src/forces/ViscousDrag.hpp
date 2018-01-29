/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 15, 2017
 *
 */

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

    // Get the viscosity
    RealType getViscosity() const { return viscosity; }

  protected:
    // Inherited private virtual functions
    virtual void _applyForce(SimData*) const;
    virtual string _summary() const;

    // The viscosity
    RealType viscosity;
  };

}
#endif // __VISCOUS_DRAG_HPP__
