/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 15, 2017
 *
 */

#ifndef __DRAG_FORCE_HPP__
#define __DRAG_FORCE_HPP__

// Includes
#include "../control/SimData.hpp"

namespace GFlow {

  /*
   * @class DragForce
   * Base class for drag forces
   */
  class DragForce {
  public:
    void applyForce(SimData* simData, int i) {
      _applyForce(simData, i);
    }

  protected:
    // Private virtual functions, purely abstract
    virtual void _applyForce(SimData*, int) const = 0;
  };

}
#endif // __DRAG_FORCE_HPP__
