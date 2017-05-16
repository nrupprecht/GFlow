/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 15, 2017
 *
 */

#ifndef __EXTERNAL_FORCE_HPP__
#define __EXTERNAL_FORCE_HPP__

// Includes
#include "../control/SimData.hpp"

namespace GFlow {

  /*
   * @class ExternalForce
   * Base class for external forces (e.g. gravity, drag, etc)
   */
  class ExternalForce {
  public:
    void applyForce(SimData* simData, int i) {
      _applyForce(simData, i);
    }

  protected:
    // Private virtual functions, purely abstract
    virtual void _applyForce(SimData*, int) const = 0;
  };

}
#endif // __EXTERNAL_FORCE_HPP__
