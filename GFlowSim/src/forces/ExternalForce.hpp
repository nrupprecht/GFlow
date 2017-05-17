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
    void applyForce(SimData* simData) {
      _applyForce(simData);
    }

  protected:
    // Private virtual functions, purely abstract
    virtual void _applyForce(SimData*) const = 0;
  };

}
#endif // __EXTERNAL_FORCE_HPP__
