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
    ExternalForce() : finished(false) {};
    virtual ~ExternalForce() {};

    void applyForce(SimData* simData) {
      _applyForce(simData);
    }

    string summary() const {
      return _summary();
    }

    // Get the finished flag
    bool getFinished() { return finished; }

    // Set the finished flag
    void setFinished(bool f) { finished = f; }

    // Reset the force (if applicable)
    virtual void reset() {};

  protected:
    // Private virtual functions, purely abstract
    virtual void _applyForce(SimData*) const = 0;
    virtual string _summary()          const = 0;

    // Whether the force is done and can be removed
    bool finished;
  };

}
#endif // __EXTERNAL_FORCE_HPP__
