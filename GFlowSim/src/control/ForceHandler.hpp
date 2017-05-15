/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 13, 2017
 *
 */

#ifndef __FORCE_HANDLER_HPP__
#define __FORCE_HANDLER_HPP__

// Includes
#include "SimData.hpp"
#include "InteractionFunctions.hpp"
#include "WallInteractionFunctions.hpp"

namespace GFlow {

  /*
   * @class ForceHandler
   * Class for handling forces, base class for all force handlers
   */
  class ForceHandler {
  public:
    // Default constructor
    ForceHandler();

    // Do inter-particle forces
    void pForces(const VListType&, SimData*) const;
    
    // Do particle-wall forces
    void wForces(const WListType&, SimData*) const;

  private:
    inline void interactP(int, int, SimData*) const;
    inline void interactW(int, int, SimData*) const;

    InteractionFunction interactionFunctions[16];
    WallInteractionFunction wallInteractionFunctions[4];
  };

}
#endif // __FORCE_HANDLER_HPP__
