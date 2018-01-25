/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 13, 2017
 *
 */

#ifndef __FORCE_HANDLER_HPP__
#define __FORCE_HANDLER_HPP__

// Includes
#include "SimData.hpp"
// #include "../forces/InteractionFunctions.hpp"
// #include "../forces/WallInteractionFunctions.hpp"

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
    void pForcesRec(const VListType&, SimData*, vector<PData>&) const;
    
    // Do particle-wall forces
    void wForces(const WListType&, SimData*) const;
    void wForcesRec(const WListType&, SimData*, vector<PData>&) const;

  private:
    inline void interactP(int, int, SimData*, RealType&, RealType&, bool = true) const;
    inline void interactW(int, int, SimData*, RealType&, RealType&, bool = true) const;

    //InteractionFunction interactionFunctions[16];
    //WallInteractionFunction wallInteractionFunctions[4];
  };

}
#endif // __FORCE_HANDLER_HPP__
