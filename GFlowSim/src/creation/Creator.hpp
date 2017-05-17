/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#ifndef __CREATOR_HPP__
#define __CREATOR_HPP__

// Includes
#include "../../include/vec2d.hpp"
#include "../control/SimData.hpp"
#include "../integrators/VelocityVerletIntegrator.hpp"
#include "../objects/Particle.hpp"
#include "../objects/Wall.hpp"
#include "../forces/ViscousDrag.hpp"
#include "../creation/DataCreator.hpp"

namespace GFlow {
  
  /*
   * @class Creator
   * Creates simulation data
   *
   */
  class Creator {
  public:
    
    // In the future this will take argments
    SimData* create();

    bool createRegion(Region&, SimData*);
    
  private:
    
  };

}
#endif // __CREATOR_HPP__
