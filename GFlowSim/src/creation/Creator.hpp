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
// #include "../creation/DataCreator.hpp"

namespace GFlow {

  /*
   * @class Region
   * Contains the data for constructing a region full of particles
   *
   */
  /*
  struct Region {
    RealType left, right, bottom, top;
    int number;
    
    DataCreator *sigma, *position
    
  };
  */

  /*
   * @class Creator
   * Creates simulation data
   *
   */
  class Creator {
  public:
    
    // In the future this will take argments
    SimData* create();

    // void createRegion(Region, SimData*);
    
  private:
    
  };

}
#endif // __CREATOR_HPP__
