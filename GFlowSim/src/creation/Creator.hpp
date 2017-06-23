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
#include "../control/StatusObject.hpp"
#include "../integrators/VelocityVerletIntegrator.hpp"
#include "../objects/Particle.hpp"
#include "../objects/Wall.hpp"
#include "../forces/ViscousDrag.hpp"
#include "../creation/DataCreator.hpp"

namespace GFlow {
  
  // Forward reference to FileParser
  class FileParser;

  /*
   * @class Creator
   * Creates simulation data
   *
   */
  class Creator {
  public:
    
    // Modify the simdata for a buoyancy study
    void createBuoyancy(SimData *&, Integrator *&, RealType, RealType, vec2, bool=false, bool=false);

    // Modify the simdata for an aerodynamics study
    void createAero(SimData *&, Integrator *&, RealType, RealType, vec2, bool);

    // FileParser is a friend class
    friend class FileParser;
  private:
    bool createRegion(Region&, SimData*);    
  };

}
#endif // __CREATOR_HPP__
