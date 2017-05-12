#ifndef __CREATOR_HPP__
#define __CREATOR_HPP__

#include "SimData.hpp"
#include "../objects/Particle.hpp"
#include "../objects/Wall.hpp"

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
    
  private:
    
  };

}
#endif // __CREATOR_HPP__
