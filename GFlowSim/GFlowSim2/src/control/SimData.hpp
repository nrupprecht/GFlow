/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#ifndef __SIM_DATA_HPP__
#define __SIM_DATA_HPP__

// Includes
#include "../utility/Utility.hpp"
#include "../control/Sectorization.hpp"

namespace GFlow {

  /*
   * @class SimData
   *
   */
  class SimData {
  public:

    // Calculate forces
    void doForces();

    // Get run time
    RealType getRunTime() { return runTime; }

  private:

    // How long the simulation should run for. --* This probably doesn't belong here *--
    RealType runTime;
    
    // A sectorization
    Sectorization sectors;
  };

}
#endif // __SIM_DATA_HPP__
