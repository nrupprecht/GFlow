/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#ifndef __SIM_DATA_HPP__
#define __SIM_DATA_HPP__

// Includes
#include "../../include/Utility.hpp"
#include "../../include/aligned_array.hpp"

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

    // Particle data
    aligned_array<RealType> px;
    aligned_array<RealType> py;
    aligned_array<RealType> vx;
    aligned_array<RealType> vy;
    aligned_array<RealType> fx;
    aligned_array<RealType> fy;
    aligned_array<RealType> th;
    aligned_array<RealType> om;
    aligned_array<RealType> tq;
    // ... Other data ...

    // How long the simulation should run for. --* This probably doesn't belong here *--
    RealType runTime;
    
  };

}
#endif // __SIM_DATA_HPP__
