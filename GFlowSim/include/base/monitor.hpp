#ifndef __MONITOR_HPP__GFLOW__
#define __MONITOR_HPP__GFLOW__

#include "../gflow.hpp"

namespace GFlowSimulation {

/*! \brief Monitors the simulation, can be used to print out data
 *         of check for errors.
 *
 *  A class of objects used to monitor_ the simulation. This could
 *  in theory all be done by a modifier, but monitors are supposed
 *  to only monitor_ the simulation, not change anything - except
 *  to exit if something catastrophic happens or some condition is
 *  fulfilled.
 *
 */
class Monitor : public Base {
 public:
  Monitor(GFlow *);
};

}
#endif // __MONITOR_HPP__GFLOW__