#ifndef __TIMER_HPP__GFLOW__
#define __TIMER_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  class Timer {
  public:
    //! \brief Start the stopwatch.
    void start();

    //! \brief Stop the stopwatch.
    void stop();

    //! \brief Clear the stopwatch.
    void clear();

    //! \brief Elapsed time of the stopwatch.
    double time();

    //! \brief Get current elapsed time (without stopping the stopwatch).
    double current();
    
  private:
    high_resolution_clock::time_point _s;
    double _t = 0;
  };

}
#endif // __TIMER_HPP__GFLOW__