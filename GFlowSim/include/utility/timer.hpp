#ifndef __TIMER_HPP__GFLOW__
#define __TIMER_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  /*
  * \brief The Timer class acts as a stopwatch, allowing high resolution timing.
  *
  * Timer is pretty much a wrapper for a couple std::chrono features. You can
  * use the timer to start and stop times (like a stopwatch), get the current elapsed
  * time (without stopping the stopwatch), get the last amount of time that elapsed 
  * between a start() and stop(), get the total elapsed time, and clear the timer.
  *
  */
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

    //! \brief Get the last amount of time that the timer recorded.
    double last();

    //! \brief Get current elapsed time (without stopping the stopwatch).
    double current();

    //! \brief Return whether the timer is running.
    bool running();
    
  private:
    //! \brief The start time.
    high_resolution_clock::time_point _s;
    //! \brief The total time recorded by the timer.
    double _t = 0;
    //! \brief The last time interval recorded.
    double _last_t = 0;
    //! \brief Flag that is true when the timer has been started, but not stopped.
    bool _running = false;
  };

}
#endif // __TIMER_HPP__GFLOW__