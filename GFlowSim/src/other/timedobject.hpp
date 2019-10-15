#ifndef __TIMED_OBJECT_HPP__GFLOW__
#define __TIMED_OBJECT_HPP__GFLOW__

#include "../utility/timer.hpp"

namespace GFlowSimulation {

  class TimedObject {
  public:
    //! \brief Start the object timer.
    void start_timer();

    //! \brief Stop the object timer.
    void stop_timer();

    //! \brief Clear the object timer.
    void clear_timer();

    //! \brief Get the time of the object timer.
    double get_time();

    //! \brief Check whether the timers are all on or off.
    static bool getTimingOn();

    //! \brief Turn on or off all the timers.
    static void setTimingOn(bool t);

  protected:
    //! \brief The actual timer.
    Timer timer;

    //! \brief Whether the timers should work.
    //!
    //! Not doing timing can actually save a lot of time depending on the system.
    static bool timing_on;
  };

}

#endif // __TIMED_OBJECT_HPP__GFLOW__
