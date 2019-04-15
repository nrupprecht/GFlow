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

  protected:
    //! \brief The actual timer.
    Timer timer;
  };

}

#endif // __TIMED_OBJECT_HPP__GFLOW__