#ifndef __SIM_DATA_HANDLER_HPP__GFLOW__
#define __SIM_DATA_HANDLER_HPP__GFLOW__

#include "simdata.hpp"

namespace GFlowSim {

  class SimDataHandler {
  public:
    SimDataHandler(GFlow*);

    //! @brief Set all data entries of the type specified by the string to a default value.
    void setAllDataF(string, RealType);

    //! @brief Set all data entries of the type specified by the string to a default value.
    void setAllDataI(string, int);

  private:
    
  };

}
#endif // __SIM_DATA_HANDLER_HPP__GFLOW__