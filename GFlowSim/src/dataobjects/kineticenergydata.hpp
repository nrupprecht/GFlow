#ifndef __KINETIC_ENERGY_DATA_HPP__GFLOW__
#define __KINETIC_ENERGY_DATA_HPP__GFLOW__

#include "../base/dataobject.hpp"

namespace GFlowSimulation {

  class KineticEnergyData : public DataObject {
  public:
    //! @brief Constructor
    KineticEnergyData(GFlow*, bool=true);

    //! @brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    //! @brief Write data to a file - if true, the string is a path, and you should use your own name as the file name.
    //!
    //! Returns true for success
    virtual bool writeToFile(string, bool=true);

  private:
    //! @brief The data
    vector<RPair> keData;
    //! @brief Use ave
    bool useAve;
  };

}
#endif // __KINETIC_ENERGY_DATA_HPP__GFLOW__