#ifndef __KINETIC_ENERGY_DATA_HPP__GFLOW__
#define __KINETIC_ENERGY_DATA_HPP__GFLOW__

#include "dataobject.hpp"

namespace GFlowSimulation {

  class KineticEnergyData : public DataObject {
  public:
    // Constructor
    KineticEnergyData(GFlow*, bool=true);

    // Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    // Write data to a file - if true, the string is a path, and you should use your own name as the file name
    // Returns true for success
    virtual bool writeToFile(string, bool=true);

  private:
    // The data
    vector<RPair> keData;
    // Use ave
    bool useAve;
  };

}
#endif // __KINETIC_ENERGY_DATA_HPP__GFLOW__