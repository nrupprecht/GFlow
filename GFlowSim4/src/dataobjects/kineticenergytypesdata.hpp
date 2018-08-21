#ifndef __KINETIC_ENERGY_TYPES_DATA_HPP__GFLOW__
#define __KINETIC_ENERGY_TYPES_DATA_HPP__GFLOW__

#include "../base/dataobject.hpp"

namespace GFlowSimulation {

  class KineticEnergyTypesData : public DataObject {
  public:
    // Constructor
    KineticEnergyTypesData(GFlow*, bool=true);

    // Constructor
    KineticEnergyTypesData(GFlow*, int, bool=true);

    // Destructor
    ~KineticEnergyTypesData();

    // Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    // Write data to a file - if true, the string is a path, and you should use your own name as the file name
    // Returns true for success
    virtual bool writeToFile(string, bool=true);

  private:
    // The data
    vector<RPair> *keData;
    // Number of particle types to track
    int ntypes;
    // Use ave
    bool useAve;
  };

}
#endif // __KINETIC_ENERGY_TYPES_DATA_HPP__GFLOW__