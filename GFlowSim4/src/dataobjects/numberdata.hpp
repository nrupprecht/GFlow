#ifndef __NUMBER_DATA_HPP__GFLOW__
#define __NUMBER_DATA_HPP__GFLOW__

#include "../base/dataobject.hpp"

namespace GFlowSimulation {

  class NumberData : public DataObject {
  public:
    NumberData(GFlow*);

    //! @brief Collect the position data from simdata --- happens during the post-step phase.
    virtual void post_step() override;

    //! @brief Write data to a file - if true, the string is a path, and you should use your own name as the file name.
    // Returns true for success.
    virtual bool writeToFile(string, bool=true) override;

  private:
    vector< vector<RPair> > numberData;
  };

}
#endif // __NUMBER_DATA_HPP__GFLOW__