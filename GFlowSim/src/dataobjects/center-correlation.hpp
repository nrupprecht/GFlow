#ifndef __CENTER_CORRELATION_HPP__GFLOW__
#define __CENTER_CORRELATION_HPP__GFLOW__

#include "../base/dataobject.hpp"

namespace GFlowSimulation {

  class CenterCorrelation : public DataObject {
  public:
    //! \brief Constructor
    CenterCorrelation(GFlow*);

    //! \brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    //! \brief Write data to a file - if true, the string is a path, and you should use your own name as the file name.
    //!
    //! Returns true for success
    virtual bool writeToFile(string, bool=true);

  private:
    //! \brief Bins
    vector<int> bins;

    //! \brief The radius of the cylinder of interest.
    RealType radius;

    //! \brief The number of desired bins
    int nbins;
  };

}
#endif // __CENTER_CORRELATION_HPP__GFLOW__