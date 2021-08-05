#ifndef __MULTIGRAPH_HPP__GFLOW__
#define __MULTIGRAPH_HPP__GFLOW__

#include "../../base/dataobject.hpp"
#include "data-base/multigraphdata.hpp"

namespace GFlowSimulation {

  class MultiGraphObject : public DataObject, public MultiGraphData {
  public:
    //! \brief Default constructor. Takes the number of entries.
    MultiGraphObject(GFlow*, const string&, int);

    //! \brief Axis name setting constructor.
    MultiGraphObject(GFlow*, const string&, const string&, const string&, int);

    //! \brief Clear the data.
    virtual void pre_integrate() override;

    //! \brief Write data to a file - if true, the string is a path, and you should use your own name as the file name.
    //!
    //! Returns true for success.
    virtual bool writeToFile(string, bool=true) override;

  };

}
#endif // __MULTIGRAPH_HPP__GFLOW__
