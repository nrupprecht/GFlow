#ifndef __GRAPH_OBJECT_HPP__GFLOW__
#define __GRAPH_OBJECT_HPP__GFLOW__

#include "../../base/dataobject.hpp"

namespace GFlowSimulation {

  /**
  *  \brief The parent class for data objects whose output is a csv file of (x, f(x)) data
  *  that can be printed in a regular, 2d graph.
  *
  */
  class GraphObject : public DataObject {
  public:
    //! \brief Default constructor.
    GraphObject(GFlow*, const string&);

    //! \brief Axis name setting constructor.
    GraphObject(GFlow*, const string&, const string&, const string&);

    //! \brief Clear the data.
    virtual void pre_integrate() override;

    //! \brief Collect the position data from simdata --- happens during the post-step phase.
    virtual void post_step() override = 0;

    //! \brief Write data to a file - if true, the string is a path, and you should use your own name as the file name.
    //!
    //! Returns true for success.
    virtual bool writeToFile(string, bool=true) override;

  protected:

    //! \brief The data, as a vector of pairs of real numbers
    vector<RPair> data;

    //! \brief The axis labels.
    string axis_x="x", axis_y="y";

    //! \brief If true, we use vistools to print the graph.
    bool print_plot = true;
  };

}
#endif // __GRAPH_OBJECT_HPP__GFLOW__