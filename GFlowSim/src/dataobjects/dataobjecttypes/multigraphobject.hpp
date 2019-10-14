#ifndef __MULTIGRAPH_HPP__GFLOW__
#define __MULTIGRAPH_HPP__GFLOW__

#include "../../base/dataobject.hpp"

namespace GFlowSimulation {

  class MultiGraphObject : public DataObject {
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

    //! \brief Compute averages of data.
    RealType ave(int);

    //! \brief Get a single data entry as a vector or RPairs.
    vector<RPair> getEntry(int);

    //! \brief Exception class for trying to access data when there is none.
    struct EmptyData : public Exception {
      EmptyData() : Exception() {};
      EmptyData(string mess) : Exception(mess) {};
    };

    //! \brief Exception class for access being out of bounds.
    struct MultiDataOutOfBounds : public Exception {
      MultiDataOutOfBounds() : Exception() {};
      MultiDataOutOfBounds(string mess) : Exception(mess) {};
    };

  protected:
    //! \brief Reset the data array.
    void resetData(int=0);

    //! \brief Add a set of zeros to all the data. This makes sure that all the data is of the same length.
    void addEntry();

    //! \brief Add a set of zeros to all the data, and sets the X value of the entry. This makes sure that all the data is of the same length.
    void addEntry(RealType);

    //! \brief Get the latest x value so it can be set.
    RealType& getX();

    //! \brief Get the different latest y values, so they can be set.
    RealType &getY(int);

    RealType& atX(int);

    RealType& atY(int, int);

    //! \brief Gather and average data on processor 0.
    void gatherAverageData(const RealType, const Vec, int);

    //! \brief The different types of data the multigraph is keeping track of.
    //!
    //! The first entry is the x value, the other [ndata] values are the data.
    vector<vector<RealType> > multi_data;

    //! \brief Controls which entries to write to files.
    vector<bool> write_data;

    //! \brief The number of types of data that the class is keeping track of.
    int ndata;

    //! \brief The number of data entries that were added.
    int ndata_points = 0;

    //! \brief The axis labels.
    string axis_x="x";
    vector<string> axis_y;

    //! \brief If true, we use vistools to print the graph.
    bool print_plot = true;
  };

}
#endif // __MULTIGRAPH_HPP__GFLOW__