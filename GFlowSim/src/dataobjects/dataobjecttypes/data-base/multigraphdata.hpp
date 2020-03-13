#ifndef __MULTI_GRAPH_DATA_HPP__GFLOW__
#define __MULTI_GRAPH_DATA_HPP__GFLOW__

#include "../../../gflow.hpp"

namespace GFlowSimulation {

  class MultiGraphData {
  public:
    //! \brief Default constructor. Takes the number of entries.
    MultiGraphData(int);

    //! \brief Axis name setting constructor.
    MultiGraphData(const string&, const string&, int);

    //! \brief Write data to a file - if true, the string is a path, and you should use your own name as the file name.
    //!
    //! Returns true for success.
    virtual bool write_to_file(const string&, const string&);

    //! \brief Compute averages of data.
    real ave(int);

    //! \brief Get a single data entry as a vector or RPairs.
    vector<RPair> getEntry(int) const;

    //! \brief Get the number of datapoints entered into the multigraph object.
    int size() const;

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

    virtual string correct_dir_name(const string& name) { return name; }

    //! \brief Reset the data array.
    void resetData(int=0);

    //! \brief Add a set of zeros to all the data. This makes sure that all the data is of the same length.
    void addEntry();

    //! \brief Add a set of zeros to all the data, and sets the X value of the entry. This makes sure that all the data is of the same length.
    void addEntry(real);

    //! \brief Get the latest x value so it can be set.
    real& getX();

    //! \brief Get the different latest y values, so they can be set.
    real &getY(int);

    real& atX(int);

    real& atY(int, int);

    //! \brief Gather and average data on processor 0.
    void gatherAverageData(const real, const Vec, int);

    //! \brief Gather sum of all data on processor 0.
    void gatherData(const real, const Vec);

    //! \brief The different types of data the multigraph is keeping track of.
    //!
    //! The first entry is the x value, the other [ndata] values are the data.
    vector<vector<real> > multi_data;

    //! \brief Controls which entries to write to files.
    vector<bool> write_data;

    //! \brief The number of types of data that the class is keeping track of.
    int ndata;

    //! \brief The axis labels.
    string axis_x="x";
    vector<string> axis_y;
  };

}
#endif // __MULTI_GRAPH_DATA_HPP__GFLOW__