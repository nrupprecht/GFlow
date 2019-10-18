#ifndef __MULTI_BIN_OBJECT_HPP__GFLOW__
#define __MULTI_BIN_OBJECT_HPP__GFLOW__

#include "../../base/dataobject.hpp"

namespace GFlowSimulation {

  class MultiBinObject : public DataObject {
  public:
    MultiBinObject(GFlow*, string, int);

    //! \brief Clear the data.
    virtual void pre_integrate() override;

    //! \brief Write data to a file - if true, the string is a path, and you should use your own name as the file name.
    //!
    //! Returns true for success.
    virtual bool writeToFile(string, bool=true) override;

    //! \brief Set the number of bins. This will destroy any existing data.
    void setNBins(int);

  protected:
    //! \brief Increment an entry in a bin.
    void increment(int, int, bool=true);

    //! \brief Increment the first entry in a bin.
    void increment(int, bool=true);

    //! \brief Clear out old data and create new entries.
    void resetData();

    //! \brief Create evenly spaced bins.
    void labelBins(RealType, RealType);

    //! \brief Vector in which to store x coordinates (bin labels).
    vector<RealType> labels;
    //! \brief Vector in which to store counts.
    vector<vector<long> > counts;

    //! \brief Whether to write each count entry to the file at the end of the run.
    vector<bool> write_data;

    //! \brief The number of bins.
    int nbins = 50;

    //! \brief The number of types of data that the class is keeping track of.
    int ndata;

    //! \brief The min and max label values.
    RealType min_label = 0, max_label = 0;

    //! \brief The axis labels.
    string axis_x="x";
    vector<string> axis_y;
  };

}
#endif // __MULTI_BIN_OBJECT_HPP__GFLOW__