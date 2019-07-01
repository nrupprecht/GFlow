#ifndef __VOLUME_PLOT_OBJECT_2D_HPP__GFLOW__
#define __VOLUME_PLOT_OBJECT_2D_HPP__GFLOW__

#include "../../base/dataobject.hpp"

namespace GFlowSimulation {

  /***
  *  \brief Data object that makes density plots.
  *
  *  
  *
  */
  class VolumePlotObject2D : public DataObject {
  public:
    VolumePlotObject2D(GFlow*, const string&, int);

    //! \brief Clear the data.
    virtual void pre_integrate() override;

    //! \brief Standard way of writing data to files.
    virtual bool writeToFile(string, bool=true) override;

    //! \brief Set the course graining of the binning.
    void setBins(int, int);

    //! \brief Set the print counts flag.
    void setPrintCounts(bool);

    // --- Data functions ---

    //! \brief Add value by bin index.
    void addToBin(int, int, int, RealType, bool inc=true);

    //! \brief Add value by position.
    void addToBin(RealType, RealType, int, RealType, bool inc=true);

    //! \brief Add a vector to the binning by position.
    void addToBin(RealType, RealType, Vec&);

    //! \brief Increment the count of a bin.
    void incrementBin(int, int);

  protected:

    //! \brief What part of the simulation to focus on. 
    //!
    //! It could be the whole thing.
    Bounds focus_bounds;

    //! \brief The size of the bin tuples.
    int data_width = 1;

    //! \brief Binning of the data.
    vector<vector<Vec> > binning;
    //! \brief Binning for recording counts.
    vector<vector<int> > counts;

    //! \brief The names of the different entries in the binning.
    vector<string> entry_names;

    //! \brief The discritization in the X direction.
    int binX = 0;
    //! \brief The discritization in the Y direction.
    int binY = 0;

    //! \brief The default max number of bins in the direction of greatest size.
    //!
    //! This helps the object choose default values of binX, binY.
    int default_max_bins = 100;

    //! \brief Number of time steps on which a recording happened.
    int recorded_frames = 0;

    //! \brief Inverse grid spacings.
    RealType invdx = 0;
    RealType invdy = 0;

    //! \brief Whether to print the counts.
    bool print_counts = true;

  };

}
#endif // __VOLUME_PLOT_OBJECT_2D_HPP__GFLOW__