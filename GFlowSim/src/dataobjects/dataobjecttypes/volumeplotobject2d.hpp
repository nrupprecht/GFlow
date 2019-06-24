#ifndef __VOLUME_PLOT_OBJECT_2D_HPP__GFLOW__
#define __VOLUME_PLOT_OBJECT_2D_HPP__GFLOW__

#include "../../base/dataobject.hpp"

namespace GFlowSimulation {

  class VolumePlotObject2D : public DataObject {
  public:
    VolumePlotObject2D(GFlow*, const string&);

  private:
    //! \brief Binning of the data.
    vector<vector<RealType> > binning;
    //! \brief Binning for recording counts.
    vector<vector<int> > counts;

    //! \brief The discritization in the X direction.
    int binX = 0;
    //! \brief The discritization in the Y direction.
    int binY = 0;

  };

}
#endif // __VOLUME_PLOT_OBJECT_2D_HPP__GFLOW__