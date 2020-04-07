#ifndef __COURSE_POSITION_DATA_HPP__GFLOW__
#define __COURSE_POSITION_DATA_HPP__GFLOW__

#include "../base/dataobject.hpp"
#include "dataobjecttypes/data-base/particle-store-data.hpp"

namespace GFlowSimulation {

  class CoursePositionData : public DataObject, public ParticleStoreData {
  public:
    CoursePositionData(GFlow*);

    //! \brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step() override;

  private:

    //! \brief The time steps of when the data was gathered
    vector<float> timeStamps;

    vector<vector<real> > grid_data;

    //! \brief The maximum dimension of the grid. This parameter is used to choose the grid dimensions.
    int max_grid_dimension = 512;

    //! \brief The width of the grid.
    int width = 0;
    //! \brief The height of the grid.
    int height = 0;

  };

}
#endif // __COURSE_POSITION_DATA_HPP__GFLOW__