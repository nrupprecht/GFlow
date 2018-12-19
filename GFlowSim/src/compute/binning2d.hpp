#ifndef __BINNING_2D_HPP__GFLOW__
#define __BINNING_2D_HPP__GFLOW__

#include "../utility/utility.hpp"

namespace GFlowSimulation {

  class Binning2d {
  public:
    Binning2d(int x, int y);

    ~Binning2d();

    void setBounds(const Bounds&);
    void setBounds(float, float, float, float);

    void bin_data(const vector<float>&, int, int);
    
    const vector<int>& at(int, int);

  private:
    //! @brief Number of bins in each dimension.
    int dimx, dimy;
    // @brief Dimensions of a bin.
    float dx, dy;
    
    //! @brief Width of the bounds.
    float min[2];
    float max[2];

    //! @brief Bins that hold the ids of data that falls into the bins
    vector<int> *bins;
  };

}
#endif // __BINNING_2D_HPP__GFLOW__