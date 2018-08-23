#ifndef __VISUALIZATION_HPP__GFLOW__
#define __VISUALIZATION_HPP__GFLOW__

#include "palette.hpp"

namespace GFlowSimulation {

  /**
  *  @brief Creates visualization from
  *
  */
  class Visualization {
  public:
    //! @brief Constructor.
    Visualization();

    //! @brief Destructor.
    ~Visualization();

    //! @brief Create a directory filled with BMP renderings of the system.
    //!
    //! This can be used to create a movie.
    void createBMPs(string, const vector<RealType*>&, const vector<int>&, int, Bounds&, int) const;

  private:
    //! @brief Creates a single frame.
    inline void createImage(string, RealType*, int, int, Bounds&, int) const;

    inline void findMaxVSqr(const vector<RealType*>&, const vector<int>&, int) const;

    //! @brief Where the position data starts.
    int pos_place;

    //! @brief Where the velocity data starts.
    int vel_place;

    //! @brief Where in the data for a particle is the radius.
    int sg_place;

    //! @brief Where in the data for a particle is its type.
    int type_place;

    //! @brief The dimensions of the image (it will be the same in x and y)
    int resolution;

    //! @brief Whether to wrap at the boundaries or not
    bool do_wrap;

    mutable RealType maxVsqr;

    RGBApixel background;

    unsigned int color_option;

    RGBApixel *colorBank;
  };

}
#endif // __VISUALIZATION_HPP__GFLOW__