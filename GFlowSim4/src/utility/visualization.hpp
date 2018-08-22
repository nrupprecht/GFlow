#ifndef __VISUALIZATION_HPP__GFLOW__
#define __VISUALIZATION_HPP__GFLOW__

#include "../EasyBMP/EasyBMP.h"
#include "../utility/utility.hpp"

namespace GFlowSimulation {

  const RGBApixel Black(0,0,0);
  const RGBApixel White(255,255,255);
  const RGBApixel Red(255,0,0);
  const RGBApixel Green(0,255,0);
  const RGBApixel Blue(0,0,255);

  inline RGBApixel randomColor() {
    int r = static_cast<int>(drand48()*255);
    int g = static_cast<int>(drand48()*255);
    int b = static_cast<int>(drand48()*255);
    return RGBApixel(r, g, b);
  }

  // @brief Returns a color based on how far a point is from the center of a circle
  typedef RGBApixel (*ColorFunction) (float, float);

  /**
  *  @brief The Palette class, used to draw images
  *
  */
  class Palette {
    //! @brief Constructor.
    Palette(int, int);

    //! @brief Destructor.
    ~Palette();

    //! @brief Move equals operator.
    Palette& operator=(const Palette&&);

    //! @brief Draw the pallate to a bmp image.
    void writeToFile(string);

    //! @brief Get a subset of the palette to draw on.
    Palette getSubPalette(int, int, int, int);

    // --- Drawing functions

    void drawCircleByFactors(float, float, float, ColorFunction, bool=true);

  private:
    //! @brief Private constructor for use in making subpalettes
    Palette(int, int, int, int, BMP*, int*, int*);

    //! @brief The image data.
    BMP *image;

    //! @brief The number of references to the image data.
    int *refs;

    //! @brief Bounds - left, right, bottom, top.
    int bounds[4];

    //! @brief Owned bounds - the subset of the pixels this palette owns.
    int owned[4];

  };

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

    inline void createImage(string, RealType*, int, int, Bounds&, int) const;

    inline void circle(BMP&, int, int, int, RGBApixel, bool) const;

    //! @brief Where in the data for a particle is the radius
    int sg_place;

    //! @brief The dimensions of the image (it will be the same in x and y)
    int resolution;

    //! @brief Whether to wrap at the boundaries or not
    bool do_wrap;

    RGBApixel background;

    RGBApixel *colorBank;
  };

}
#endif // __VISUALIZATION_HPP__GFLOW__