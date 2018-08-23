#ifndef __PALETTE_HPP__GFLOW__
#define __PALETTE_HPP__GFLOW__

#include "../EasyBMP/EasyBMP.h"
#include "../utility/utility.hpp"
#include <functional>

namespace GFlowSimulation {

  // Some common colors
  const RGBApixel RGB_Black(0,0,0);
  const RGBApixel RGB_White(255,255,255);
  const RGBApixel RGB_Red(255,0,0);
  const RGBApixel RGB_Green(0,255,0);
  const RGBApixel RGB_Blue(0,0,255);

  //! @brief Function that creates a random color
  inline RGBApixel randomColor() {
    int r = static_cast<int>(drand48()*255);
    int g = static_cast<int>(drand48()*255);
    int b = static_cast<int>(drand48()*255);
    return RGBApixel(r, g, b);
  }

  //! @brief Create a color with the same strength of R,G,B
  inline RGBApixel makePixel(int s) {
    return RGBApixel(s, s, s);
  }

  template<typename T> inline T clamp(T x) {
    return x<T(0) ? T(0) : x;
  }

  // @brief Returns a color based on a two number parameterization.
  typedef std::function<RGBApixel(float, float, bool&)> ColorFunction;

  //! @brief Color something white.
  inline RGBApixel colorWhite(float, float, bool doColor) {
    doColor = true;
    return RGB_White;
  }

  // @brief Dot product of colors.
  inline RGBApixel operator*(const RGBApixel a, const RGBApixel b) {
    return RGBApixel(a.Red*b.Red/255, a.Green*b.Green/255, a.Blue*b.Blue/255);
  }

  //! @brief Scale a color by a float.
  inline RGBApixel operator*(float a, const RGBApixel c) {
    return RGBApixel(
      static_cast<ebmpBYTE>(a*c.Red), 
      static_cast<ebmpBYTE>(a*c.Green), 
      static_cast<ebmpBYTE>(a*c.Blue)
    );
  }

  struct PixL {
    //! @brief Default constructor.
    PixL() {};

    //! @brief The colors that have been set for this PixL
    vector<RGBApixel> colors;
    //! @brief The depths of the colors.
    vector<float> depths;
    //! @brief The object that set the color.
    vector<int> objects;
  };

  struct Vec3 {
    Vec3() : x(0), y(0), z(0) {};
    Vec3(RealType _x, RealType _y, RealType _z) : x(_x), y(_y), z(_z) {};

    Vec3 operator+(const Vec3 v) const {
      return Vec3(x+v.x, y+v.y, z+v.z);
    }

    Vec3 operator*(RealType s) const {
      return Vec3(s*x, s*y, s*z);
    }

    RealType x, y, z;
  };

  struct GraphOptions {
    GraphOptions() : useMaxX(false), useMinX(false), useMaxY(false), useMinY(false) {};
    bool useMaxX, useMinX, useMaxY, useMinY;
    float minX, maxX, minY, maxY;
  };

  inline RGBApixel colorAngle(float theta) {
    const static Vec3 C(1./3., 1./3., 1./3.);
    const static Vec3 N1(0.4*sqrt(2./3.), -0.2*sqrt(2./3.), 0.2*sqrt(2./3.));
    const static Vec3 N2(0., 0.2*sqrt(2.), -0.2*sqrt(2.));
    // Calculate the color
    Vec3 cVec = C + N1*cos(theta) + N2*sin(theta);
    return RGBApixel( 
      static_cast<int>(255*1.5*cVec.x), 
      static_cast<int>(255*1.5*cVec.y), 
      static_cast<int>(255*1.5*cVec.z)
    );
  }

  /**
  *  @brief The Palette class, used to draw images
  *
  */
  class Palette {
  public:
    //! @brief Constructor.
    Palette(int, int);

    //! @brief Copy constructor.
    Palette(const Palette&);

    //! @brief Move constructor.
    Palette(const Palette&&);

    //! @brief Destructor.
    ~Palette();

    //! @brief Move equals operator.
    Palette& operator=(const Palette&&);

    //! @brief Draw the pallate to a bmp image.
    void writeToFile(string);

    //! @brief Get a subset of the palette to draw on, by pixel coordinates.
    Palette getSubPalette(int, int, int, int);

    //! @brief Get a subset of the palette to draw on, by scaled coordinates.
    Palette getSubPalette(float, float, float, float);

    // --- Drawing functions

    //! @brief Draw a circle by giving scaled coordinates and radius.
    void drawCircleByFactors(float, float, float, ColorFunction, bool=true);

    // @brief Draw a line by giving scaled coordinates.
    void drawLineByFactors(float, float, float, float, ColorFunction, float);

    //! @brief Cover the entire palette with a single color.
    void coverPalette(RGBApixel);

    //! @brief Draw a 2d plot
    void drawGraph2d(vector<pair<float, float> >&, GraphOptions&);

    // --- Accessors

    //! @brief Get the width (in pixels) of the area the palette covers.
    int getWidth() const;

    //! @brief Get the height (in pixels) of the area the palette covers.
    int getHeight() const;

    // --- Exception class

    class PaletteOutOfBounds {};

  private:
    //! @brief Private constructor for use in making subpalettes
    Palette(int, int, int, int, BMP*, int*, int*);

    //! @brief Gives the factor of a pixel
    inline std::pair<float, float> pixelFactor(int, int) const;

    inline float pixelFactorX(int) const;
    inline float pixelFactorY(int) const;

    //! @brief Set a pixel relative to the owned coordinates (origin is {owned[0], owned[2]})
    inline void setPixel(int, int, RGBApixel);

    //! @brief Draw a line using Bresenham's algorithm.
    inline void drawLine_BresenhamAlgorithm(int, int, int, int, RGBApixel);

    //! @brief Draw a line using Wu's algorithm.
    inline void drawLine_WuAlgorithm(int, int, int, int, RGBApixel);

    //! @brief Draw a circle using the midpoint algorithm
    inline void drawCircle_MidpointAlgorithm(int, int, int, RGBApixel);

    //! @brief The image data.
    BMP *image;

    //! @brief The number of references to the image data.
    int *refs;

    //! @brief Bounds - left, right, bottom, top. 
    //! I guess left and bottom are always 0.
    int bounds[4];

    //! @brief Owned bounds - the subset of the pixels this palette owns.
    int owned[4];

    //! TEST
    PixL *pixelData;

  };


}
#endif // __PALETTE_HPP__GFLOW__