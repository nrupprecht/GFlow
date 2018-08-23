#include "palette.hpp"
// Other files
#include "vectormath.hpp"

namespace GFlowSimulation {

  Palette::Palette(int width, int height) {
    // Set the bounds
    owned[0] = bounds[0] = 0;
    owned[1] = bounds[1] = width;
    owned[2] = bounds[2] = 0;
    owned[3] = bounds[3] = height;
    // Create the image
    image = new BMP;
    image->SetSize(width, height);
    // References
    refs = new int(1);
  }
    
  Palette::~Palette() {
    // Clean up if we are the last palette to reference the object
    if (refs[0]==1) {
      --refs[0];
      delete image;
      delete refs;
    }
    else --refs[0];
  }

  Palette& Palette::operator=(const Palette&& p) {
    for (int i=0; i<4; ++i) {
      owned[i] = p.owned[i];
      bounds[i] = p.bounds[i];
    }
    // Image and refs
    image = p.image;
    refs = p.refs;
    ++refs[0];
    // Return
    return *this;
  }

  void Palette::writeToFile(string fileName) {
    image->WriteToFile(fileName.c_str());
  }

  Palette Palette::getSubPalette(int mx, int Mx, int my, int My) {
    return Palette(owned[0]+mx, owned[0]+Mx, owned[2]+my, owned[2]+My, image, refs, bounds);
  }

  Palette Palette::getSubPalette(float mfx, float Mfx, float mfy, float Mfy) {
    return Palette(owned[0]+mfx*getWidth(), owned[0]+Mfx*getWidth(), owned[2]+mfy*getHeight(), owned[2]+Mfy*getHeight(), image, refs, bounds);
  }

  void Palette::drawCircleByFactors(float xf, float yf, float rf, ColorFunction colorF, bool wrap) {
    // Width of owned region
    int wx = getWidth();
    int wy = getHeight();
    float iwx = 1.f/wx, iwy = 1.f/wy;
    // Find the center (in pixel coordinates)
    int px = xf*wx;
    int py = yf*wy;
    // Find dp's (radius in pixel coordinates)
    int dpx = ceil(rf*wx);
    int dpy = ceil(rf*wy);

    // Radius factor squared
    float rsqr = sqr(rf);
    bool doColor = true;
    // Go through pixels
    for (int dy=-dpy; dy<=dpy; ++dy) {
      for (int dx=-dpx; dx<=dpx; ++dx) {
        if (sqr(iwx*dx) + sqr(iwy*dy) <= rsqr) {
          // The potential pixel
          int w_px = px + dx, w_py = py + dy;
          // Get the pixel factor
          std::pair<float, float> pF = pixelFactor(w_px, w_py);
          float pfx = (pF.first - xf)/rf;
          float pfy = (pF.second - yf)/rf;
          // Wrap pixel
          if (wrap) {
            w_px = w_px >= wx ? w_px-wx : w_px;
            w_px = w_px <  0 ? w_px+wx : w_px;
            w_py = w_py >= wy ? w_py-wy : w_py;
            w_py = w_py <  0 ? w_py+wy : w_py;
            // Set pixel
            RGBApixel color = colorF(pfx, pfy, doColor);
            if (doColor) setPixel(w_px, w_py, color);
          }
          else if (0<=w_px && w_px<wx && 0<=w_py && w_py<wy) {
            // Set pixel
            RGBApixel color = colorF(pfx, pfy, doColor);
            if (doColor) setPixel(w_px, w_py, color);
          }
        }
      }
    }
  }

  void Palette::coverPalette(RGBApixel color) {
    for (int y=owned[2]; y<owned[3]; ++y)
      for (int x=owned[0]; x<owned[1]; ++x) 
        image->SetPixel(x, y, color);
  }

  int Palette::getWidth() const {
    return owned[1] - owned[0];
  }

  int Palette::getHeight() const {
    return owned[3] - owned[2];
  }

  Palette::Palette(int mx, int Mx, int my, int My, BMP *img, int *rfs, int *bnds) {
    // Set sub-bounds
    owned[0] = mx;
    owned[1] = Mx;
    owned[2] = my;
    owned[3] = My;
    // Copy full bounds
    for (int i=0; i<4; ++i) bounds[i] = bnds[i];
    // Image and refs
    image = img;
    refs = rfs;
    ++refs[0];
  }

  std::pair<float, float> Palette::pixelFactor(int x, int y) const {
    // Create a pair to set and return
    std::pair<float, float> pf;
    // We need to find where the centers of pixels are
    int X = owned[1]-owned[0]+1;
    int Y = owned[3]-owned[2]+1;
    // Set pixel factor
    pf.first  = static_cast<float>(x)/X;
    pf.second = static_cast<float>(y)/Y;
    // Return
    return pf;
  }

  //! @brief Set a pixel relative to the owned coordinates (origin is {owned[0], owned[2]})
  //!
  //! EasyBMP has the top left corner as its origin, it would be easier to use the bottom left, so
  //! this function also takes care of that.
  void Palette::setPixel(int x, int y, RGBApixel color) {
    int X = x+owned[0], Y = owned[3]-y-1;
    // Check the bounds
    if (X<0 || bounds[1]<=X || Y<0 || bounds[3]<=Y) throw PaletteOutOfBounds();
    // Set the pixel
    image->SetPixel(X, Y, color);
  }

}