#include "visualization.hpp"
#include "vectormath.hpp"

#include <opencv2/opencv.hpp>

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
    return Palette(mx, Mx, my, My, image, refs, bounds);
  }

  void Palette::drawCircleByFactors(float xf, float yf, float rf, ColorFunction colorF, bool wrap) {
    // Width of owned region
    int wx = owned[1] - owned[0];
    int wy = owned[3] - owned[2];
    // Find the center
    int px = xf*wx + owned[0];
    int py = yf*wy + owned[2];
    // Find dp's
    int dpx = ceil(rf*(owned[1] - owned[0]));
    int dpy = ceil(rf*(owned[3] - owned[2]));

    // Radius factor squared
    float rsqr = sqr(rf);

    for (int dy=-dpy; dy<dpy; ++dy) {
      for (int dx=-dpx; dx<dpx; ++dx) {
        if (sqr(static_cast<float>(dx)/wx) + sqr(static_cast<float>(dy)/wy) <= rsqr) {
          // The potential pixel
          int w_px = px + dx, w_py = py + dy;
          // Wrap pixel
          if (wrap) {
            w_px = w_px >= owned[1] ? w_px-wx : w_px;
            w_px = w_px <  owned[0] ? w_px+wx : w_px;
            w_py = w_py >= owned[3] ? w_py-wy : w_py;
            w_py = w_py <  owned[2] ? w_py+wy : w_py;
            // Set pixel
            image->SetPixel(w_px, w_py, 
              colorF(static_cast<float>(dx)/wx, static_cast<float>(dy)/wy)
            );
          }
          else if (owned[0]<=w_px && w_px<owned[1] && owned[2]<=w_py && w_py<owned[3]) {
            image->SetPixel(w_px, w_py, 
              colorF(static_cast<float>(dx)/wx, static_cast<float>(dy)/wy)
            );
          }
        }
      }
    }
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

  //-------------

  Visualization::Visualization() : sg_place(DIMENSIONS), resolution(1024), do_wrap(true), background(Black) {
    colorBank = new RGBApixel[10];
    for (int i=0; i<10; ++i)
      colorBank[i] = randomColor();
  };

  Visualization::~Visualization() {
    if (colorBank) delete [] colorBank;
  }

  void Visualization::createBMPs(string dirName, const vector<RealType*>& data, const vector<int>& number, 
    int dataWidth, Bounds& bounds, int dimensions) const 
  {
    for (int i=0; i<data.size(); ++i) {
      string fileName = dirName + "/frame" + toStr(i) + ".bmp";
      createImage(fileName, data[i], number[i], dataWidth, bounds, dimensions);
    }
  }

  inline void Visualization::createImage(string fileName, RealType *data, int number, int dataWidth, Bounds& bounds, int dimensions) const {
    // Get some data from the bounds
    double wx = bounds.wd(0)/resolution;
    double wy = bounds.wd(1)/resolution;
    double left = bounds.min[0];
    double bott = bounds.min[1];

    BMP image;
    image.SetSize(resolution, resolution);

    // Set background
    for (int j=0; j<resolution; ++j)
      for (int i=0; i<resolution; ++i)
        image.SetPixel(i, j, background);

    // Print all particles
    for (int i=0; i<number; ++i) {
      // Get the data for this particle
      RealType *pdata = &data[i*dataWidth];
      // Find the center of the particle
      int px = (pdata[0] - left)/wx;
      int py = (pdata[1] - bott)/wy;
      int dp = pdata[sg_place]/wx;
      // Draw a circle
      circle(image, px, py, dp, colorBank[i%10], do_wrap);
    }

    // Save image
    image.WriteToFile(fileName.c_str());
  }

  inline void Visualization::circle(BMP& image, int x, int y, int r, RGBApixel color, bool wrap) const {
    for (int dy=-r; dy<r; ++dy) {
      for (int dx=-r; dx<r; ++dx) {
        if (sqr(dx) + sqr(dy) <= sqr(r)) {
          // The potential pixel
          int w_px = x + dx, w_py = y + dy;
          // Wrap pixel
          if (wrap) {
            w_px = w_px>=resolution ? w_px-resolution : w_px;
            w_px = w_px<0 ? w_px+resolution : w_px;
            w_py = w_py>=resolution ? w_py-resolution : w_py;
            w_py = w_py<0 ? w_py+resolution : w_py;
            // Set pixel
            image.SetPixel(w_px, w_py, color);
          }
          else if (0<=w_px && w_px<resolution && 0<=w_py && w_py<resolution) {
            image.SetPixel(w_px, w_py, color);
          }
        }
      }
    }
  }


}