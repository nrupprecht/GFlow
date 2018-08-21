#include "visualization.hpp"
#include "vectormath.hpp"

namespace GFlowSimulation {

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