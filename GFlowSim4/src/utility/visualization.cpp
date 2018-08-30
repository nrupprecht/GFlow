#include "visualization.hpp"
// Other files
#include "vectormath.hpp"

namespace GFlowSimulation {

  Visualization::Visualization() : pos_place(0), vel_place(DIMENSIONS), sg_place(2*DIMENSIONS), type_place(2*DIMENSIONS+1), 
    resolution(1024*1.5), do_wrap(true), background(RGB_Black), maxVsqr(0), color_option(0)
  {
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
    // Find the maximum velocity (if needed)
    if (color_option==2)
      findMaxVSqr(data, number, dataWidth);
    // Create frames
    for (int i=0; i<data.size(); ++i) {
      string fileName = dirName + "/frame" + toStr(i) + ".bmp";
      createImage(fileName, data[i], number[i], dataWidth, bounds, dimensions);
    }
  }

  void Visualization::createVideo2d(string fileName, const vector<RealType*>& data, int dataWidth) {
    vector<int> data_size = { 2, 2, 1 };
    vector<string> data_name = { "x", "v", "sg" };
  }

  inline void Visualization::createImage(string fileName, RealType *data, int number, int dataWidth, Bounds& bounds, int dimensions) const {
    // Get some data from the bounds
    float wx = bounds.wd(0);
    float wy = bounds.wd(1);
    float left = bounds.min[0];
    float bott = bounds.min[1];

    // Figure out the needed resolution
    int res_x = resolution, res_y = resolution;
    if (wx>wy) res_y = wy/wx*resolution;
    else if (wy>wx) res_x = wx/wy*resolution;

    // Create the main palette
    Palette palette(res_x, res_y);
    palette.coverPalette(RGB_Black);
    // Get a subpallete on which to write the image
    Palette image = palette.getSubPalette(10, res_x-10, 10, res_y-10);

    // Set background
    image.coverPalette(background);

    // Print all particles
    for (int i=0; i<number; ++i) {
      // Get the data for this particle
      RealType *pdata = &data[i*dataWidth];
      // Get individual entries
      RealType *vel  = vel_place>-1 ? &pdata[vel_place] : nullptr;
      RealType sigma = sg_place>-1 ? pdata[sg_place] : 0;
      int type       = type_place > -1 ? static_cast<int>(pdata[type_place]) : 0;
      // Find the center of the particle
      float xf = (pdata[0] - left)/wx;
      float yf = (pdata[1] - bott)/wy;
      float rf = sigma/wx;

      // --- Determine the color
      RGBApixel color = RGB_Green;
      if (colorBank) {
        switch (color_option) {
          default:
          case 0: { // Color by type
            color = colorBank[type];
            break;
          }
          case 1: { // Color randomly
            int id = i%10;
            color = colorBank[id];
            break;
          }
          case 2: { // Color by velocity
            if (vel) {
              float V = magnitudeVec(vel)/maxVsqr;
              color = RGBApixel(floor(255*V), 0, 200*(1-V));
            }
            break;
          }
          case 3: { // Color by orientation
            if (vel) {
              float theta = atan2(vel[1], vel[0]);
              color = colorAngle(theta);
            }
            break;
          }
        }
      }

      // --- Create the color function.
      // xf, yf \in [-1, 1] (roughly - could be a bit larger)
      std::function<RGBApixel(float, float, bool&)> colorF = 
        [&] (float xf, float yf, bool &doColor)  {
          float s = 255.f*(1.f - 0.5f*(sqr(xf) + sqr(yf)));
          if (s<0) {
            doColor = false;
            return RGB_White;
          }
          else {
            doColor = true;
            return color*makePixel(floor(s)); 
          }
        };

      // Draw the particle
      image.drawCircleByFactors(xf, yf, rf, colorF, do_wrap); 
    }

    // Save image
    palette.writeToFile(fileName);
  }

  inline void Visualization::findMaxVSqr(const vector<RealType*>& dataVector, const vector<int>& numbers, int dataWidth) const {
    // Reset
    maxVsqr = 0;
    // Look for max vsqr
    for (int iter=0; iter<dataVector.size(); ++iter) {
      RealType *data = dataVector.at(iter);
      int number = numbers.at(iter);
      for (int i=0; i<number; ++i) {
        float V = sqr(data[i*dataWidth + vel_place]);
        if (V>maxVsqr) maxVsqr = V;
      }
    }
    // Take the sqrt
    maxVsqr = sqrt(maxVsqr);
  }

}