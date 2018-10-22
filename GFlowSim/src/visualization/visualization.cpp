#include "visualization.hpp"
// Other files
#include "../utility/vectormath.hpp"

namespace GFlowSimulation {

  Visualization::Visualization() : pos_place(0), vel_place(DIMENSIONS), sg_place(2*DIMENSIONS), type_place(2*DIMENSIONS+1), 
    resolution(1.5*1024), do_wrap(true), background(RGB_Black), maxVsqr(0), color_option(0)
  {
    // Default size - 10
    createColorBank(10);
  };

  bool Visualization::load_and_create(string loadName, string saveName) {
    // Open file
    std::ifstream fin(loadName);
    if (fin.fail()) return false;

    // Get the data width, dimensions, number of times sampled, and the number of types
    int dataWidth, dimensions, samples, ntypes;
    dataWidth = getNextNumber<int>(fin);
    dimensions = getNextNumber<int>(fin);
    samples = getNextNumber<int>(fin);
    ntypes = getNextNumber<int>(fin);

    // Get the dimensions - mins first, then maxes
    BoundsPack bounds(dimensions);
    for (int i=0; i<dimensions; ++i) 
      bounds.min[i] = getNextNumber<RealType>(fin);
    for (int i=0; i<dimensions; ++i)
      bounds.max[i] = getNextNumber<RealType>(fin);

    // Vector to hold data
    vector<vector<RealType> > data;

    // The first number in each line is the number of particles to expect
    for (int iter=0; iter<samples; ++iter) {
      // Get the length of data to expect
      int data_length = getNextNumber<int>( fin);
      // Get this iter's data
      vector<RealType> pdata;
      for (int i=0; i<data_length; ++i) {
        RealType datum = getNextNumber<RealType>(fin);
        pdata.push_back(datum);
      }
      // Store this iter's data vector
      data.push_back(pdata);
    }

    // Close the file stream
    fin.close();

    // If we are coloring by type, we need a large enough color bank
    if (color_option==0 && colorBank.size()<ntypes) setColorBankSize(ntypes);

    // Create a video from the data
    if (dimensions==2) createVideo2d(saveName, data, dataWidth, bounds, dimensions);
    if (dimensions>2)  createVideo3d(saveName, data, dataWidth, bounds, dimensions);

    return true;
  }

  void Visualization::setColorBankSize(int cbs) {
    if (cbs!=colorBank.size() && cbs>0) {
      colorBank = vector<RGBApixel>(cbs);
      for (int i=0; i<cbs; ++i)
      colorBank[i] = randomColor();
    }
  }

  void Visualization::setRadiusMultiple(RealType r) {
    radius_multiple = r;
  }

  void Visualization::setColorOption(int opt) {
    color_option = opt;
  }

  void Visualization::setResolution(int r) {
    resolution = r;
  }

  void Visualization::createVideo2d(string dirName, const vector<vector<RealType> >& data, int dataWidth, BoundsPack& bounds, int dimensions) const 
  {
    // Find the maximum velocity (if needed)
    if (color_option==2)
      findMaxVSqr(data, dataWidth);
    // Create frames
    for (int i=0; i<data.size(); ++i) {
      string fileName = dirName + "/frame" + toStr(i) + ".bmp";
      createImage(fileName, data[i], dataWidth, bounds, dimensions);
    }
  }

  void Visualization::createVideo3d(string dirName, const vector<vector<RealType> >& data, int dataWidth, BoundsPack& bounds, int dimensions) const {
    for (int i=0; i<data.size(); ++i) {
      string fileName = dirName + "/frame" + toStr(i) + ".bmp";
      createImage3d(fileName, data[i], dataWidth, bounds, dimensions);
    }
  }

  inline void Visualization::createImage(string fileName, const vector<RealType>& data, int dataWidth, BoundsPack& bounds, int dimensions) const {
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
    Palette palette(resolution, resolution);
    palette.coverPalette(RGB_Black);

    // Get a centered subset
    Palette sub1 = palette.getSubPaletteCentered(0.01);
    sub1.coverPalette(RGB_White);

    // Get a subpallete on which to write the image
    Palette image = sub1.getMaxCenteredSubPalette(wx, wy);
    image.setSpaceBounds(wx, wy);

    // Set background
    image.coverPalette(background);
    
    // A vector for holding data
    RealType *pdata = new RealType[dataWidth];
    // Print all particles
    for (int i=0; i<data.size(); i+=dataWidth) {
      // Extract data for the i-th particle
      for (int j=0; j<dataWidth; ++j)
        pdata[j] = data[i+j];
      // Get individual entries
      RealType *pos  = pdata;
      RealType *vel  = vel_place>-1    ? &pdata[vel_place] : nullptr; // Point to start of velocity data
      RealType sigma = sg_place>-1     ?  pdata[sg_place]  : 0; // Get sigma
      int type       = type_place > -1 ? static_cast<int>(pdata[type_place]) : 0; // Get type
      // If type<0, continue
      if (type<0) continue;
      // Find the center of the particle
      float xf = (pos[0] - left)/wx;
      float yf = (pos[1] - bott)/wy;
      float rf = sigma/wx;

      // --- Determine the color
      RGBApixel color = RGB_Green;
      if (!colorBank.empty()) {
        switch (color_option) {
          default:
          case 0: { // Color by type
            color = getColor(type);
            break;
          }
          case 1: { // Color randomly
            color = getColor(i);
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
      //image.drawCircleByFactors(xf, yf, rf*radius_multiple, colorF, do_wrap); 

      image.drawCircle(pos[0] - left, pos[1] - bott, sigma, colorF, do_wrap);
    }

    // Clean up pdata
    delete [] pdata;

    // Save image
    palette.writeToFile(fileName);
  }

  inline void Visualization::createImage3d(string fileName, const vector<RealType>& data, int dataWidth, BoundsPack& bounds, int dimensions) const {
    
  }

  inline void Visualization::findMaxVSqr(const vector<vector<RealType> >& dataVector, int dataWidth) const {
    // Reset
    maxVsqr = 0;
    // Look for max vsqr
    for (int iter=0; iter<dataVector.size(); ++iter) {
      if (dataVector[iter].empty()) continue;
      const RealType *data = &dataVector[iter][0];
      int number = dataVector[iter].size()/dataWidth;
      for (int i=0; i<number; ++i) {
        float V = sqr(data[i*dataWidth + vel_place]);
        if (V>maxVsqr) maxVsqr = V;
      }
    }
    // Take the sqrt
    maxVsqr = sqrt(maxVsqr);
  }

  inline void Visualization::createColorBank(int size) {
    if (size>0) {
      colorBank.clear();
      for (int i=0; i<size; ++i)
        colorBank.push_back(randomColor());
    }
  }

  inline RGBApixel Visualization::getColor(int c) const {
    return colorBank[ c%colorBank.size() ];
  }

}