#include "visualization.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "../compute/load_data.hpp"

namespace GFlowSimulation {

  Visualization::Visualization() : pos_place(0), vel_place(-1), sg_place(-1), type_place(-1), 
    distance_place(-1), resolution(1.5*1024), do_wrap(true), background(RGB_Black), maxVsqr(-1), maxDistance(-1), color_option(0)
  {
    // Default size - 10
    createColorBank(10);
  };

  bool Visualization::load_and_create(string loadName, string saveName) {
    // Create a loader to load the data
    LoadData loader;
    if (!loader.load(loadName)) return false;
    // Extract data from the loader
    int dataWidth = loader.getDataWidth(), ntypes = loader.getNTypes();
    dimensions = loader.getDimensions();
    Bounds bounds = loader.getBounds();
    vector_data = loader.get_vector_data();
    scalar_data = loader.get_scalar_data();
    integer_data = loader.get_integer_data();
    // Find data positions
    findPlaces();
    // Get the bulk of the data
    const vector<vector<double> > &data = loader.getData();

    // If we are coloring by type, we need a large enough color bank
    if (color_option==0 && colorBank.size()<ntypes) setColorBankSize(ntypes);

    // Create a video from the data
    if (dimensions==2) createVideo2d(saveName, data, dataWidth, bounds, dimensions);
    if (dimensions>2)  createVideo3d(saveName, data, dataWidth, bounds, dimensions);

    // Reset places so they will be gathered again next time
    resetPlaces();

    return true;
  }

  void Visualization::setColorBankSize(int cbs) {
    if (cbs!=colorBank.size() && cbs>0) {
      colorBank = vector<RGBApixel>(cbs);
      for (int i=0; i<cbs; ++i)
      colorBank[i] = randomColor();
    }
  }

  void Visualization::setRadiusMultiple(double r) {
    radius_multiple = r;
  }

  void Visualization::setColorOption(int opt) {
    color_option = opt;
  }

  void Visualization::setResolution(int r) {
    resolution = r;
  }

  void Visualization::createVideo2d(string dirName, const vector<vector<double> >& data, int dataWidth, Bounds& bounds, int dimensions) const 
  {
    // Find the maximum velocity (if needed)
    if (color_option==2)
      findMaxVSqr(data, dataWidth);
    else if (color_option==4)
      findMaxDistance(data, dataWidth);
    // Create frames
    for (int i=0; i<data.size(); ++i) {
      string fileName = dirName + "/frame" + toStr(i) + ".bmp";
      createImage(fileName, data[i], dataWidth, bounds, dimensions);
    }
  }

  void Visualization::createVideo3d(string dirName, const vector<vector<double> >& data, int dataWidth, Bounds& bounds, int dimensions) const {
    for (int i=0; i<data.size(); ++i) {
      string fileName = dirName + "/frame" + toStr(i) + ".bmp";
      createImage3d(fileName, data[i], dataWidth, bounds, dimensions);
    }
  }

  void Visualization::createImage(string fileName, const vector<double>& data, int dataWidth, Bounds& bounds, int dimensions) const {
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

    // Make sure we have the appropriate normalization factors.
    // This will generally only happen when we are making a single image.
    if (color_option==2 && maxVsqr<0) {
      auto pack_data = vector<vector<double> >(1, data);
      findMaxVSqr(pack_data, dataWidth);
    }
    else if (color_option==4 && maxDistance<0) {
      auto pack_data = vector<vector<double> >(1, data);
      findMaxDistance(pack_data, dataWidth);
    }
    if (stripex_place<0) {
      color_option = 0;
    }
    
    // A vector for holding data
    double *pdata = new double[dataWidth];
    // Print all particles
    for (int i=0; i<data.size(); i+=dataWidth) {
      // Extract data for the i-th particle
      for (int j=0; j<dataWidth; ++j)
        pdata[j] = data[i+j];
      // Get individual entries
      double *pos  = pdata;
      double *vel  = vel_place>-1    ? &pdata[vel_place] : nullptr; // Point to start of velocity data
      double sigma = sg_place>-1     ?  pdata[sg_place]  : 0; // Get sigma
      int type       = type_place > -1 ? static_cast<int>(pdata[type_place]) : 0; // Get type
      double distance = distance_place > -1 ? pdata[distance_place] : 0; // Get distance traveled
      double stripex = stripex_place<0 ? 0 : pdata[stripex_place];
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
              float V = magnitudeVec(vel, dimensions)/maxVsqr;
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
          case 4: { // Color by distance
            float D = maxDistance==0 ? 1. : log(1.f + distance)/log(1.f + maxDistance);
            color = RGBApixel(floor(255*D), floor(255*(1-D)), 0);
            break;
          }
          case 5: { // Color by xstripe
            RealType width = bounds.wd(1);
            const int nstripes = 20;
            int s = (stripex - bott)/wy * nstripes;
            if (s%2) color = RGB_Green;
            else color = RGB_White;
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
      image.drawCircle(pos[0] - left, pos[1] - bott, sigma, colorF, do_wrap);
    }

    // Clean up pdata
    delete [] pdata;

    // Reset
    maxVsqr = -1;
    maxDistance = -1;

    // Save image
    palette.writeToFile(fileName);
  }

  void Visualization::createImage3d(string fileName, const vector<double>& data, int dataWidth, Bounds& bounds, int dimensions) const {
    // @todo Implement this.
    throw false; 
  }

  inline void Visualization::findMaxVSqr(const vector<vector<double> >& dataVector, int dataWidth) const {
    // Reset
    maxVsqr = 0;
    // Look for max vsqr
    for (int iter=0; iter<dataVector.size(); ++iter) {
      if (dataVector[iter].empty()) continue;
      const double *data = &dataVector[iter][0];
      int number = dataVector[iter].size()/dataWidth;
      for (int i=0; i<number; ++i) {
        float V = sqr(data[i*dataWidth + vel_place]);
        if (V>maxVsqr) maxVsqr = V;
      }
    }
    // Take the sqrt
    maxVsqr = sqrt(maxVsqr);
  }

  inline void Visualization::findMaxDistance(const vector<vector<double> >& dataVector, int dataWidth) const {
    // Reset
    maxDistance = 0;
    // Look for max distance - start after first iteration.
    for (int iter=0; iter<dataVector.size(); ++iter) {
      if (dataVector[iter].empty()) continue;
      const double *data = &dataVector[iter][0];
      int number = dataVector[iter].size()/dataWidth;
      for (int i=0; i<number; ++i) {
        float D = data[i*dataWidth + distance_place];
        if (D>maxDistance) maxDistance = D;
      }
    }
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

  inline void Visualization::resetPlaces() {
    pos_place = -1;
    vel_place = -1;
    sg_place = -1;
    type_place = -1;
    distance_place = -1;
    stripex_place = -1;
  }

  inline void Visualization::findPlaces() {
    auto find = [] (string key, vector<string>& vec) -> int {
      for (int i=0; i<vec.size(); ++i)
        if (key==vec[i]) return i;
      return -1;
    };

    // Vector data
    pos_place = find("X", vector_data)*dimensions;
    vel_place = find("V", vector_data)*dimensions;

    // Scalar data
    int shift = vector_data.size()*dimensions;
    sg_place = find("Sg", scalar_data) + shift;
    stripex_place = find("StripeX", scalar_data) + shift;    

    // Integer data
    shift += scalar_data.size();
    type_place = find("Type", integer_data) + shift;
  }

}