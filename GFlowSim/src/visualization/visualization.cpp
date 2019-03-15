#include "visualization.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "../compute/load_data.hpp"
#include "../compute/store_data.hpp"

namespace GFlowSimulation {

  Visualization::Visualization() {
    createColorBank(10); // Default size - 10
    resolution = 1536;
  };

  bool Visualization::load_and_create(string loadName, string saveName) {
    // Create a loader to load the data
    auto start_load = current_time();
    cout << "Starting load...\n";
    LoadData loader;
    if (!loader.load(loadName)) return false;
    auto  end_load = current_time();
    // Print load time
    cout << "Load time: " << time_span(end_load, start_load) << endl;
    // Extract data from the loader
    dataWidth = loader.getDataWidth();
    ntypes = loader.getNTypes();
    dimensions = loader.getDimensions();
    bounds = loader.getBounds();
    // Get entries
    vector_data = loader.get_vector_data();
    scalar_data = loader.get_scalar_data();
    integer_data = loader.get_integer_data();
    // Find data positions
    findPlaces();
    // Get the bulk of the data
    const vector<vector<float> > &data = loader.getData();

    // If we are coloring by type, we need a large enough color bank
    if (color_option==0 && colorBank.size()<ntypes) setColorBankSize(ntypes);

    // Create a video from the data
    if (dimensions==2) createVideo2d(saveName, data);
    if (dimensions>2)  createVideo3d(saveName, data);

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

  void Visualization::setRadiusMultiple(float r) {
    radius_multiple = r;
  }

  void Visualization::setColorOption(int opt) {
    color_option = opt;
  }

  void Visualization::setResolution(int r) {
    resolution = r;
  }

  void Visualization::setMetaParameters(const StoreData& storeData) {
    dataWidth = storeData.getDataWidth();
    ntypes = storeData.getNTypes();
    dimensions = storeData.getDimensions();
    bounds = storeData.getBounds();
    // Get entries
    vector_data = storeData.get_vector_data();
    scalar_data = storeData.get_scalar_data();
    integer_data = storeData.get_integer_data();
    // Find data positoins
    findPlaces();
  }

  void Visualization::createVideo2d(string dirName, const vector<vector<float> >& data) {
    cout << "Starting image write.\n";
    // Timing
    auto start_time = current_time();
    // Find the maximum velocity (if needed)
    if (color_option==2)
      findMaxVSqr(data);
    // Find the maximum distance (if needed)
    else if (color_option==4)
      findMaxDistance(data);
    // Create frames
    for (int i=0; i<data.size(); ++i) {
      string fileName = dirName + "/frame" + toStr(i) + ".bmp";
      createImage(fileName, data[i]);
    }
    // Timing
    auto end_time = current_time();
    // Print timing.
    cout << "Image creation time: " << time_span(end_time, start_time) << endl;
  }

  void Visualization::createVideo3d(string dirName, const vector<vector<float> >& data) {
    cout << "Starting image write.\n";
    // Timing
    auto start_time = current_time();
    // Find the maximum velocity (if needed)
    if (color_option==2)
      findMaxVSqr(data);
    // Find the maximum distance (if needed)
    else if (color_option==4)
      findMaxDistance(data);

    // Set up the ray tracer - it should be empty
    tracer.setBounds(bounds);
    // Set the tracer's camera
    float bounds_center[3];
    float scale = 0.5*max_width(bounds);
    float camera[3] = {1.3f, 2.4f, 2.f};
    // Scale the camera placement vector
    scalarMultVec(scale, camera, 3);
    // Find the center of the bounds
    addVec(bounds.min, bounds.max, bounds_center, 3); 
    scalarMultVec(0.5f, bounds_center, 3); 
    // Shift the camera placement vector
    plusEqVec(camera, bounds_center, 3);

    // Have the orientation point towards the center of the bounds
    float orientation[3];
    subtractVec(bounds_center, camera, orientation, 3);
    normalizeVec(orientation, 3);

    // Set camera and orientation vectors
    tracer.setCameraPlacement(camera);
    tracer.setCameraOrientation(orientation);

    // Create all the images
    for (int i=0; i<data.size(); ++i) {
      string fileName = dirName + "/frame" + toStr(i) + ".bmp";
      createImage3d(fileName, data[i]);
      // projectImage(fileName, data[i]); // For now, we just draw a projected image
    }
    // Clean up the ray tracer's kd tree structure.
    tracer.empty();
    // Timing
    auto end_time = current_time();
    // Print timing.
    cout << "Image creation time: " << time_span(end_time, start_time) << endl;
  }

  void Visualization::createImage(string fileName, const vector<float>& data) {
    // Get some data from the bounds
    float wx = bounds.wd(0), wy = bounds.wd(1), left = bounds.min[0], bott = bounds.min[1];
    // Figure out the needed resolution
    int res_x = resolution, res_y = resolution;
    if (wx>wy) res_y = wy/wx*resolution;
    else if (wy>wx) res_x = wx/wy*resolution;
    // Create the main palette
    Palette palette(resolution, resolution);
    palette.coverPalette(RGB_Black);
    Palette sub1 = palette.getSubPaletteCentered(0.01); // Get a centered subset
    sub1.coverPalette(RGB_White);
    Palette image = sub1.getMaxCenteredSubPalette(wx, wy); // Get a subpallete on which to write the image
    image.setSpaceBounds(wx, wy);
    image.coverPalette(background);
    // Do checks of positions
    if (!do_checks(data)) return;
    // Draw all particles
    for (int i=0; i<data.size(); i+=dataWidth) {
      // Get valudes
      const float *pos, *vel; float sigma, distance, stripex; int type;
      get_values(&data[i], pos, vel, sigma, type, distance, stripex);
      // If type<0, continue
      if (type<0) continue;
      // Find the center of the particle
      float xf = (pos[0] - left)/wx;
      float yf = (pos[1] - bott)/wy;
      float rf = sigma/wx;
      // --- Determine the color
      RGBApixel color = RGB_Green;
      determine_color(color, i, pos, vel, type, distance, stripex);
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
      image.drawCircle(pos[0] - left, pos[1] - bott, radius_multiple*sigma, colorF, do_wrap);
    }
    // Reset
    maxVsqr = -1;
    maxDistance = -1;
    // Save image
    palette.writeToFile(fileName);
  }

  void Visualization::projectImage(string fileName, const vector<float>& data) {
    // Get some data from the bounds
    float wx = bounds.wd(0), wy = bounds.wd(1), left = bounds.min[0], bott = bounds.min[1];
    // Figure out the needed resolution
    int res_x = resolution, res_y = resolution;
    if (wx>wy) res_y = wy/wx*resolution;
    else if (wy>wx) res_x = wx/wy*resolution;
    // Create the main palette
    Palette palette(resolution, resolution);
    palette.coverPalette(RGB_Black);
    Palette sub1 = palette.getSubPaletteCentered(0.01); // Get a centered subset
    sub1.coverPalette(RGB_White);
    Palette image = sub1.getMaxCenteredSubPalette(wx, wy); // Get a subpallete on which to write the image
    image.setSpaceBounds(wx, wy);
    image.coverPalette(background);
    // Do checks of positions
    if (!do_checks(data)) return;
    // Draw all particles
    for (int i=0; i<data.size(); i+=dataWidth) {
      // Get valudes
      const float *pos, *vel; float sigma, distance, stripex; int type;
      get_values(&data[i], pos, vel, sigma, type, distance, stripex);
      // If type<0, continue
      if (type<0) continue;
      // Find the center of the particle
      float xf = (pos[0] - left)/wx;
      float yf = (pos[1] - bott)/wy;
      float d_tangent_sqr = 0;
      for (int d=2; d<dimensions; ++d) d_tangent_sqr += sqr(pos[d]);
      // If particle does not intersect with the plane of projection
      if (d_tangent_sqr>=sqr(sigma)) continue;
      float rf = sqrt(1. - d_tangent_sqr/sqr(sigma));

      // --- Determine the color
      RGBApixel color = RGB_Green;
      determine_color(color, i, pos, vel, type, distance, stripex);
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
      image.drawCircle(pos[0] - left, pos[1] - bott, rf*radius_multiple*sigma, colorF, do_wrap);
    }
    // Reset
    maxVsqr = -1;
    maxDistance = -1;
    // Save image
    palette.writeToFile(fileName);
  }

  void Visualization::createImage3d(string fileName, const vector<float>& data) {
    // Add all objects to the ray tracer.
    tracer.reserve(data.size()/dataWidth);
    tracer.setResolution(resolution);
    RGBApixel color;
    // Do checks of positions
    if (!do_checks(data)) return;
    // Add all particles to the ray tracer
    for (int i=0; i<data.size(); i+=dataWidth) {
      // Get values
      const float *pos, *vel; float sigma, distance, stripex; int type;
      get_values(&data[i], pos, vel, sigma, type, distance, stripex);
      // Determine the color
      determine_color(color, i, pos, vel, type, distance, stripex);
      // Add the sphere to the ray tracer.
      tracer.addSphere(pos, sigma, color);
    }
    // Tells the ray tracer we are done adding objects. The tracer constructs a KD tree, and gets ready to render.
    tracer.initialize();
    // Render via ray tracing and produce an image.
    tracer.render();
    tracer.empty();
    // Write the image to the file
    tracer.saveImage(fileName);
  }

  inline void Visualization::findMaxVSqr(const vector<vector<float> >& dataVector) {
    // Check that we have velocity data
    if (vel_place<0) return;
    // Reset
    maxVsqr = 0;
    // Look for max vsqr
    for (int iter=0; iter<dataVector.size(); ++iter) {
      if (dataVector[iter].empty()) continue;
      const float *data = &dataVector[iter][0];
      int number = dataVector[iter].size()/dataWidth;
      for (int i=0; i<number; ++i) {
        float V = sqr(data[i*dataWidth + vel_place]);
        if (V>maxVsqr) maxVsqr = V;
      }
    }
    // Take the sqrt
    maxVsqr = sqrt(maxVsqr);
  }

  inline void Visualization::findMaxDistance(const vector<vector<float> >& dataVector) {
    // Check that we have distance data
    if (distance_place<0) return;
    // Reset
    maxDistance = 0;
    // Look for max distance - start after first iteration.
    for (int iter=0; iter<dataVector.size(); ++iter) {
      if (dataVector[iter].empty()) continue;
      const float *data = &dataVector[iter][0];
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

  inline bool Visualization::do_checks(const vector<float>& data) {
    // Check what data we have, make sure we can make the requested image
    if (dataWidth==0 || pos_place<0) {
      cout << "No position data detected. This is essential. Exiting.\n";
      return false;
    }
    else if ((color_option==2 || color_option==3) && vel_place<0) {
      cout << "No velocity data detected. Switching to color option 0.\n";
      color_option = 0;
    }
    else if (color_option==4 && distance_place<0) {
      cout << "No distance data detected. Switching to color option 0.\n";
      color_option = 0;
    }
    else if (color_option==5 && stripex_place<0) {
      cout << "No stripe-x data detected. Switching to color option 0.\n";
      color_option = 0;
    }
    // Make sure we have the appropriate normalization factors.
    // This will generally only happen when we are making a single image.
    if (color_option==0 && ntypes>colorBank.size()) {
      createColorBank(ntypes);
    }
    if (color_option==2 && maxVsqr<0) {
      auto pack_data = vector<vector<float> >(1, data);
      findMaxVSqr(pack_data);
    }
    else if (color_option==4 && maxDistance<0) {
      auto pack_data = vector<vector<float> >(1, data);
      findMaxDistance(pack_data);
    }

    // Return success
    return true;
  }

  inline void Visualization::get_values(const float * pdata, const float * &pos, const float * &vel, float& sigma, int& type, float& distance, float& stripex) {
    // Get individual entries
    pos = pos_place<0 ? nullptr : &pdata[pos_place];
    vel = vel_place<0 ? nullptr : &pdata[vel_place]; // Point to start of velocity data
    sigma = sg_place<0 ? 0 : pdata[sg_place]; // Get sigma
    type = type_place<0 ? 0 : static_cast<int>(pdata[type_place]); // Get type
    distance = distance_place<0 ? 0 : pdata[distance_place]; // Get distance traveled
    stripex  = stripex_place<0  ? 0 : pdata[stripex_place];
  }

  inline void Visualization::determine_color(RGBApixel& color, int i, const float *pos, const float *vel, int type, float distance, float stripex) {
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
          const int nstripes = 5;
          int s = (stripex - bounds.min[1])/bounds.wd(1) * nstripes;
          int c = s%4;
          switch (c) {
            case 0:
              color = RGB_Green;
              break;
            case 1:
              color = RGB_Blue;
              break;
            case 2:
              color = RGB_White;
              break;
            case 3:
              color = RGB_Red;
              break;
          }
        }
      }
    }
  }

}