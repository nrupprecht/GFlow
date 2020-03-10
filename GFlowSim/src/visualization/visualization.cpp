#include "visualization.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "../compute/load_data.hpp"
#include "../compute/store_data.hpp"

namespace GFlowSimulation {

  Visualization::Visualization() : do_wrap(true), v_scale_average(2) {
    setColorBankSize(10);
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
    v_scale_average = Vec(dimensions);
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
    // I just picked a random large enough number. There isn't really any overhead to having lots of colors.
    if (color_option==6) setColorBankSize(64); 

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

  void Visualization::setColorSelectionMethod(const unsigned int method) {
    color_selection_method = method;
  }

  void Visualization::setSelectionName(const string name) {
    selection_name = name;
  }

  void Visualization::setColorFunctionSelection(const unsigned int f_number) {
    color_function_selection = f_number;
  }

  void Visualization::setRadiusMultiple(float r) {
    radius_multiple = r;
  }

  void Visualization::setColorOption(int opt) {
    color_selection_method = 0;
    color_option = opt; 
    findPlaces();
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
    // Find data positions
    findPlaces();
  }

  void Visualization::createVideo2d(string dirName, const vector<vector<float> >& data) {
    cout << "Starting image write.\n";
    // Timing
    auto start_time = current_time();
    
    // Create frames
    for (int i=0; i<data.size(); ++i) {
      string fileName = dirName + "/frame" + toStr(i) + ".bmp";
      frame = i;
      createImage(fileName, data[i]);
    }
    frame = -1;
    // Timing
    auto end_time = current_time();
    // Print timing.
    cout << "Image creation time: " << time_span(end_time, start_time) << endl;
  }

  void Visualization::createVideo3d(string dirName, const vector<vector<float> >& data) {
    cout << "Starting image write.\n";
    // Timing
    auto start_time = current_time();

    if (choice_3d) {
      // Set up the camera
      standard_camera_setup();

      // Create all the images
      for (int i=0; i<data.size(); ++i) {
        string fileName = dirName + "/frame" + toStr(i) + ".bmp";
        frame = i;
        createImage3d(fileName, data[i]);
        // projectImage(fileName, data[i]); // For now, we just draw a projected image
      }
      // Clean up the ray tracer's kd tree structure.
      tracer.empty();
    }
    else {
      // Create all the images
      for (int i=0; i<data.size(); ++i) {
        string fileName = dirName + "/frame" + toStr(i) + ".bmp";
        frame = i;
        projectImage(fileName, data[i]); // For now, we just draw a projected image
      }
    }
    // Timing
    auto end_time = current_time();
    // Print timing.
    cout << "Image creation time: " << time_span(end_time, start_time) << endl;
  }

  void Visualization::createImage(string fileName, const vector<float>& data) {
    // Get some data from the bounds
    float buffer = 0.f;
    float wx = bounds.wd(0) + 2*buffer, wy = bounds.wd(1) + 2*buffer, left = bounds.min[0] - buffer, bott = bounds.min[1] - buffer;
    // Figure out the needed resolution
    int res_x = resolution, res_y = resolution;
    if (wx>wy) res_y = wy/wx*resolution;
    else if (wy>wx) res_x = wx/wy*resolution;
    // Create the main palette
    Palette palette(res_x, res_y);
    palette.coverPalette(RGB_Black);
    Palette sub1 = palette.getSubPaletteCentered(0.01); // Get a centered subset
    sub1.coverPalette(RGB_White);
    Palette image = sub1.getMaxCenteredSubPalette(wx, wy); // Get a subpallete on which to write the image
    image.setSpaceBounds(wx, wy);
    image.coverPalette(background);
    // Do checks of positions
    if (!do_checks(data)) return;
    determineScaleFactors(data);
    // Draw all particles
    for (int i=0; i<data.size(); i+=dataWidth) {
      // Get values
      const float *pos, *v_data; 
      float sigma, s_data; 
      int i_data;
      get_values(&data[i], pos, sigma, v_data, s_data, i_data);     
      // Find the center of the particle
      float xf = (pos[0] - left)/wx;
      float yf = (pos[1] - bott)/wy;
      float rf = sigma/wx;
      // --- Determine the color
      RGBApixel color = RGB_Green;
      determine_color(color, i, v_data, s_data, i_data);
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
    determineScaleFactors(data);
    // Draw all particles
    for (int i=0; i<data.size(); i+=dataWidth) {
      // Get values
      const float *pos, *v_data; 
      float sigma, s_data; 
      int i_data;
      get_values(&data[i], pos, sigma, v_data, s_data, i_data);
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
      determine_color(color, i, v_data, s_data, i_data);
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
    // Save image
    palette.writeToFile(fileName);
  }

  void Visualization::createImage3d(string fileName, const vector<float>& data) {
    // Set up the camera if it hasn't been set up before.
    if (!camera_set) standard_camera_setup();
    // Add all objects to the ray tracer.
    tracer.reserve(data.size()/dataWidth);
    tracer.setResolution(resolution);
    RGBApixel color;
    // Do checks of positions
    if (!do_checks(data)) return;
    determineScaleFactors(data);
    // Add all particles to the ray tracer
    const float *pos, *v_data; 
    float sigma, s_data; 
    int i_data;
    for (int i=0; i<data.size(); i+=dataWidth) {
      get_values(&data[i], pos, sigma, v_data, s_data, i_data);
      // Determine the color
      determine_color(color, i, v_data, s_data, i_data);   
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

  inline void Visualization::findVectorScaleFactors(const vector<float>& data_frame) {
    if (v_select<0 || dataWidth<=v_select) {
      cout << "Error: invalid vector selection.\n";
      return;
    }
    v_scale_max = 0.f;
    v_scale_average = Vec(dimensions);
    for (int i=v_select; i<data_frame.size(); i+=dataWidth) {
      float magnitude = magnitudeVec(&data_frame[i], dimensions);
      if (magnitude>v_scale_max) v_scale_max = magnitude;
      plusEqVec(v_scale_average.data, &data_frame[i], dimensions);
    }
    v_scale_average *= (dataWidth/data_frame.size());
  }

  inline void Visualization::findScalarScaleFactors(const vector<float>& data_frame) {
    if (s_select<0 || dataWidth<=s_select) {
      cout << "Error: invalid scalar selection.\n";
      return;
    }
    s_scale_min = data_frame[s_select], s_scale_max = data_frame[s_select];
    s_scale_average = 0;
    for (int i=s_select; i<data_frame.size(); i+=dataWidth) {
      real scalar_value = data_frame[i];
      if (scalar_value<s_scale_min) s_scale_min = scalar_value;
      else if (scalar_value>s_scale_max) s_scale_max = scalar_value;
      s_scale_average += scalar_value;
    }
    s_scale_average *= (dataWidth/data_frame.size());
  }

  inline void Visualization::determineScaleFactors(const vector<float>& data_frame) {
    switch (color_selection_method) {
      case 0: {
        break;
      }
      case 1: {
        findVectorScaleFactors(data_frame);
        break;
      }
      case 2: {
        findScalarScaleFactors(data_frame);
        break;
      }
      case 3: {
        break;
      }
    }
  }

  inline RGBApixel Visualization::getColor(int c) const {
    return colorBank[ c%colorBank.size() ];
  }

  inline void Visualization::resetPlaces() {
    pos_place = -1;
    sg_place = -1;
    v_select = -1;
    s_select = -1;
    i_select = -1;
  }

  inline void Visualization::findPlaces() {
    auto find = [] (string key, vector<string>& vec) -> int {
      for (int i=0; i<vec.size(); ++i)
        if (key==vec[i]) return i;
      return -1;
    };

    if (color_selection_method==0) {
      // Use a predefined method to color the particles.
      // 0 - Type.
      // 1 - Random coloring.
      // 2 - Velocity magnitude.
      // 3 - Orientation.
      // 4 - Distance traveled.
      // 5 - x-stripe.
      // 6 - processor #.
      // 7 - angular velocity.
      // 8 - Force magnitude.
      switch (color_option) {
        case 0: {
          chooseSelectionType(3, "Type");
          break;
        }
        case 1: {
          color_selection_method = 0;
          break;
        }
        case 2: {
          if (!chooseSelectionType(2, "V-M")) 
            chooseSelectionType(1, "V");
          break;
        }
        case 3: {
          chooseSelectionType(1, "V");
          color_function_selection = 1;
          break;
        }
        case 4: {
          chooseSelectionType(2, "Distance");
          break;
        }
        case 5: {
          chooseSelectionType(2, "StripeX");
          break;
        }
        case 6: { // Color by proc #.
          chooseSelectionType(3, "Proc");
          break;
        }
        case 7: { // Color by angular velocity.
          chooseSelectionType(2, "Om");
          break;
        }
        case 8: { // Color by force.
          chooseSelectionType(2, "F-M");
          break;
        }
      }
    }

    // Vector data
    pos_place = find("X", vector_data)*dimensions;
    if (color_selection_method==1) v_select = find(selection_name, vector_data)*dimensions;
    else v_select = -1;

    // Scalar data
    int shift = vector_data.size()*dimensions;
    sg_place = find("Sg", scalar_data) + shift;
    if (color_selection_method==2) s_select = find(selection_name, scalar_data) + shift;
    else s_select = -1;

    // Integer data
    shift += scalar_data.size();
    if (color_selection_method==3) i_select = find(selection_name, integer_data) + shift;
    else i_select = -1;
  }

  inline bool Visualization::chooseSelectionType(const unsigned try_selection_method, const string& try_name) {
    switch (try_selection_method) {
      case 0:
        return true;
      case 1: {
        auto it = std::find(vector_data.begin(), vector_data.end(), try_name);
        if (it!=vector_data.end()) {
          color_selection_method = try_selection_method;
          selection_name = try_name;
          return true;
        }
        else return false;
      }
      case 2: {
        auto it = std::find(scalar_data.begin(), scalar_data.end(), try_name);
        if (it!=scalar_data.end()) {
          color_selection_method = try_selection_method;
          selection_name = try_name;
          return true;
        }
        else return false;
      }
      case 3: {
        auto it = std::find(integer_data.begin(), integer_data.end(), try_name);
        if (it!=integer_data.end()) {
          color_selection_method = try_selection_method;
          selection_name = try_name;
          return true;
        }
        else return false;
      }
      default:
        return false;
    }
  }

  inline bool Visualization::do_checks(const vector<float>& data) {
    // Check what data we have, make sure we can make the requested image
    if (dataWidth==0) {
      cout << "Data width is zero. This means there is no data, or something is very wrong. Exiting.\n";
      return false;
    }
    if (pos_place<0) {
      cout << "No position data detected. This is essential. Exiting.\n";
      return false;
    }
    if (sg_place<0) {
      cout << "No radius data detected. This is essential. Exiting.\n";
      return false;
    }
    return true;
  }

  inline void Visualization::get_values(
    const float * pdata,
    const float * &pos, 
    float& sigma, 
    const float * &v_data, 
    float &s_data, 
    int &i_data
  ) {
    // Essential data.
    pos = pos_place<0 ? nullptr : &pdata[pos_place];
    sigma = sg_place<0 ? 0 : pdata[sg_place]; // Get sigma
    // Data for coloring particles.
    v_data = v_select>=0 ? &pdata[v_select] : nullptr;
    s_data = s_select>=0 ? pdata[s_select] : 0.f;
    i_data = i_select>=0 ? static_cast<int>(pdata[i_select]) : 0;
  }

  inline void Visualization::determine_color(RGBApixel& color, int i, const float *v_data, float &s_data, int &i_data) {
    switch (color_selection_method) {
      case 0: {
        // Color randomly, by order in the data. This will cause particles to switch colors when exchanging particles.
        color = getColor(i);
        break;
      }
      case 1: {
        switch (color_function_selection) {
          default:
          case 0: {
            // Color by magnitude.
            float value = magnitudeVec(v_data, dimensions) / v_scale_max;
            color = RGBApixel(floor(255*value), 0, 200*(1-value));
            break;
          }
          case 1: {
            // Color by direction - only for 2+ dimensions, and only considers x,y directions.
            float theta = atan2(v_data[1] - v_scale_average[1], v_data[0] - v_scale_average[0]);
            color = colorAngle(theta);
            break;
          }
        }
        break;
      }
      case 2: {
        float value = 0.f;
        // Color by magnitude.
        switch (color_function_selection) {
          default:
          case 0: {
            value = (s_data - s_scale_min) / (s_scale_max - s_scale_min);
            break;
          }
          case 1: {
            value = log( 1. + 100.*s_data / s_scale_max) / log(101.);
            break;
          }
        }
        color = RGBApixel(floor(255*value), 0, ceil(200*(1-value)));
        break;
      }
      case 3: {
        // Color by integer.
        color = getColor(i_data);
        break;
      }
    }
  }

  inline void Visualization::standard_camera_setup() {
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
    // We have set up the camera
    camera_set = true;
  }

}
