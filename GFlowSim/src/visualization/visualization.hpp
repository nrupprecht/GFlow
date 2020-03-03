#ifndef __VISUALIZATION_HPP__GFLOW__
#define __VISUALIZATION_HPP__GFLOW__

#include "palette.hpp"
#include "raytrace.hpp"
#include <map>

namespace GFlowSimulation {

  /**
  *  \brief Creates visualizations from data.
  *
  *  This class can create both images and videos. Data is passed in as a vector<float> (for images), or 
  *  a vector<vector<float>> (for movies). Meta-data such as the physical bounds, data width, number of
  *  particle types, etc. must be set from the outside if using the createVideo or create image functions,
  *  but will be taken from the load object if using the load_and_create function.
  */
  class Visualization {
  public:
    //! \brief Constructor.
    Visualization();

    //! Read a csv file and create images from all the data in the file.
    bool load_and_create(string, string);

    //! \brief Set the size of the color bank.
    void setColorBankSize(int);

    //! \brief Choose the color selection method.
    void setColorSelectionMethod(const unsigned int);

    //! \brief Choose the name of the data entry that should be processed to get the color data.
    void setSelectionName(const string);

    //! \brief Select what function should be used to turn data into a color.
    void setColorFunctionSelection(const unsigned int);

    //! \brief Set a radius multiplier factor.
    //!
    //! This option was used before the variable cutoff system was used, and is therefore mostly obsolete.
    void setRadiusMultiple(float);

    //! \brief Set the coloring option.
    void setColorOption(int);

    //! \brief Set the resolution.
    void setResolution(int);

    //! \brief Set the "meta-parameters" of this visualization class based on the store data object.
    // This would include the number of dimensions, number of types, the data entry names, etc.
    void setMetaParameters(const class StoreData&);

    //! \brief Create a directory filled with BMP renderings of the system.
    //!
    //! This can be used to create a movie.
    void createVideo2d(string, const vector<vector<float> >&);

    void createVideo3d(string, const vector<vector<float> >&);

    //! \brief Creates a single frame of 2D data.
    void createImage(string, const vector<float>&);

    //! \brief Create a single frame from 2+ dimensional data projected into the first two dimensions.
    void projectImage(string, const vector<float>&);

    //! \brief Creates a single frame of 3D data.
    void createImage3d(string, const vector<float>&);

  private:

    //! \brief Find the scale factors for the requested vector data.
    inline void findVectorScaleFactors(const vector<float>&);
    //! \brief Find the scale factors for the requested scalar data.
    inline void findScalarScaleFactors(const vector<float>&);

    //! \brief Find the scale factors for whatever type of data was requested.
    inline void determineScaleFactors(const vector<float>&);

    //! \brief Access a color from the color bank.
    inline RGBApixel getColor(int) const;

    //! \brief Reset all place pointers.
    inline void resetPlaces();

    //! \brief Find all the requisite place pointers.
    inline void findPlaces();

    //! \brief Do whatever checks are needed to make sure everything is in order to render a frame.
    inline bool do_checks(const vector<float>&);

    //! \brief Get the needed values from a particle data entry. What is needed will include the position and radius,
    //! and whatever data is needed for rendering.
    inline void get_values(const float*, const float*&, float&, const float*&, float&, int&);

    //! \brief Determine what color a particle should be.
    inline void determine_color(RGBApixel&, int, const float*, float&, int&);

    //! \brief Set up the camera to be in the standard position.
    inline void standard_camera_setup();

    //! \brief The current frame being animated.
    int frame = -1;

    //! \brief The dimensions of the image (it will be the same in x and y)
    int resolution = 1536;

    //! \brief Where the position data starts.
    int pos_place = -1;

    //! \brief Where in the data for a particle is the radius.
    int sg_place = -1;

    //! \brief The width of the data for a single particle.
    int dataWidth = 0;

    //! \brief The number of particle types.
    int ntypes = 0;

    //! \brief The dimensionality of vector data.
    int dimensions = 0;

    //! \brief The physical bounds of the data.
    Bounds bounds = Bounds(2);

    //! \brief A ray tracer, for 3d image creation.
    RayTrace tracer;

    //! \brief Whether to wrap at the boundaries or not
    bool do_wrap = true;

    //! \brief Whether to create a 3D image, or just a projection when rendering in 3D
    //!
    //! True - 3D image, False - projection.
    bool choice_3d = true;

    //! \brief Whether the (3D) camera has been set up by anyone or anything.
    bool camera_set = false;

    //! \brief How much larger/smaller the drawn radius should be compared to the radius specified in the data.
    //! The default is, of course, 1.
    float radius_multiple = 1.;

    //! \brief How the particles should be colored.
    unsigned int color_option = 0;

    //! \brief A bank of colors.
    vector<RGBApixel> colorBank;

    //! \brief The background color.
    RGBApixel background = RGB_Black;

    // The names of the data
    vector<string> vector_data;
    vector<string> scalar_data;
    vector<string> integer_data;

    //! \brief Choice of how to color the particles. 
    //!
    //! 0 - Other
    //! 1 - Use vector data.
    //! 2 - Use scalar data.
    //! 3 - Use integer data.
    unsigned int color_selection_method = 0;

    //! \brief The name of the entry that should be used to color the particles.
    string selection_name;

    //! \brief A flag that determines which function should be applied to the chosen data to produce a color.
    unsigned int color_function_selection = 0;

    // Only one of these selectors will actually be used. Which one is used is determined by the color_selection_method flag.

    //! \brief Pointer to the place in a particle data entry where the required vector data starts.
    int v_select = -1;
    //! \brief Pointer to the place in a particle data entry where the required scalar data is.
    int s_select = -1;
    //! \brief Pointer to the place in a particle data entry where the required integer data is.
    int i_select = -1;

    float v_scale_max = 0.f;
    Vec v_scale_average;
    float s_scale_min = 0.f, s_scale_max = 1.f, s_scale_average;
    int i_scale_min = 0, i_scale_max = 1;
  };

}
#endif // __VISUALIZATION_HPP__GFLOW__