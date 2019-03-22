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

    bool load_and_create(string, string);

    void setColorBankSize(int);

    void setRadiusMultiple(float);

    void setColorOption(int);

    void setResolution(int);

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

    inline void findMaxVSqr(const vector<vector<float> >&);

    inline void findMaxDistance(const vector<vector<float> >&);

    inline void createColorBank(int);

    inline RGBApixel getColor(int) const;

    inline void resetPlaces();

    inline void findPlaces();

    inline bool do_checks(const vector<float>&);

    inline void get_values(const float*, const float* &, const float* &, float&, int&, float&, float&);

    //! \brief Determine what color a particle should be.
    inline void determine_color(RGBApixel&, int, const float*, const float*, int, float, float);

    //! \brief Set up the camera to be in the standard position.
    inline void standard_camera_setup();

    //! \brief The dimensions of the image (it will be the same in x and y)
    int resolution = 1536;

    //! \brief Where the position data starts.
    int pos_place = -1;

    //! \brief Where the velocity data starts.
    int vel_place = -1;

    //! \brief Where in the data for a particle is the radius.
    int sg_place = -1;

    //! \brief Where in the data for a particle is its type.
    int type_place = -1;

    //! \brief Where in the data for a particle is the distance traveled.
    int distance_place = -1;

    //! \brief Where in the data the stripe data is.
    int stripex_place = -1;

    //! \brief The width of the data.
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

    //! \brief Whether the (3D) camera has been set up by anyone or anything.
    bool camera_set = false;

    float maxVsqr = -1.;
    float maxDistance = -1.;
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
  };

}
#endif // __VISUALIZATION_HPP__GFLOW__