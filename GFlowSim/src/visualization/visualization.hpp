#ifndef __VISUALIZATION_HPP__GFLOW__
#define __VISUALIZATION_HPP__GFLOW__

#include "palette.hpp"
#include <map>

namespace GFlowSimulation {

  /**
  *  @brief Creates visualizations from data.
  *
  */
  class Visualization {
  public:
    //! @brief Constructor.
    Visualization();

    bool load_and_create(string, string);

    void setColorBankSize(int);

    void setRadiusMultiple(double);

    void setColorOption(int);

    void setResolution(int);

    //! @brief Create a directory filled with BMP renderings of the system.
    //!
    //! This can be used to create a movie.
    void createVideo2d(string, const vector<vector<double> >&, int, Bounds&, int) const;

    void createVideo3d(string, const vector<vector<double> >&, int, Bounds&, int) const;

    //! @brief Creates a single frame.
    void createImage(string, const vector<double>&, int, Bounds&, int) const;

    void createImage3d(string, const vector<double>&, int, Bounds&, int) const;

  private:

    inline void findMaxVSqr(const vector<vector<double> >&, int) const;

    inline void findMaxDistance(const vector<vector<double> >&, int) const;

    inline void createColorBank(int);

    inline RGBApixel getColor(int) const;

    inline void setPlaces(const int) const;

    inline void resetPlaces();

    template<typename T> inline T getNextNumber(std::ifstream& fin) const {
      char c;
      fin.get(c);
      stringstream stream;
      while (!fin.eof() && c!=',' && c!='\n') {
        if (isdigit(c) || c=='.' || c=='-' || c=='e' || c=='E') stream << c;
        else;
        // Get the next character
        fin.get(c);
      }
      // Convert to a type T and return
      T val;
      stream >> val;
      return val;
    }

    //! @brief Where the position data starts.
    mutable int pos_place;

    //! @brief Where the velocity data starts.
    mutable int vel_place;

    //! @brief Where in the data for a particle is the radius.
    mutable int sg_place;

    //! @brief Where in the data for a particle is its type.
    mutable int type_place;

    //! @brief Where in the data for a particle is the distance traveled.
    mutable int distance_place;

    //! @brief The dimensions of the image (it will be the same in x and y)
    int resolution;

    //! @brief Whether to wrap at the boundaries or not
    bool do_wrap;

    mutable float maxVsqr;

    mutable float maxDistance;

    RGBApixel background;

    unsigned int color_option;

    vector<RGBApixel> colorBank;

    double radius_multiple = 1.;
  };

}
#endif // __VISUALIZATION_HPP__GFLOW__