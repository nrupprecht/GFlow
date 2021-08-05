#ifndef __FIELD_PROPERTIES_HPP__GFLOW__
#define __FIELD_PROPERTIES_HPP__GFLOW__

#include "binning2d.hpp"

namespace GFlowSimulation {

  struct vec2d {
    vec2d() : x(0), y(0) {};
    vec2d(float v) : x(v), y(v) {};
    vec2d(float a, float b) : x(a), y(b) {};

    float x, y;
  };

  inline string toStr(vec2d v) {
    return toStr(v.x) + "," + toStr(v.y);
  }

  template<typename T> class Field2d {
  public:
    Field2d() : dimx(0), dimy(0) {};

    Field2d(int x, int y) : dimx(x), dimy(y) {
      data = vector<T>(x*y, T(0));
    }

    void setSize(int x, int y) {
      if (x<1 || y<1) throw BadDimension(); // @todo This is not actually the right error to use.
      dimx = x; dimy = y;
      data = vector<T>(x*y, T(0));
    }

    void setBounds(const Bounds& bounds) {
      min[0] = bounds.min[0];
      min[1] = bounds.min[1];
      max[0] = bounds.max[0];
      max[1] = bounds.max[1];
      // Set up dims
      dx = (max[0]-min[0])/dimx;
      dy = (max[1]-min[1])/dimy;
    }

    //! @brief Get at entry.
    T& at(int x, int y) {
      return data.at(y*dimx+x);
    }

    //! @brief Set all entries to zero.
    void clear() {
      for (auto &v : data) v = T(0);
    }

    string toMMA() {
      string str = "{";
      for (int y=0; y<dimy; ++y) {
        for (int x=0; x<dimx; ++x) {
          str += "{" + toStr(x*dx+min[0]) + "," + toStr(y*dy+min[1]) + "," + toStr(data.at(y*dimx+x)) + "}";
          if (!(y==dimy-1 && x==dimx-1)) str += ",";
        }
      }
      str += "}";
      return str;
    }

    string toCSV() {
      string str;
      for (int y=0; y<dimy; ++y)
        for (int x=0; x<dimx; ++x)
          str += toStr(x*dx+min[0]) + "," + toStr(y*dy+min[1]) + "," + toStr(data.at(y*dimx+x)) + '\n';
      return str;
    }

  private:
    //! @brief Number of data points in each dimension.
    int dimx, dimy;

    // @brief Dimensions of a bin.
    float dx, dy;
    //! @brief Width of the bounds.
    float min[2];
    float max[2];

    //! @brief The data.
    vector<T> data;


  };


  class FieldProperties {
  public:

    void create_field(const vector<float>&, const Bounds&, int, int);

    void load_and_create(string);

    string toStrMMA_phiField();
    string toStrMMA_aveVField();
    string toStrMMA_varVField();

    string toCSV_phiField();
    string toCSV_aveVField();
    string toCSV_varVField();
    string toCSV_aveVVectorField();

  private:

    inline float vel(int, int, const vector<float>&);

    inline void setPlaces(const int);

    //! @brief Where the position data starts.
    int pos_place = -1;

    //! @brief Where the velocity data starts.
    int vel_place = -1;

    //! @brief Where in the data for a particle is the radius.
    int sg_place = -1;

    //! @brief Where in the data for a particle is its type.
    int type_place = -1;

    //! @brief Where in the data for a particle is the distance traveled.
    int distance_place = -1;

    //! @brief The data width.
    int data_width = -1;

    // --- Fields

    Field2d<float> phi_field;
    Field2d<float> ave_v_field;
    Field2d<float> var_v_field;

    Field2d<vec2d> ave_v_vector_field;
  };

}
#endif // __FIELD_PROPERTIES_HPP__GFLOW__