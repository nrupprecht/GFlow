#ifndef __VISUALIZATION_HPP__GFLOW__
#define __VISUALIZATION_HPP__GFLOW__

#include "palette.hpp"
#include <map>

namespace GFlowSimulation {

  class DataLayout {
  public:
    //! @brief Constructor. Pass in data sizes and data names.
    DataLayout(vector<int>& ds, vector<string>& dn) { 
      if (ds.size()!=dn.size()) throw false;
      data_width = 0;
      for (int i=0; i<ds.size(); ++i) {
        name_place_mapping.insert(std::pair<string, int>(dn.at(i), ds.at(i)));
        data_width += ds[i];
      }
    }

    RealType *access(RealType *data, string &name) {
      auto f = name_place_mapping.find(name);
      if (f==name_place_mapping.end()) return nullptr;
      else return &data[f->second];
    } 
  private:
    //! @brief The total length of a data entry.
    int data_width;
    //! @brief A mapping between data names, and the place where that data entry starts.
    std::map<string, int> name_place_mapping;
  };

  /**
  *  @brief Creates visualization from
  *
  */
  class Visualization {
  public:
    //! @brief Constructor.
    Visualization();

    //! @brief Destructor.
    ~Visualization();

    bool load_and_create(string, string);

    void setColorBankSize(int);

    void setRadiusMultiple(RealType);

    void setColorOption(int);

    void setResolution(int);

    //! @brief Create a directory filled with BMP renderings of the system.
    //!
    //! This can be used to create a movie.
    void createVideo2d(string, const vector<vector<RealType> >&, int, BoundsPack&, int) const;

    void createVideo3d(string, const vector<vector<RealType> >&, int, BoundsPack&, int) const;

  private:
    //! @brief Creates a single frame.
    inline void createImage(string, const vector<RealType>&, int, BoundsPack&, int) const;

    inline void createImage3d(string, const vector<RealType>&, int, BoundsPack&, int) const;

    inline void findMaxVSqr(const vector<vector<RealType> >&, int) const;

    inline void createColorBank(int);

    inline RGBApixel getColor(int) const;

    template<typename T> inline T getNextNumber(std::ifstream& fin) const {
      char c;
      fin.get(c);
      stringstream stream;
      while (!fin.eof() && c!=',' && c!='\n') {
        if (isdigit(c) || c=='.' || c=='-') stream << c;
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
    int pos_place;

    //! @brief Where the velocity data starts.
    int vel_place;

    //! @brief Where in the data for a particle is the radius.
    int sg_place;

    //! @brief Where in the data for a particle is its type.
    int type_place;

    //! @brief The dimensions of the image (it will be the same in x and y)
    int resolution;

    //! @brief Whether to wrap at the boundaries or not
    bool do_wrap;

    mutable RealType maxVsqr;

    RGBApixel background;

    unsigned int color_option;

    vector<RGBApixel> colorBank;

    double radius_multiple = 1.;
  };

}
#endif // __VISUALIZATION_HPP__GFLOW__