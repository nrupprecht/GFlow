#ifndef __WALL_HPP__
#define __WALL_HPP__

// Includes
#include "../../include/vec2d.hpp"
#include "../../include/CSVUtility.hpp"
#include "../../include/DefaultConstants.hpp"

namespace GFlow {

  /*
   * @class Wall
   *
   */
  struct Wall {
    Wall();
    Wall(RealType, RealType, RealType, RealType);

    // Get right and left
    vec2 getLeft()  const { return left; }
    vec2 getRight() const { return left+length*normal; }

    // Wall data
    vec2 left;
    RealType length;
    vec2 normal;
    RealType repulsion, dissipation, coeff;
  };

  inline string toCSV(const Wall& w) {
    stringstream stream;
    string str;
    stream << toCSV(w.getLeft()) << "," << toCSV(w.getRight());
    stream >> str;
    return str;
  }

}
#endif // __WALL_HPP__
