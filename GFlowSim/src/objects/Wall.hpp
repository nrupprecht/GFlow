/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 12, 2017
 *
 */

#ifndef __WALL_HPP__
#define __WALL_HPP__

// Includes
#include "../../include/vec2d.hpp"
#include "../../include/CSVUtility.hpp"
#include "../../include/DefaultConstants.hpp"

namespace GFlow {

  /*
   * @class Wall
   * Defines a stationary wall
   *
   */
  struct Wall {
    Wall();
    Wall(vec2, vec2);
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
  
  // Adjust displacement to do the proper minimum displacement between a point and a wall
  inline void wallDisplacement(vec2& displacement, const RealType sigma, const Wall &w) {
    // We are given displacement = p.position - w.left;
    double l_par = displacement*w.normal;
    vec2 d_par = l_par*w.normal;
    vec2 d_perp = displacement - d_par;
    // Check whether the particle is between the start and end of the wall
    if (l_par>=0) { // Located forward of the origin
      if (w.length>l_par) displacement = d_perp;  // The particle is above the wall (in the perp. direction)
      else displacement -= w.length*w.normal; // Displacement from the nearest end (the far end) of the wall
    }
  }
  
}
#endif // __WALL_HPP__
