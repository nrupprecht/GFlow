#ifndef __WALL_HPP__
#define __WALL_HPP__

// Includes
#include "../../include/vec2d.hpp"
#include "../../include/DefaultConstants.hpp"

namespace GFlow {

  /*
   * @class Wall
   *
   */
  struct Wall {
    Wall();
    Wall(RealType, RealType, RealType, RealType);

    vec2 left;
    RealType length;
    vec2 normal;
    RealType repulsion, dissipation, coeff;
  };

}
#endif // __WALL_HPP__
