#ifndef __WALL_HPP__
#define __WALL_HPP__

namespace GFlow {

  /*
   * @class Wall
   *
   */
  struct Wall {
    Wall() {}; // STUB

    vec2 left;
    RealType length;
    vec2 normal;
    RealType repulsion, dissipation, coeff;
  };

}
#endif // __WALL_HPP__
