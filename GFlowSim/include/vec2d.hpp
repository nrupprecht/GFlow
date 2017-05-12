#ifndef __VEC_2D_HPP__
#define __VEC_2D_HPP__

#include "Utility.h"

namespace GFlow {
  
  /*
   * @class vec2
   *
   */
  struct vec2 {
    // Default constructor
    vec2() : x(0), y(0) {};

    // Initialized constructor
    vec2(RealType x, RealType y) : x(x), y(y) {};

    // The actual vector data
    RealType x, y;
  };

  const vec2 Zero(0,0);

}
#endif // __VEC_2D_HPP__
