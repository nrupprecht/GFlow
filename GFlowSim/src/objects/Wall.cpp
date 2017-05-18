#include "Wall.hpp"

namespace GFlow {
  
  Wall::Wall() : left(Zero) {
    length = 0;
    normal = Zero;

    repulsion = default_wall_repulsion;
    dissipation = default_wall_dissipation;
    coeff = default_wall_coeff;
  }

  Wall::Wall(vec2 left, vec2 right) {
    this->left = left;
    length = sqrt(sqr(right-left));
    normal = (1./length)*(right-left);

    repulsion = default_wall_repulsion;
    dissipation = default_wall_dissipation;
    coeff = default_wall_coeff;
  }
  
  Wall::Wall(RealType lx, RealType ly, RealType rx, RealType ry) {
    left = vec2(lx, ly); vec2 right(rx, ry);
    length = sqrt(sqr(right-left));
    normal = (1./length)*(right-left);

    repulsion = default_wall_repulsion;
    dissipation = default_wall_dissipation;
    coeff = default_wall_coeff;
  }

}
