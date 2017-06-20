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
    // Make sure left is actually on the left
    if (right.x<left.x) swap(left, right);
    // If left and right x values coincide, left is the bottom corner
    if (left.x==right.x && right.y<left.y) swap(left, right);
    // Set positions
    this->left = left;
    length = sqrt(sqr(right-left));
    normal = (1./length)*(right-left);
    // Set values
    repulsion = default_wall_repulsion;
    dissipation = default_wall_dissipation;
    coeff = default_wall_coeff;
  }
  
  Wall::Wall(RealType lx, RealType ly, RealType rx, RealType ry) {
    // Create vectors
    left = vec2(lx, ly); vec2 right(rx, ry);
    // Make sure left is actually on the left
    if (right.x<left.x) swap(left, right);
    // If left and right x values coincide, left is the bottom corner
    if (left.x==right.x && right.y<left.y) swap(left, right);
    // Set positions
    length = sqrt(sqr(right-left));
    normal = (1./length)*(right-left);
    // Set values
    repulsion = default_wall_repulsion;
    dissipation = default_wall_dissipation;
    coeff = default_wall_coeff;
  }

  void Wall::setRight(vec2 right) {
    *this = Wall(left, right);
  }

}
