#include "Wall.hpp"
#include "../control/SimData.hpp"

namespace GFlow {
  
  Wall::Wall() : px(0), py(0), vx(0), vy(0), fx(0), fy(0), th(0), om(0), sg(0), im(0), iI(0) {
    rp = default_wall_repulsion;
    ds = default_wall_dissipation;
    cf = default_wall_coeff;
    wd = default_wall_width;
  }

  Wall::Wall(vec2 left, vec2 right) {
    RealType lx = left.x, ly = left.y, rx = right.x, ry = right.y;
    // Center of mass position
    px = 0.5*(lx + rx);
    py = 0.5*(ly + ry);
    // Angle
    th = atan2(ry-ly, rx-lx);
    // Sigma (characteristic length - half the total length of the wall)
    sg = sqrt(sqr(right - left));
    // Set velocities, forces, torques to zero
    vx = vy = fx = fy =om = tq = 0;
    // Set default values
    rp = default_wall_repulsion;
    ds = default_wall_dissipation;
    cf = default_wall_coeff;
    // Set mass and moment of inertia to be infinite
    im = iI = 0;
    /*
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
    */
  }
  
  Wall::Wall(RealType lx, RealType ly, RealType rx, RealType ry) {
    // Center of mass position
    px = 0.5*(lx + rx);
    py = 0.5*(ly + ry);
    // Angle
    th = atan2(ry-ly, rx-lx);
    // Sigma (characteristic length - half the total length of the wall)
    sg = sqrt(sqr(rx - lx) + sqr(ry-ly));
    // Set velocities, forces, torques to zero
    vx = vy = fx = fy = om = tq = 0;
    // Set default values
    rp = default_wall_repulsion;
    ds = default_wall_dissipation;
    cf = default_wall_coeff;
    // Set mass and moment of inertia to be infinite
    im = iI = 0;

    /*
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
    */    
  }

  /*
  void Wall::setRight(vec2 right) {
    *this = Wall(left, right);
  }
  */

  void Wall::clearForceTorque() {
    fx = fy = tq = 0;
  }

  /*
  // The displacement between a point and the closest point on a wall
  inline vec2 wallDisplacement(const vec2& pos, const Wall& w, SimData *simData) {
    // Calculate displacement and normal vectors
    vec2 displacement = simData->getDisplacement(pos.x, pos.y, w.px, w.py);
    vec2 normal(cos(w.th), sin(w.th)), perpendicular(-normal.y, normal.x);
    // Get components of the displacement in the "normal/perpendicular" framc
    RealType dn = displacement*normal, dp = displacement*perpendicular;
    // Calculate the displacement between the center of the particle and the closest point on the wall
    vec2 minimalDisplacement = dp*perpendicular;
    if (fabs(dn)>w.sg) minimalDisplacement += sign(dn)*(dn-w.sg)*normal;
    // Return the displacement
    return minimalDisplacement;
  }
  */

}
