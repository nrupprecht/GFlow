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
    wd = default_wall_width;
    // Set mass and moment of inertia to be infinite
    im = iI = 0;
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
    wd = default_wall_width;
    // Set mass and moment of inertia to be infinite
    im = iI = 0;
  }

  void Wall::clearForceTorque() {
    fx = fy = tq = 0;
  }
}
