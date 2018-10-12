#include "Characteristic.hpp"
#include "../control/SimData.hpp"

namespace GFlow {

  void ConstantVelocity::modify(SimData* simData, int id, RealType dt, bool modP) {
    // Checks
    if (!active) return;
    // Load either particle or wall data
    RealType &px = modP ? simData->getPxPtr() [id] : simData->getWalls() [id].px;
    RealType &py = modP ? simData->getPyPtr() [id] : simData->getWalls() [id].py;
    RealType &vx = modP ? simData->getVxPtr() [id] : simData->getWalls() [id].vx;
    RealType &vy = modP ? simData->getVyPtr() [id] : simData->getWalls() [id].vy;
    // Keep velocity and omega constant
    if (useV) {
      // Update time
      time += dt;
      // Set position
      RealType x = pos.x + velocity.x*time, y =  pos.y + velocity.y*time;
      // Make sure to wrap the position
      simData->wrap(x,y);
      // Manually set the position
      px = x;
      py = y;
      // Set velocities so dissipative forces work correctly
      vx = velocity.x;
      vy = velocity.y;
    }
    if (useOm) {
      if (modP) simData->getOmPtr() [id] = omega;
      else simData->getWalls() [id].om = omega;
    }
    // Check if we should deactivate -- temporary code (?)
    RealType close = (simData->getPy(id) - simData->getSimBounds().bottom)/simData->getSg(id);
    // Stop if we are supposed to
    if (stop && close<5) {
      active = false;
      simData->getVx(id) = 0;
      simData->getVy(id) = 0;
    }
  }
  
  void ConstantVelocity::reset() {
    time = 0;
    active = true;
    stop = false;
  }

  void Insertion::modify(SimData* simData, int id, RealType dt, bool modP) {
    if (!active) return;
    // Keep velocity and omega constant
    if (useV) {
      simData->getVx(id) = velocity.x;
      simData->getVy(id) = velocity.y;
    }
    if (useOm) simData->getOm(id) = omega;
    // Check if we should reverse our direction
    RealType close = (simData->getPy(id) - simData->getSimBounds().bottom)/simData->getSg(id);
    if (forward && close<5) {
      forward = false;
      if (useV)  velocity = -velocity;
      if (useOm) omega = -omega;
    }
    // Keep track of how far we are from our starting point
    distance += (forward ? speed*dt : -speed*dt);
    // Terminate once we're back (actually, when we have gone past being back)
    if (distance<-3*simData->getSg(id)) simData->setTerminate(true);
  }

  void Circulate::modify(SimData* simData, int id, RealType dt, bool modP) {
    // Checks
    if (modP==false) throw false;
    if (!active) return;
    if (first) {
      RealType px = simData->getPx(id), py = simData->getPy(id);
      // Set center
      px -= radius*cos(theta); py -= radius*sin(theta);
      center = vec2(px, py);
      // Done
      first = false;
    }
    // Update angle
    simData->getPx(id) = center.x + radius*cos(theta);
    simData->getPy(id) = center.y + radius*sin(theta);
    theta += dt*omega;
  }

  Fixed::Fixed(int id, SimData* sd, bool modP) {
    // If modifying a particle
    if (modP) {
      pos.x = sd->getPxPtr() [id];
      pos.y = sd->getPyPtr() [id];
    }
    // If modifying a wall
    else {
      pos.x = sd->getWalls() [id].px;
      pos.y = sd->getWalls() [id].py;
    } 
  }

  void Fixed::modify(SimData* simData, int id, RealType dt, bool modP) {
    if (!active) return;
    // If modifying a particle
    if (modP) {
      simData->getFxPtr() [id] = 0;
      simData->getFyPtr() [id] = 0;
      simData->getVxPtr() [id] = 0;
      simData->getVyPtr() [id] = 0;
      simData->getPxPtr() [id] = pos.x;
      simData->getPyPtr() [id] = pos.y;
    }
    // If modifying a wall
    else {
      auto &wall = simData->getWalls() [id];
      wall.fx = 0;
      wall.fy = 0;
      wall.vx = 0;
      wall.vy = 0;
      wall.px = pos.x;
      wall.py = pos.y;
    }
  }

  ApplyForce::ApplyForce(vec2 df) : F0(Zero), F(Zero), dF(df) {};

  ApplyForce::ApplyForce(vec2 f0, vec2 df) : F0(f0), F(f0), dF(df) {};

  void ApplyForce::modify(SimData* simData, int id, RealType dt, bool modP) {
    if (!active) return;
    F += dt*dF;
    // Modify either particle or wall
    RealType &fx = modP ? simData->getFxPtr() [id] : simData->getWalls() [id].fx;
    RealType &fy = modP ? simData->getFyPtr() [id] : simData->getWalls() [id].fy;
    fx = F.x;
    fy = F.y;
  }

  void ApplyForce::reset() {
    F = F0;
  }

}