#include "Characteristic.hpp"
#include "../control/SimData.hpp"

namespace GFlow {

  void ConstantVelocity::modify(SimData* simData, int id, RealType dt) {
    // Keep velocity and omega constant
    if (active) {
      if (useV) {
	// Update time
	time += dt;
	// Set position
	RealType x = pos.x + velocity.x*time, y =  pos.y + velocity.y*time;
	simData->wrap(x,y);
	simData->getPxPtr() [id] = x;
	simData->getPyPtr() [id] = y;
	// Set velocities
	simData->getVxPtr() [id] = velocity.x;
	simData->getVyPtr() [id] = velocity.y;
      }
      if (useOm) {
	simData->getOm(id) = omega;
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
  }

  void Insertion::modify(SimData* simData, int id, RealType dt) {
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

  void Circulate::modify(SimData* simData, int id, RealType dt) {
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

  void Fixed::modify(SimData* simData, int id, RealType dt) {
    simData->getFxPtr() [id] = 0;
    simData->getFyPtr() [id] = 0;
    simData->getVxPtr() [id] = 0;
    simData->getVyPtr() [id] = 0;
    simData->getPxPtr() [id] = pos.x;
    simData->getPyPtr() [id] = pos.y;
  }

  void ApplyForce::modify(SimData* simData, int id, RealType dt) {
    F += dt*dF;
    simData->getFxPtr() [id] = F.x;
    simData->getFyPtr() [id] = F.y;
  }

}
