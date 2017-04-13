#include "Characteristic.h"
#include "Sectorization.h"

Bacteria::Bacteria() {
  orient = randV();
  reorient = 1./default_bacteria_reorient;
  strength = default_bacteria_strength;
  delay = 0.05;
  timer = drand48()*delay;
}

void Bacteria::modify(double& vx, double& vy, double& fx, double& fy, double& torque, Sectorization *sectors) {
  // (Possibly) reorient
  if (delay<timer) {
    timer = 0;
    double p = exp(-delay*reorient);
    if (p<drand48()) orient = randV();
  }
  // Run
  if (sqr(vx)+sqr(vy)<default_bacteria_target_velocity_sqr) {
    fx += strength*orient.x;
    fy += strength*orient.y;
  }
  // Update timer
  timer += sectors->getEpsilon();
}

Characteristic* Bacteria::create() {
  return new Bacteria;
}
