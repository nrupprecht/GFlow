#include "Characteristic.h"
#include "Sectorization.h"

Bacteria::Bacteria() {
  orient = randV();
  reorient = 1./default_bacteria_reorient;
  strength = default_bacteria_strength;
  secretion = default_bacteria_secretion;
  delay = 0.05;
  timer = drand48()*delay;
  fitness = 0.5;
  reproduction = 1./default_bacteria_reproduction_const;
  death = 1./default_bacteria_death_const;
}

void Bacteria::modify(double **pdata, Sectorization *sectors, int id) {
  // Set up convenience pointers
  double *px=pdata[0], *py=pdata[1], *vx=pdata[2], *vy=pdata[3], *fx=pdata[4], *fy=pdata[5], *th=pdata[6], *om=pdata[7], *tq=pdata[8], *sg=pdata[9], *im=pdata[10], *iI=pdata[11], *rp=pdata[12], *ds=pdata[13], *cf=pdata[14];
  // (Possibly) die because of lack of fitness
  if (fitness<0) {
    sectors->removeAt(id);
    return;
  }
  // (Possibly) reorient or reproduce or die
  if (delay<timer) {
    timer = 0;
    // Die
    static double gg = exp(-delay*death);
    if (gg<drand48()) {
      sectors->removeAt(id);
      return;
    }
    // Reorientation
    static double po = exp(-delay*reorient);
    if (po<drand48()) orient = randV();
    // Reproduction
    static double pr = exp(-delay*reproduction);
    if (pr<drand48()) {
      if (!sectors->isFull()) sectors->insertParticle(Particle(px[id], py[id], sg[id]), this->create());
    }
  }
  // Run
  if (sqr(vx[id])+sqr(vy[id])<default_bacteria_target_velocity_sqr) {
    fx[id] += strength*orient.x;
    fy[id] += strength*orient.y;
  }
  // Update timer
  timer += sectors->getEpsilon();
}

Characteristic* Bacteria::create() {
  return new Bacteria;
}

ConstantVelocity::ConstantVelocity(vec2 v) : targetVelocity(v), strength(10.) {};

void ConstantVelocity::modify(double **pdata, Sectorization*, int id) {
  double *px=pdata[0], *py=pdata[1], *vx=pdata[2], *vy=pdata[3], *fx=pdata[4], *fy=pdata[5], *th=pdata[6], *om=pdata[7], *tq=pdata[8], *sg=pdata[9], *im=pdata[10], *iI=pdata[11], *rp=pdata[12], *ds=pdata[13], *cf=pdata[14];
  // Apply a force to keep the object moving at a constant velocity
  fx[id] += strength*(targetVelocity.x-vx[id])/im[id];
  fy[id] += strength*(targetVelocity.y-vy[id])/im[id];
}

Characteristic* ConstantVelocity::create() {
  return new ConstantVelocity(targetVelocity);
}
