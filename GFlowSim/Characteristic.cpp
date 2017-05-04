#include "Characteristic.h"
#include "Sectorization.h"

Bacteria::Bacteria() {
  orient = randV();
  reorient = 1./default_bacteria_reorient;
  strength = default_bacteria_strength;
  secretion = default_bacteria_secretion;
  velocity = default_bacteria_target_velocity;
  delay = 0.05;
  timer = drand48()*delay;
  fitness = 0.5;
  reproduction = 1./default_bacteria_reproduction_const;
  death = 1./default_bacteria_death_const;
}

Bacteria::Bacteria(const Bacteria& b) {
  orient = randV();
  reorient = b.reorient;
  strength = b.strength;
  secretion = b.secretion;
  velocity = b.velocity;
  delay = 0.05;
  timer = drand48()*delay;
  fitness = b.fitness;
  reproduction = b.reproduction;
  death = b.death;
}

void Bacteria::modify(double **pdata, Sectorization *sectors, int id) {
  // Set up convenience pointers
  double *px=pdata[0], *py=pdata[1], *vx=pdata[2], *vy=pdata[3], *fx=pdata[4], *fy=pdata[5], *th=pdata[6], *om=pdata[7], *tq=pdata[8], *sg=pdata[9], *im=pdata[10], *iI=pdata[11], *rp=pdata[12], *ds=pdata[13], *cf=pdata[14];

  // Update timer
  timer += sectors->getEpsilon();

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
    if (pr<drand48() && !sectors->isFull()) 
      sectors->insertParticle(Particle(px[id], py[id], sg[id]), this->create());
  }
  // Run
  double dvx = velocity*orient.x-vx[id], dvy = velocity*orient.y-vy[id];
  double mass = 1./im[id];
  fx[id] += strength*dvx*mass;
  fy[id] += strength*dvy*mass;
}

Characteristic* Bacteria::create() {
  return new Bacteria;
}

ConstantVelocity::ConstantVelocity() {
  // Default as a stationary object
  targetVelocity = 0;
  targetOmega = 0;
  useV = useOm = true;
}

ConstantVelocity::ConstantVelocity(const ConstantVelocity& cv) {
  targetVelocity = cv.targetVelocity;
  targetOmega    = cv.targetOmega;
  useV = cv.useV; useOm = cv.useOm;
}

ConstantVelocity::ConstantVelocity(vec2 v) : targetVelocity(v), targetOmega(0), useV(true), useOm(false) {};

ConstantVelocity::ConstantVelocity(vec2 v, double om) : targetVelocity(v), targetOmega(om), useV(true), useOm(true) {};

ConstantVelocity::ConstantVelocity(vec2 v, bool uv, double om, bool uo) : targetVelocity(v), targetOmega(om), useV(uv), useOm(uo) {};

void ConstantVelocity::modify(double **pdata, Sectorization*, int id) {
  double *px=pdata[0], *py=pdata[1], *vx=pdata[2], *vy=pdata[3], *fx=pdata[4], *fy=pdata[5], *th=pdata[6], *om=pdata[7], *tq=pdata[8], *sg=pdata[9], *im=pdata[10], *iI=pdata[11], *rp=pdata[12], *ds=pdata[13], *cf=pdata[14];
  // Keep the object moving at a constant velocity or omega
  if (useV) {
    vx[id] = targetVelocity.x;
    vy[id] = targetVelocity.y;
  }
  if (useOm) om[id] = targetOmega;
}

Characteristic* ConstantVelocity::create() {
  return new ConstantVelocity(targetVelocity);
}
