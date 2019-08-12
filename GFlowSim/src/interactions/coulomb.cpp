#include "coulomb.hpp"

namespace GFlowSimulation {

  Coulomb2d::Coulomb2d(GFlow *gflow) : Interaction2d(gflow), repulsion(0.0025) {
    cutoff = 5;
  };

  void Coulomb2d::setCutoff(RealType cut) {
    // Change cutoff.
    cutoff = cut>0 ? cut : cutoff;
    // Calculate PE shift.
    potential_energy_shift = repulsion/cutoff;
  }

  void Coulomb2d::kernel(int id1, int id2, RealType R1, RealType R2, RealType rsqr, RealType *dr, RealType **f) const {
    // Square root and inverse
    RealType r = sqrt(rsqr);
    RealType invr = 1./r;

    // Create a normal vector
    dr[0] *= invr;
    dr[1] *= invr;

    // Calculate the magnitude of the force
    RealType magnitude = repulsion/rsqr - repulsion/sqr(cutoff), fr;

    // Update forces
    fr = magnitude * dr[0];
    f[id1][0] += fr;
    f[id2][0] -= fr;
    fr = magnitude * dr[1];
    f[id1][1] += fr;
    f[id2][1] -= fr;

    // Calculate potential
    if (do_potential) {
      potential += repulsion/r - potential_energy_shift;
    }
    // Calculate virial
    if (do_virial) {
      virial += magnitude*r;
    }
  }

}