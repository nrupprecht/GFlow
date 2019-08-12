#include "hertz_force__2d.hpp"

namespace GFlowSimulation {

  void HertzForce2d::kernel(int id1, int id2, RealType R1, RealType R2, RealType rsqr, RealType *dr, RealType **f) const {
    /*
    // Square root and inverse
    RealType r = sqrt(rsqr);
    RealType invr = 1./r;

    // Create a normal vector
    dr[0] *= invr;
    dr[1] *= invr;
    // Effective mass.
    RealType meff = 1./(1./im[id1] + 1./im[id2]);

    // Normal force.
    RealType Fn = sqrt(r/(R1+R2))* (K_n * r - 0.5*meff*gamma_n * Vn);
    // Tangential force.

    // Update forces
    fr = Fn * dr[0];
    f[id1][0] += fr;
    f[id2][0] -= fr;
    fr = Fn * dr[1];
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
  */
  }

}