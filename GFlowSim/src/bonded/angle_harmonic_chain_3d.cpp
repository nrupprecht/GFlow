#include "angle_harmonic_chain_3d.hpp"
// Other files
#include "../utility/generic-dimension.hpp"

namespace GFlowSimulation {

  AngleHarmonicChain_3d::AngleHarmonicChain_3d(GFlow *gflow) : AngleHarmonicChain(gflow) {};

  AngleHarmonicChain_3d::AngleHarmonicChain_3d(GFlow *gflow, RealType K) : AngleHarmonicChain(gflow, K) {};

  void AngleHarmonicChain_3d::interact() const {
    // Call parent class.
    AngleHarmonicChain::interact();
    // The dimensionality is constant.
    constexpr int dims = 3;
    // Get arrays.
    auto x = simData->X();
    auto f = simData->F();
    // If there are not enough particles in the buffer
    if (local_ids.size()<3) return;
    // Variables and buffers
    RealType dx1, dy1, dz1, dx2, dy2, dz2, rsqr1, rsqr2, r1, r2, angle_cos, a11, a12, a22;
    RealType X1[dims], X2[dims];
    RealType f1[dims], f2[dims];

    // Get the bounds and boundary conditions.
    BCFlag bcs[dims];
    RealType widths[dims];
    // Fill the vectors.
    for (int d=0; d<dims; ++d) {
      bcs[d] = gflow->getBC(d);
      widths[d] = gflow->getBounds().wd(d);
    }

    // If there are fewer than three particles, this loop will not execute
    for (int i=0; i+2<local_ids.size(); ++i) {
      // Get the global, then local ids of the particles.
      int id1 = local_ids[i], id2 = local_ids[i+1], id3=local_ids[i+2];
      // First bond
      subtract_vec<dims>(x(id1), x(id2), X1);
      subtract_vec<dims>(x(id3), x(id2), X2);

      // Harmonic corrections to distance.
      harmonic_correction<dims>(bcs, X1, widths);
      harmonic_correction<dims>(bcs, X2, widths);

      // Get distance squared, distance
      rsqr1 = dot_vec<dims>(X1, X1);
      r1 = sqrt(rsqr1);
      rsqr2 = dot_vec<dims>(X2, X2);
      r2 = sqrt(rsqr2);

      // Find the cosine of the angle between the atoms
      //angle_cos = dx1*dx2 + dy1*dy2 + dz1*dz2;
      angle_cos = dot_vec<dims>(X1, X2);
      angle_cos /= (r1*r2);
      // Calculate the force
      a11 = angleConstant*angle_cos / rsqr1;
      a12 = -angleConstant / (r1*r2);
      a22 = angleConstant*angle_cos / rsqr2;
      // Force buffers - this is specific to 3 dimensions
      f1[0] = a11*X1[0] + a12*X2[0];
      f1[1] = a11*X1[1] + a12*X2[1];
      f1[2] = a11*X1[2] + a12*X2[2];
      f2[0] = a22*X2[0] + a12*X1[0];
      f2[1] = a22*X2[1] + a12*X1[1];
      f2[2] = a22*X2[2] + a12*X1[2];

      // First atom
      sum_eq_vec<dims>(f(id1), f1);
      // Second atom
      f(id2, 0) -= (f1[0] + f2[0]);
      f(id2, 1) -= (f1[1] + f2[1]);
      f(id2, 2) -= (f1[2] + f2[2]);
      // Third atom
      sum_eq_vec<dims>(f(id2), f2);
    }
    // \todo For two atoms, there is no angle. Just compute the bond force.
  }

}