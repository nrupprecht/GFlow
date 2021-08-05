#include <bonded/angle_harmonic_chain_2d.hpp>

using namespace GFlowSimulation;

AngleHarmonicChain_2d::AngleHarmonicChain_2d(GFlow *gflow)
    : AngleHarmonicChain(gflow) {};

AngleHarmonicChain_2d::AngleHarmonicChain_2d(GFlow *gflow, RealType K)
    : AngleHarmonicChain(gflow, K) {};

void AngleHarmonicChain_2d::interact() const {
  // Call parent class.
  AngleHarmonicChain::interact();
  // The dimensionality is constant.
  constexpr int dims = 2;
  // Get arrays
  auto x = simData->X();
  auto f = simData->F();

  // If there are not enough particles in the buffer
  if (local_ids.size() < 3) {
    return;
  }
  // Variables and buffers
  RealType rsqr1, rsqr2, r1, r2, angle_cos;
  real dx1[dims], dx2[dims], f1[dims], f3[dims], pa[dims], pb[dims];

  // Get the bounds and boundary conditions.
  BCFlag bcs[dims];
  RealType widths[dims];
  // Fill the vectors.
  for (int d = 0; d < dims; ++d) {
    bcs[d] = gflow->getBC(d);
    widths[d] = gflow->getBounds().wd(d);
  }

  // Set forces to be the negative of their current value, so we can add the resulting forces at the end.
  for (int i = 0; i < local_ids.size(); ++i) {
    forces[i][0] = -f(local_ids[i], 0);
    forces[i][1] = -f(local_ids[i], 1);
  }

  // If there are fewer than three particles, this loop will not execute
  for (int i = 0; i + 2 < local_ids.size(); ++i) {
    // Get the global, then local ids of the particles.
    int id1 = local_ids[i], id2 = local_ids[i + 1], id3 = local_ids[i + 2];
    // id1     id2     id3
    //  A <---- B ----> C
    //     dx1    dx2
    //
    // First bond
    if (i == 0) {
      // dx1[0] = x[id1][0] - x[id2][0];
      // dx1[1] = x[id1][1] - x[id2][1];
      subtract_vec<dims>(x(id1), x(id2), dx1);
      // Minimum image under harmonic b.c.s
      //gflow->minimumImage(dx1.data);
      harmonic_correction<dims>(bcs, dx1, widths);
      // Get distance squared, distance
      rsqr1 = dot_vec<dims>(dx1, dx1);
      // rsqr1 = dx1[0]*dx1[0] + dx1[1]*dx1[1];
      r1 = sqrt(rsqr1);
      // Normalize dx1, dx2
      dx1[0] /= r1;
      dx1[1] /= r1;
    }
    else {
      dx1[0] = -dx2[0];
      dx1[1] = -dx2[1];
      r1 = r2;
    }
    // Second bond
    dx2[0] = x(id3, 0) - x(id2, 0);
    dx2[1] = x(id3, 1) - x(id2, 1);
    // Minimum image under harmonic b.c.s
    //gflow->minimumImage(dx2.data);
    harmonic_correction<dims>(bcs, dx2, widths);

    // Get distance squared, distance
    rsqr2 = dx2[0] * dx2[0] + dx2[1] * dx2[1];
    r2 = sqrt(rsqr2);
    // Normalize dx2
    dx2[0] /= r2;
    dx2[1] /= r2;
    // Find the cosine of the angle between the atoms
    angle_cos = (dx1[0] * dx2[0] + dx1[1] * dx2[1]);

    // If angle_cos is too large or small (probably due to inexactness of floating point), we can get nan-s.
    RealType theta = acos(0.999 * angle_cos);
    RealType dTh = (theta - PI);
    //RealType dTh = -1. - angle_cos;

    tripleProduct2(dx1, dx1, dx2, pa);
    tripleProduct2(dx2, dx1, dx2, pb);
    negate_vec<dims>(pb);
//    // Normalize pa, pb
//    RealType pa_norm = sqrtf(pa[0] * pa[0] + pa[1] * pa[1]);
//    RealType pb_norm = sqrtf(pb[0] * pb[0] + pb[1] * pb[1]);

    if (!normalize_vec<dims>(pa) || !normalize_vec<dims>(pb)) {
      continue;
    }

    /*
    // In case vectors (dx1, dx2) are very close to being parallel or antiparallel.
    if (pa_norm==0 || pb_norm==0) continue;
    // Normalize pa, pb
    pa[0] /= pa_norm; pa[1] /= pa_norm;
    pb[0] /= pb_norm; pb[1] /= pb_norm;
    */

    // Calculate force strength due to angle
    RealType str = -angleConstant * dTh;
    RealType strength1 = str / r1;
    RealType strength2 = str / r2;

    // Calculate force strength due to harmonic bond A-B
    RealType dr = r1 - distance[i];
    RealType repl = springConstant * dr;

    // Set force buffers
    f1[0] = strength1 * pa[0];
    f1[1] = strength1 * pa[1];
    f3[0] = strength2 * pb[0];
    f3[1] = strength2 * pb[1];

    // First atom
    f(id1, 0) += (f1[0] - repl * dx1[0]);
    f(id1, 1) += (f1[1] - repl * dx1[1]);
    // Second atom
    f(id2, 0) -= (f1[0] + f3[0] - repl * dx1[0]);
    f(id2, 1) -= (f1[1] + f3[1] - repl * dx1[1]);
    // Third atom
    f(id3, 0) += f3[0];
    f(id3, 1) += f3[1];

    if (do_potential) {
      potential += 0.5 * angleConstant * sqr(theta - PI) + springConstant * sqr(dr);
    }

    // If this is the last angle, compute the B-C harmonic bond as well.
    if (i + 3 == local_ids.size()) {
      dr = r2 - distance[i + 1];
      repl = springConstant * dr;
      // Second atom
      f(id2, 0) += repl * dx2[0];
      f(id2, 1) += repl * dx2[1];
      // Third atom
      f(id3, 0) -= repl * dx2[0];
      f(id3, 1) -= repl * dx2[1];
      // Extra potential from the last bond.
      if (do_potential) {
        potential += springConstant * sqr(dr);
      }
    }

  }

  // Add forces to buffer
  for (int i = 0; i < local_ids.size(); ++i) {
    forces[i][0] += f(local_ids[i], 0);
    forces[i][1] += f(local_ids[i], 1);
  }

  // For two atoms, there is no angle. Just compute the bond force.
}
