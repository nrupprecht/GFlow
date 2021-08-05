#include <bonded/harmonicbond_2d.hpp>

using namespace GFlowSimulation;

HarmonicBond_2d::HarmonicBond_2d(GFlow *gflow)
    : HarmonicBond(gflow) {};

HarmonicBond_2d::HarmonicBond_2d(GFlow *gflow, RealType K)
    : HarmonicBond(gflow, K) {};

void HarmonicBond_2d::interact() const {
  // Call parent class.
  HarmonicBond::interact();
  // This object only acts in two dimensions.
  if (sim_dimensions != 2) {
    return;
  }
  // Get the number of bonds. left and right will have the same size - we checked this last time
  // updateLocalIDs was called.
  int nbonds = left.size();

  // Check if local ids need updating.
  if (simData->getNeedsRemake()) {
    updateLocalIDs();
  }

  // Get simdata, check if the local ids need updating
  auto x = simData->X();
  auto f = simData->F();
  auto type = simData->Type();
  // Displacement vector
  RealType dx[2], r, invr, dr, scale;
  // Get the bounds and boundary conditions
  Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
  BCFlag boundaryConditions[2];
  copyVec(Base::gflow->getBCs(), boundaryConditions, 2); // Keep a local copy of the bcs

  // Go through all bonds, calculating forces.
  for (int i = 0; i < nbonds; ++i) {
    // Get the global, then local ids of the particles.
    int id1 = left[i], id2 = right[i];
    // Check if the types are good
    if (type(id1) < 0 || type(id2) < 0) {
      continue;
    }
    // Calculate displacement
    dx[0] = x(id1, 0) - x(id2, 0);
    dx[1] = x(id1, 1) - x(id2, 1);
    gflow->minimumImage(dx);
    // Get the distance, inverse distance
    r = sqrtf(dotVec(dx, dx, 2));
    invr = 1. / r;
    // Calculate displacement from equilibrium
    dr = r - distance[i];
    // Makes dx into the force. The springConstant*dr comes from the force, the invr comes from normalizing dx
    scale = springConstant * dr * invr;
    dx[0] *= scale;
    dx[1] *= scale;
    // Add forces to particles
    f(id1, 0) -= dx[0];
    f(id1, 1) -= dx[1];
    f(id2, 0) += dx[0];
    f(id2, 1) += dx[1];

    // Potential energy
    if (do_potential) {
      potential += springConstant * sqr(dr);
    }
    if (do_virial) {
      virial += springConstant * dr * r;
    }
  }
}
