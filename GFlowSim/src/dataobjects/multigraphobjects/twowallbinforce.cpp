#include "twowallbinforce.hpp"
// Other files
#include "../graphobjects/kineticenergydata.hpp"
#include "../../base/simdata.hpp"

namespace GFlowSimulation {

  TwoWallBinForce::TwoWallBinForce(GFlow *gflow, WallSlideBody *wa, WallSlideBody *wb) : MultiGraphObject(gflow, "TwoWallBinForce", "x", "<F>", 2), wallA(wa), wallB(wb) {
    // Data
    resetData(nbins);
    // Set counts to not be written
    write_data[1] = false;
  }

  void TwoWallBinForce::pre_integrate() {
    // Clear data vector
    MultiGraphObject::pre_integrate();
    // Clear data
    resetData(nbins);
  }

  void TwoWallBinForce::post_forces() {
    // Only record if enough time has gone by
    if (!DataObject::_check() || wallA==nullptr || wallB==nullptr) return;

    // Compute the temperature of the system.
    int n_solvent = simData->number(); //  - wallA->size() = wallB->size();
    RealType KB = gflow->getKB();
    // Get KE per (non infinitely massive) particle.
    RealType ke = KineticEnergyData::calculate_kinetic(simData, true);
    RealType T = 2./static_cast<RealType>(sim_dimensions) * ke / KB;
    // Compute the number density. Do not include your own particles. Assumes the other line has roughly the same number of particles.
    RealType rho = n_solvent / gflow->getBounds().vol();
    // Compute the normalization - HOPEFULLY sg(0) represents the size of solvent molecules.
    RealType length = wallA->getLength() + 2*simData->Sg(0);
    RealType norm = 1./(rho * KB * T * length);

    // Get distance between walls
    RealType x1 = wallA->getPosition();
    RealType x2 = wallB->getPosition();
    RealType dx = x2 - x1;
    // Correct for harmonic bc's
    gflow->minimumImage(dx, 0);

    // Check cutoff
    if (min_distance <= dx && dx<max_distance) {
      // Distance bin
      RealType dr = (max_distance - min_distance)/nbins;
      int bx = static_cast<int>((dx-min_distance)/dr);

      // Get the forces for left wall
      RealType F = wallA->getFnet()*sign(dx);
      atY(0, bx) += F*norm;
      // Get the forces for the right wall
      F = -wallB->getFnet()*sign(dx);
      atY(0, bx) += F*norm;

      // Increment counts
      atY(1, bx) += 2; // We added two force data points.
    }
  }

  void TwoWallBinForce::post_integrate() {
    // Only process if the run is in simulation mode.
    if (gflow->getRunMode()!=RunMode::SIM) return;

    // Process data
    RealType dr = (max_distance - min_distance)/nbins;
    // Process data
    for (int i=0; i<nbins; ++i) {
      // Set distance
      atX(i) = (i+0.5)*dr + min_distance;
      // Normalize forces
      if (atY(1, i)>0) atY(0, i) /= atY(1, i);
      else atY(0, i) = 0;
    }
  }

  void TwoWallBinForce::setMaxDistance(RealType md) {
    max_distance = md;
  }

  void TwoWallBinForce::setMinDistance(RealType md) {
    min_distance = md;
  }

  void TwoWallBinForce::setBins(int b) {
    if (b>0) {
      nbins = b;
      resetData(nbins);
    }
  }

}