#include "twowallbinforce.hpp"
// Other files
#include "kineticenergydata.hpp"
#include "../../base/simdata.hpp"

namespace GFlowSimulation {

  TwoWallBinForce::TwoWallBinForce(GFlow *gflow, WallSlideBody *wa, WallSlideBody *wb) : GraphObject(gflow, "TwoWallBinForce", "x", "<F>"), wallA(wa), wallB(wb) {
    // Data
    data = vector<RPair>(nbins, RPair(0,0));
  }

  void TwoWallBinForce::pre_integrate() {
    // Clear data vector
    GraphObject::pre_integrate();
    // Clear data
    data = vector<RPair>(nbins, RPair(0,0));
    counts = vector<int>(nbins, 0);
  }

  void TwoWallBinForce::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check() || wallA==nullptr || wallB==nullptr) return;

    // Bin forces
    int n = Base::simData->number();

    // Compute the temperature of the system.
    int n_solvent = simData->number(); //  - wallA->size() = wallB->size();
    RealType KB = gflow->getKB();
    // Get KE per (non infinitely massive) particle.
    RealType ke = KineticEnergyData::calculate_kinetic(simData, true);
    RealType T = 2./static_cast<RealType>(sim_dimensions) * ke / KB;
    // Compute the number density. Do not include your own particles. Assumes the other line has roughly the same number of particles.
    RealType rho = n_solvent / gflow->getBounds().vol();
    // Compute the normalization
    RealType length = wallA->getLength();
    RealType norm = 1./(rho * KB * T * length);

    // Get distance between walls
    RealType x1 = wallA->getPosition();
    RealType x2 = wallB->getPosition();
    RealType dx = x2 - x1;

    // Correct for harmonic bc's
    RealType width = gflow->getBounds().wd(0);
    if (gflow->getBC(0)==BCFlag::WRAP && fabs(dx)>0.5*width)
      dx = dx>0 ? dx - 0.5*width : dx + 0.5*width;

    // Check cutoff
    if (min_distance <= dx && dx<max_distance) {
      // Distance bin
      int bx = static_cast<int>((dx-min_distance)/(max_distance-min_distance)*nbins);

      // Get the forces for left wall
      RealType F = wallA->getFnet()*sign(dx);
      data[bx].second += F*norm;
      // Get the forces for the right wall
      F = wallB->getFnet()*sign(dx);
      data[bx].second -= F*norm;

      // Increment counts
      counts[bx] += 2;
    }
  }

  bool TwoWallBinForce::writeToFile(string fileName, bool useName) {
    // Process data
    RealType dr = (max_distance - min_distance)/(nbins-1);
    for (int i=0; i<nbins; ++i) {
      data[i].first = (i+0.5)*dr + min_distance;
      if (counts[i]>0) data[i].second /= counts[i];
    }

    // Write to file
    GraphObject::writeToFile(fileName, useName);
  }

}