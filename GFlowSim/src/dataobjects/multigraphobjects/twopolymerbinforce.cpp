#include "twopolymerbinforce.hpp"
// Other files
#include "../graphobjects/kineticenergydata.hpp"
#include "../../base/simdata.hpp"
#include "../../bonded/angle-harmonic-chain.hpp"

namespace GFlowSimulation {

  TwoPolymerBinForce::TwoPolymerBinForce(GFlow *gflow) : MultiGraphObject(gflow, "TwoPolymerBinForce", "x", "<F>", 3) {
    // Data
    resetData(nbins);
    // Set label.
    axis_y[1] = "<F> chain";
    axis_y[2] = "counts";
  }

  TwoPolymerBinForce::TwoPolymerBinForce(GFlow *gflow, Group& ga, Group& gb) : MultiGraphObject(gflow, "TwoPolymerBinForce", "x", "<F>", 3), 
    polyA(ga), polyB(gb) 
  {
    // Data
    resetData(nbins);
    // Set label.
    axis_y[1] = "<F> chain";
    axis_y[2] = "counts";
  }

  void TwoPolymerBinForce::pre_integrate() {
    // Clear data vector
    MultiGraphObject::pre_integrate();
    // Clear data
    resetData(nbins);
  }

  void TwoPolymerBinForce::post_forces() {
    // Only record if enough time has gone by
    if (!DataObject::_check() || polyA.empty() || polyB.empty()) return;

    // Compute the temperature of the system.
    int n_solvent = simData->number();
    RealType KB = gflow->getKB();
    // Get KE per (non infinitely massive) particle.
    RealType ke = KineticEnergyData::calculate_kinetic(simData, true);
    RealType T = 2. * ke / (sim_dimensions * KB) ;
    // Compute the number density. Do not include your own particles. Assumes the other line has roughly the same number of particles.
    RealType rho = n_solvent / gflow->getBounds().vol();
    // Compute the normalization - HOPEFULLY sg(0) represents the size of solvent molecules.
    RealType norm = 1./(rho * KB * T);

    // Make sure local ids are up to date.
    if (locals_changed) {
      polyA.update_local_ids(simData);
      polyB.update_local_ids(simData);
    }
    locals_changed = false;
    // Find stats for both polymers.
    find_forces(0, norm);
    find_forces(1, norm);
  }

  void TwoPolymerBinForce::post_integrate() {
    // Only process if the run is in simulation mode.
    if (gflow->getRunMode()!=RunMode::SIM) return;

    // Process data
    RealType dr = (max_distance - min_distance)/nbins;
    // Process data
    for (int i=0; i<nbins; ++i) {
      // Set distance
      atX(i) = (i+0.5)*dr + min_distance;
      // Normalize forces
      if (atY(2, i)>0) {
        atY(0, i) /= atY(2, i);
        atY(1, i) /= atY(2, i);
      }
    }
  }

  void TwoPolymerBinForce::setMaxDistance(RealType md) {
    max_distance = md;
  }

  void TwoPolymerBinForce::setMinDistance(RealType md) {
    min_distance = md;
  }

  void TwoPolymerBinForce::setBins(int b) {
    if (b>0) {
      nbins = b;
      resetData(nbins);
    }
  }

  void TwoPolymerBinForce::setCType(int t) {
    c_type = t;
  }

  void TwoPolymerBinForce::setFirstPolymer(Group& g) {
    polyA = g;
  }

  void TwoPolymerBinForce::setFirstChain(AngleHarmonicChain *ch) {
    chainA = ch; 
  }

  void TwoPolymerBinForce::setSecondPolymer(Group& g) {
    polyB = g;
  }

  void TwoPolymerBinForce::setSecondChain(AngleHarmonicChain *ch) {
    chainB = ch; 
  }

  inline void TwoPolymerBinForce::find_forces(int choice, RealType norm) {
    RealType **x = simData->X();
    Vec X1(sim_dimensions), X2(sim_dimensions), dX(sim_dimensions), mdX(sim_dimensions);

    // Select objects
    Group& first  = (choice==0) ? polyA : polyB;
    Group& second = (choice==0) ? polyB : polyA;
    AngleHarmonicChain *firstCh  = (choice==0) ? chainA : chainB;
    // For holding the force on a particle.
    Vec force(sim_dimensions);

    // Check polymer A.
    for (int i=0; i<first.size(); ++i) {
      // Get id of particle from first group.
      int id1 = first.at(i), min_id = -1;
      // Only look at forces on the primary monomers, not the chain.
      if (simData->Type(id1)==c_type) continue;
      // Minimum distance to any particle in the other group. Start with a huge number
      RealType minD = 10000.;
      // Store position of first particle.
      X1 = x[id1];
      // Iterate through particles in second polymer
      for (int j=0; j<second.size(); ++j) {
        int id2 = second.at(j);
        // Copy vector
        X2 = x[id2];
        dX = X2 - X1;
        gflow->minimumImage(dX.data);
        RealType d = distance(dX);
        // Check if this is the new min
        if (d<minD) {
          minD = d;
          min_id = id2;
          mdX = dX;
        }
      }

      // Check cutoff
      if (min_distance <= minD && minD<max_distance) {
        // Distance bin
        RealType dr = (max_distance - min_distance)/nbins;
        int bx = static_cast<int>((minD-min_distance)/dr);
        // Get the net force on the particle.
        force = simData->F(id1);
        // If we can, subtract away the force of the polymer on itself.
        if (chainA && chainB) {
          Vec fr = firstCh->getForceByID(id1);
          force -= fr;
          atY(1, bx) += mdX*fr;
        }
        // Project force on particle from polymer A in the mdX direction.
        mdX.normalize();
        RealType F = mdX*force;
        atY(0, bx) += F*norm;
        // Increment counts
        ++atY(2, bx); // We added two force data points.
      }
    }
  }

}