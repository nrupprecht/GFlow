#ifndef __GENERIC_HANDLERS_HPP__GFLOW__
#define __GENERIC_HANDLERS_HPP__GFLOW__

#include "../base/interaction.hpp"
#include "../utility/generic-dimension.hpp"

namespace GFlowSimulation {

  template<int dims, class ForceType> class VerletListPairs : public Interaction {
  public:
    //! \brief Constructor.
    VerletListPairs(GFlow *gflow) : Interaction(gflow) {};

    //! \brief Iterate through all 
    virtual void interact() const override {
      // Common tasks 
      Interaction::interact();

      // Do dimensional check.
      // \todo Should probably have some sort of global error message system.
      if (sim_dimensions!=dims) return;

      // Get the data pointers.
      RealType **x = simData->X();
      RealType **f = simData->F();
      RealType *rd = simData->Sg();
      int    *type = simData->Type();

      // Make sure all needed pointers are non null.
      // \todo Should probably have some sort of global error message system.
      if (x==nullptr || f==nullptr || rd==nullptr || type==nullptr) return;

      // Needed constants
      RealType R1, R2, rsqr, r, X[dims];

      // --- Go through all particles in verlet.
      for (int i=0; i<verlet.size(); i+=2) {
        int id1 = verlet[i];
        int id2 = verlet[i+1];
        // Check if the types are good
        if (type[id1]<0 || type[id2]<0) continue;
        // Calculate displacement.
        subtract_vec<dims>(x[id1], x[id2], X);
        // Calculate squared distance
        rsqr = dot_vec<dims>(X, X);
        // Get radii
        R1 = rd[id1];
        R2 = rd[id2];
        // If close, interact.
        if (rsqr < sqr((R1 + R2)*cutoff)) 
          force(id1, id2, R1, R2, rsqr, X, f);
      }

      // --- Go through all particles in verlet wrap.
      if (verlet_wrap.empty()) return;

      // Get the bounds and boundary conditions
      vector<BCFlag> bcs;
      vector<RealType> widths;
      // Fill the vectors.
      for (int d=0; d<dims; ++d) {
        bcs.push_back(gflow->getBC(d));
        widths.push_back(gflow->getBounds().wd(d));
      }

      // --- Go through all particles in verlet_wrap.
      for (int i=0; i<verlet_wrap.size(); i+=2) {
        int id1 = verlet_wrap[i];
        int id2 = verlet_wrap[i+1];
        // Check if the types are good
        if (type[id1]<0 || type[id2]<0) continue;
        // Calculate displacement.
        subtract_vec<dims>(x[id1], x[id2], X);
        // Harmonic corrections to distance.
        harmonic_correction<dims>(bcs.data(), X, widths.data());
        // Calculate squared distance
        rsqr = dot_vec<dims>(X, X);
        // Get radii
        R1 = rd[id1];
        R2 = rd[id2];
        // If close, interact.
        if (rsqr < sqr((R1 + R2)*cutoff)) 
          force(id1, id2, R1, R2, rsqr, X, f);
      }
    }

    virtual void force(const int id1, const int id2, const RealType R1, const RealType R2, const RealType rsqr, RealType *X, RealType **f) const = 0;
  };

}
#endif // __GENERIC_HANDLERS_HPP__GFLOW__