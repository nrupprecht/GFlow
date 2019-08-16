#ifndef __HARD_SPHERE_HPP__GFLOW__
#define __HARD_SPHERE_HPP__GFLOW__

#include "../base/interaction.hpp"

#include "../utility/generic-dimension.hpp"

namespace GFlowSimulation {

  /**
  *  \brief Template class for hard sphere interactions.
  */
  template<int dims> class HardSphere : public Interaction {
  public:
    //! \brief Default constructor.
    HardSphere(GFlow *gflow) : Interaction(gflow), repulsion(DEFAULT_HARD_SPHERE_REPULSION) {};
  
    //! \brief Set the repulsion parameter.
    void setRepulsion(RealType r) { repulsion = r; }

    //! \brief Suggests a safe timescale given the minimum mass of a particle that has this interaction.
    virtual RealType suggest_timescale(RealType mass) const override {
      return 2*PI/sqrt(2*repulsion/mass);
    }

    virtual void interact() const override {
      // Common tasks
      Interaction::interact();

      // Do dimensional check.
      // \todo Should probably have some sort of global error message system.
      if (sim_dimensions!=dims) return;

      // Get the data pointers.
      RealType **x = Base::simData->X();
      RealType **f = Base::simData->F();
      RealType *sg = Base::simData->Sg();
      int    *type = Base::simData->Type();

      // Make sure all needed pointers are non null.
      // \todo Should probably have some sort of global error message system.
      if (x==nullptr || f==nullptr || sg==nullptr || type==nullptr) return;

      // Needed constants
      RealType sg1, sg2, rsqr, r, invr, magnitude, X[dims];

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
        sg1 = sg[id1];
        sg2 = sg[id2];
        // If close, interact.
        if (rsqr < sqr(sg1 + sg2)) {
          // Calculate distance, inverse distance.
          r = sqrt(rsqr);
          invr = 1./r;
          // Create a normal vector
          scalar_mult_eq_vec<dims>(X, invr);
          // Calculate the magnitude of the force
          magnitude = repulsion*(sg1 + sg2 - r);
          // Update forces
          sum_eq_vec_scaled<dims>(f[id1], X, magnitude);
          sum_eq_vec_scaled<dims>(f[id2], X, -magnitude);

          // Calculate potential
          if (do_potential) {
            potential += 0.5*repulsion*sqr(r - sg1 - sg2);
          }
          // Calculate virial
          if (do_virial) {
            virial += magnitude*r;
          }
        }
      }

      // --- Do verlet wrap part.
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
        sg1 = sg[id1];
        sg2 = sg[id2];
        // If close, interact.
        if (rsqr < sqr(sg1 + sg2)) {
          // Calculate distance, inverse distance.
          r = sqrt(rsqr);
          invr = 1./r;
          // Create a normal vector
          scalar_mult_eq_vec<dims>(X, invr);
          // Calculate the magnitude of the force
          magnitude = repulsion*(sg1 + sg2 - r);
          // Update forces
          sum_eq_vec_scaled<dims>(f[id1], X, magnitude);
          sum_eq_vec_scaled<dims>(f[id2], X, -magnitude);

          // Calculate potential
          if (do_potential) {
            potential += 0.5*repulsion*sqr(r - sg1 - sg2);
          }
          // Calculate virial
          if (do_virial) {
            virial += magnitude*r;
          }
        }
      }
    }

  protected:
    //! \brief The hard sphere repulsion parameter.
    RealType repulsion;
  };

}

#endif // __HARD_SPHERE_HPP__GFLOW__