#ifndef __HARD_SPHERE_CF_HPP__GFLOW__
#define __HARD_SPHERE_CF_HPP__GFLOW__

#include "hard_sphere.hpp"

#include "../utility/simd_generic.hpp" // For un_clamp
#include "../integrators/angularvelocityverlet2d.hpp"

namespace GFlowSimulation {

  /** 
  *  \brief Template class for dissipative hard sphere interactions.
  */
  template<int dims> class HardSphereCf : public HardSphere<2> {
  public:
    //! \brief Default constructor.
    HardSphereCf(GFlow *gflow) 
      : HardSphere<dims>(gflow), dissipation(DEFAULT_HARD_SPHERE_DISSIPATION), angular_dissipation(DEFAULT_HARD_SPHERE_DISSIPATION), mu(0.25) {

      if (dims==2) {
        simData->requestScalarData("Om");
        simData->requestScalarData("Tq");
        // Add an angular integrator.
        gflow->addIntegrator(new AngularVelocityVerlet2d(gflow));
      }

    };

    //! \brief Set the dissipation constant.
    void setDissipation(RealType d) { dissipation = d>=0 ? d : dissipation; }

    //! \brief Set the dissipation constant.
    void setAngularDissipation(RealType d) { angular_dissipation = d>=0 ? d : angular_dissipation; }

    //! \brief Set the coefficient of friction.
    void setMu(RealType c) { mu = c>=0 ? c : mu; }

    virtual void interact() const override {
      // Common tasks - DONT call HardSpher::interact().
      Interaction::interact();

      // Do dimensional check.
      // \todo Should probably have some sort of global error message system.
      if (sim_dimensions!=dims) return;

      // Get the data pointers.
      RealType **x = simData->X();
      RealType **v = simData->V();
      RealType **f = simData->F();
      RealType *sg = simData->Sg();
      int om_add = simData->getScalarData("Om");
      int tq_add = simData->getScalarData("Tq");
      RealType *om = simData->ScalarData(om_add);
      RealType *tq = simData->ScalarData(tq_add);
      int    *type = simData->Type();

      // Make sure all needed pointers are non null.
      // \todo Should probably have some sort of global error message system.
      if (x==nullptr || f==nullptr || sg==nullptr || type==nullptr) return;

      // Needed constants
      RealType sg1, sg2, rsqr, r, invr, X[dims], V[dims], Vt[dims], Vn[dims];

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

          // Calculate relative velocity.
          subtract_vec<dims>(v[id2], v[id1], V);
          
          // Calculate normal velocity.
          RealType vn = dot_vec<dims>(V, X);
          // Repulsion + dissipation, which only occurs on loading the spring.
          RealType Fn = repulsion*(sg1 + sg2 - r) + dissipation * un_clamp(vn);

          // Total velocity minus velocity in the normal direction is the tangential velocity.
          subtract_vec<dims>(V, Vn, Vt);
          // Magnitude of the tangential velocity, and its inverse.
          RealType vt = sqrt(dot_vec<dims>(Vt, Vt));
          RealType invvt = 1./vt;

          // Normalize tangential velocity.
          scalar_mult_eq_vec<dims>(Vt, invvt);
          
          RealType Ft = mu*Fn;

          // Angular velocities.
          // if (dims==2) {
          //   // Calculate difference in velocities at the point of intersection.
          //   RealType VT = vt - om[id1]*sg[id1] - om[id2]*sg[id2]; // Relative surface velocity due to angular velocity.
          // }

          // Update forces - normal force.
          sum_eq_vec_scaled<dims>(f[id1], X, Fn);
          sum_eq_vec_scaled<dims>(f[id1], Vt, Ft);
          // Update forces - tangential force.
          sum_eq_vec_scaled<dims>(f[id2], X, -Fn);
          sum_eq_vec_scaled<dims>(f[id2], Vt, -Ft);

          // Update torques
          if (dims==2) {
            tq[id1] += Ft*sg[id1];
            tq[id2] += Ft*sg[id2];
          }

          // Calculate potential
          if (do_potential) {
            potential += 0.5*repulsion*sqr(r - sg1 - sg2);
          }
          // Calculate virial
          if (do_virial) {
            virial += Fn*r;
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

          // Calculate relative velocity.
          subtract_vec<dims>(v[id2], v[id1], V);
          
          // Calculate normal velocity.
          RealType vn = dot_vec<dims>(V, X);
          // Repulsion + dissipation, which only occurs on loading the spring.
          RealType Fn = repulsion*(sg1 + sg2 - r) + dissipation * un_clamp(vn);

          // Total velocity minus velocity in the normal direction is the tangential velocity.
          subtract_vec<dims>(V, Vn, Vt);
          // Magnitude of the tangential velocity, and its inverse.
          RealType vt = sqrt(dot_vec<dims>(Vt, Vt));
          RealType invvt = 1./vt;

          // Normalize tangential velocity.
          scalar_mult_eq_vec<dims>(Vt, invvt);

          RealType Ft = mu*Fn;

          // Angular velocities.
          // if (dims==2) {
          //   // Calculate difference in velocities at the point of intersection.
          //   RealType VT = vt - om[id1]*sg[id1] - om[id2]*sg[id2]; // Relative surface velocity due to angular velocity.
          // }

          // Update forces - normal force.
          sum_eq_vec_scaled<dims>(f[id1], X, Fn);
          sum_eq_vec_scaled<dims>(f[id1], Vt, Ft);
          // Update forces - tangential force.
          sum_eq_vec_scaled<dims>(f[id2], X, -Fn);
          sum_eq_vec_scaled<dims>(f[id2], Vt, -Ft);

          // Update torques
          if (dims==2) {
            tq[id1] += Ft*sg[id1];
            tq[id2] += Ft*sg[id2];
          }

          // Calculate potential
          if (do_potential) {
            potential += 0.5*repulsion*sqr(r - sg1 - sg2);
          }
          // Calculate virial
          if (do_virial) {
            virial += Fn*r;
          }
        }
      }
    }

  protected:
    //! \brief The dissipation constant.
    RealType dissipation;

    //! \brief The angular dissipation constant.
    RealType angular_dissipation;

    //! \brief The coefficient of friction.
    RealType mu;
  };

}

#endif // __HARD_SPHERE_DS_HPP__GFLOW__