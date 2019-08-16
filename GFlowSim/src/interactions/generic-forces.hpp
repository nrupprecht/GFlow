#ifndef __HARD_SPHERE_GENERIC_HPP__GFLOW__
#define __HARD_SPHERE_GENERIC_HPP__GFLOW__

#include "generic-handlers.hpp"

namespace GFlowSimulation {

  template<int dims, template<int, class> class handler> class HardSphereGeneric : public handler<dims, HardSphereGeneric<dims, handler> > {
  public:
    //! \brief Constructor.
    //HardSphereGeneric(GFlow *gflow) : VerletListPairs<dims, HardSphereGeneric<dims> >(gflow) {};
    HardSphereGeneric(GFlow *gflow) : handler<dims, HardSphereGeneric<dims, handler> >(gflow) {};

    //! \brief Constructor, sets repulsion
    //HardSphereGeneric(GFlow *gflow, RealType rp) : VerletListPairs<dims, HardSphereGeneric<dims> >(gflow), repulsion(rp) {};
    HardSphereGeneric(GFlow *gflow, RealType rp) : handler<dims, HardSphereGeneric<dims, handler> >(gflow), repulsion(rp) {};

    //! \brief The force kernel for hard sphere forces.
    virtual void force(const int id1, const int id2, const RealType R1, const RealType R2, const RealType rsqr, RealType *X, RealType **f) const override {
      // Calculate distance, inverse distance.
      RealType r = sqrt(rsqr);
      RealType invr = 1./r;
      // Create a normal vector
      scalar_mult_eq_vec<dims>(X, invr);

      // Calculate the magnitude of the force
      RealType Fn = repulsion*(R1 + R2 - r);
      // Update forces
      sum_eq_vec_scaled<dims>(f[id1], X,  Fn);
      sum_eq_vec_scaled<dims>(f[id2], X, -Fn);

      // Calculate potential
      if (Interaction::do_potential) {
        Interaction::potential += 0.5*repulsion*sqr(r - R1 - R2);
      }
      // Calculate virial
      if (Interaction::do_virial) {
        Interaction::virial += Fn*r;
      }
    }

    virtual RealType suggest_timescale(RealType mass) const override {
      return 2*PI/sqrt(2*repulsion/mass);
    }

  private:
    //! \brief HThe strength of the force.
    RealType repulsion = DEFAULT_HARD_SPHERE_REPULSION;
  };

  // --- Define specific forces

  //! \brief Define the HardSphereVLP class to be hard sphere force using verlet list pairs container.
  template<int dims> 
  using HardSphereVLP = HardSphereGeneric<dims, VerletListPairs>;


  #include "../utility/simd_generic.hpp" // For un_clamp.

  /*
  * \brief Generic class for hard sphere forces with dissipation.
  *
  */
  template<int dims, template<int, class> class handler> class HardSphereDsGeneric : public handler<dims, HardSphereDsGeneric<dims, handler> > {
  public:
    //! \brief Constructor.
    HardSphereDsGeneric(GFlow *gflow) : handler<dims, HardSphereDsGeneric<dims, handler> >(gflow) {};

    //! \brief Constructor, sets repulsion
    HardSphereDsGeneric(GFlow *gflow, RealType rp) : handler<dims, HardSphereDsGeneric<dims, handler> >(gflow), repulsion(rp) {};

    //! \brief The force kernel for hard sphere forces.
    virtual void force(const int id1, const int id2, const RealType R1, const RealType R2, const RealType rsqr, RealType *X, RealType **f) const override {
      // Calculate distance, inverse distance.
      RealType r = sqrt(rsqr);
      RealType invr = 1./r;
      
      // Create a normal vector
      scalar_mult_eq_vec<dims>(X, invr);
      
      // Calculate relative velocity.
      RealType V[dims];
      subtract_vec<dims>(Base::simData->V(id2), Base::simData->V(id1), V);
      // Calculate normal velocity.
      RealType vn = dot_vec<dims>(V, X);
      // Dissipation only occurs on loading the spring.
      //Fn += dissipation * un_clamp(vn);
      // Calculate the magnitude of the force
      RealType Fn = max(RealType(0), repulsion*(R1 + R2 - r) + dissipation*vn);
      // Update forces
      sum_eq_vec_scaled<dims>(f[id1], X,  Fn);
      sum_eq_vec_scaled<dims>(f[id2], X, -Fn);

      // Calculate potential
      if (Interaction::do_potential) {
        Interaction::potential += 0.5*repulsion*sqr(r - R1 - R2);
      }
      // Calculate virial
      if (Interaction::do_virial) {
        Interaction::virial += Fn*r;
      }
    }

    virtual RealType suggest_timescale(RealType mass) const override {
      return 2*PI/sqrt(2*repulsion/mass);
    }

  private:
    //! \brief The strength of the force.
    RealType repulsion = DEFAULT_HARD_SPHERE_REPULSION;

    //! \brief The dissipation constant.
    RealType dissipation = DEFAULT_HARD_SPHERE_DISSIPATION;
  };

  // --- Define specific forces

  //! \brief Define the HardSphereVLP class to be hard sphere force using verlet list pairs container.
  template<int dims> 
  using HardSphereDsVLP = HardSphereDsGeneric<dims, VerletListPairs>;


  /*
  * \brief Generic class for lennard jones forces.
  *
  */
  template<int dims, template<int, class> class handler> class LennardJonesGeneric : public handler<dims, LennardJonesGeneric<dims, handler> > {
  public:
    //! \brief Constructor.
    //LennardJonesGeneric(GFlow *gflow) : VerletListPairs<dims, LennardJonesGeneric<dims> >(gflow) {
    LennardJonesGeneric(GFlow *gflow) : handler<dims, LennardJonesGeneric<dims, handler> >(gflow) {
      setCutoff(2.5);
    };

    //! \brief Constructor, sets strength.
    //LennardJonesGeneric(GFlow *gflow, RealType str) : VerletListPairs<dims, LennardJonesGeneric<dims> >(gflow), strength(str) {
    LennardJonesGeneric(GFlow *gflow, RealType str) : handler<dims, LennardJonesGeneric<dims, handler> >(gflow), strength(str) {
      setCutoff(2.5);
    };

    //! \brief Set the cutoff factor.
    void setCutoff(RealType c) { 
      Interaction::cutoff = c>0 ? c : Interaction::cutoff; 
      // Calculate the potential when r is the cutoff radius
      RealType gamma = 1./Interaction::cutoff;
      RealType g3  = gamma*gamma*gamma; 
      RealType g6  = g3*g3;
      RealType g12 = g6*g6;
      // Set the potential energy shift
      potential_energy_shift = 4.*strength*(g12 - g6);
    }

    virtual void force(const int id1, const int id2, const RealType R1, const RealType R2, const RealType rsqr, RealType *X, RealType **f) const override {
      // Calculate distance, inverse distance.
      RealType r = sqrt(rsqr);
      RealType invr = 1./r;
      // Create a normal vector
      scalar_mult_eq_vec<dims>(X, invr);

      // Calculate the magnitude of the force.
      RealType gamma = (R1+R2)*invr;
      RealType g3  = gamma*gamma*gamma; 
      RealType g6  = g3*g3;
      RealType g12 = g6*g6;
      // Calculate magnitude
      RealType Fn = 24.*strength*(2.*g12 - g6)*invr;
      // Update forces
      sum_eq_vec_scaled<dims>(f[id1], X,  Fn);
      sum_eq_vec_scaled<dims>(f[id2], X, -Fn);

      // Calculate potential
      if (Interaction::do_potential) {
        Interaction::potential += 4.*strength*(g12 - g6) - potential_energy_shift;
      }
      // Calculate virial
      if (Interaction::do_virial) {
        Interaction::virial += Fn*r;
      }      
    }

  private:
    //! \brief The strength of the force.
    RealType strength = DEFAULT_LENNARD_JONES_STRENGTH;

    //! \brief Energy shift caused by cutoff.
    RealType potential_energy_shift = 0;
  };

  // --- Define specific forces

  //! \brief Define the LennardJonesVLP class to be lennard jones force using verlet list pairs container.
  template<int dims> 
  using LennardJonesVLP = LennardJonesGeneric<dims, VerletListPairs>;


  /*
  * \brief Generic class for (short range) coulombic forces.
  *
  */
  template<int dims, template<int, class> class handler> class CoulombGeneric : public handler<dims, CoulombGeneric<dims, handler> > {
  public:
    //! \brief Constructor.
    CoulombGeneric(GFlow *gflow) : handler<dims, CoulombGeneric<dims, handler> >(gflow) {
      setCutoff(3.);
    };

    //! \brief Constructor, sets strength.
    CoulombGeneric(GFlow *gflow, RealType str) : handler<dims, CoulombGeneric<dims, handler> >(gflow), strength(str) {
      setCutoff(3.);
    };

    //! \brief Set the cutoff factor.
    void setCutoff(RealType c) { 
      Interaction::cutoff = c>0 ? c : Interaction::cutoff; 
      // Calculate the potential when r is the cutoff radius
      potential_energy_shift = strength/Interaction::cutoff;
    }

    virtual void force(const int id1, const int id2, const RealType R1, const RealType R2, const RealType rsqr, RealType *X, RealType **f) const override {
      // Calculate distance, inverse distance.
      RealType r = sqrt(rsqr);
      RealType invr = 1./r;
      // Create a normal vector
      scalar_mult_eq_vec<dims>(X, invr);

      // Calculate the magnitude of the force.
      RealType Fn = strength/rsqr - strength/sqr(Interaction::cutoff);

      // Update forces
      sum_eq_vec_scaled<dims>(f[id1], X,  Fn);
      sum_eq_vec_scaled<dims>(f[id2], X, -Fn);

      // Calculate potential
      if (Interaction::do_potential) {
        Interaction::potential += strength/r - potential_energy_shift;
      }
      // Calculate virial
      if (Interaction::do_virial) {
        Interaction::virial += Fn*r;
      }      
    }

  private:
    //! \brief The strength of the force.
    RealType strength = 0.0025;

    //! \brief Energy shift caused by cutoff.
    RealType potential_energy_shift = 0;
  };

  // --- Define specific forces

  //! \brief Define the LennardJonesVLP class to be lennard jones force using verlet list pairs container.
  template<int dims> 
  using CoulombVLP = CoulombGeneric<dims, VerletListPairs>;



  /*
  * \brief Generic class for (short range) coulombic forces.
  *
  */
  template<int dims, template<int, class> class handler> class BuckinghamGeneric : public handler<dims, BuckinghamGeneric<dims, handler> > {
  public:
    //! \brief Constructor.
    BuckinghamGeneric(GFlow *gflow) : handler<dims, BuckinghamGeneric<dims, handler> >(gflow) {
      setCutoff(3.);
    };

    //! \brief Constructor, sets strength.
    BuckinghamGeneric(GFlow *gflow, RealType str) : handler<dims, BuckinghamGeneric<dims, handler> >(gflow), strength(str) {
      setCutoff(3.);
    };

    //! \brief Set the cutoff factor.
    void setCutoff(RealType c) { 
      Interaction::cutoff = c>0 ? c : Interaction::cutoff;
      // Calculate the potential when r is the cutoff radius
      potential_energy_shift = 0; //\todo Calculate this.
    }

    virtual void force(const int id1, const int id2, const RealType R1, const RealType R2, const RealType rsqr, RealType *X, RealType **f) const override {
      // Calculate distance, inverse distance.
      RealType r = sqrt(rsqr);
      RealType invr = 1./r;
      // Create a normal vector
      scalar_mult_eq_vec<dims>(X, invr);

      // Calculate the magnitude of the force
      RealType radius = (R1+R2)*invr;
      RealType exp1 = expf(-ratio/radius);
      RealType sigma2 = radius*radius;
      RealType sigma6 = sigma2*sigma2*sigma2;
      RealType Fn = strength*invr*(ratio*exp1 - 6*sigma6);
      // Apply cutoff
      //Fn = min(inner_force, Fn);

      // Apply cutoff
      Fn = Fn<inner_force ? inner_force : Fn;

      // Update forces
      sum_eq_vec_scaled<dims>(f[id1], X,  Fn);
      sum_eq_vec_scaled<dims>(f[id2], X, -Fn);

      // Calculate potential
      if (Interaction::do_potential) {
        Interaction::potential += strength*(exp1 - sigma6) - potential_energy_shift;
      }
      // Calculate virial
      if (Interaction::do_virial) {
        Interaction::virial += Fn*r;
      }      
    }

  private:
    //! \brief The strength of the force.
    RealType strength = 0.0025;

    //! \brief The ratio r0/r1 factor, where r0 is the sum of the particle radii.
    RealType ratio = 1.;

    //! \brief The minimum allowable force (a large negative force). The force that occurs at the inner cutoff
    RealType inner_force = -1.;

    //! \brief Energy shift caused by cutoff.
    RealType potential_energy_shift = 0.;
  };

  // --- Define specific forces

  //! \brief Define the LennardJonesVLP class to be lennard jones force using verlet list pairs container.
  template<int dims> 
  using BuckinghamVLP = BuckinghamGeneric<dims, VerletListPairs>;

};

#endif // __HARD_SPHERE_GENERIC_HPP__GFLOW__