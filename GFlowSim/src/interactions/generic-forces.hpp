#ifndef __HARD_SPHERE_GENERIC_HPP__GFLOW__
#define __HARD_SPHERE_GENERIC_HPP__GFLOW__

#include "generic-handlers.hpp"
// Other files
#include "../allintegrators.hpp"

namespace GFlowSimulation {

  /*
  * \brief Generic class for basic hard sphere forces.
  *
  */
  template<int dims, template<int, class> class handler> class HardSphereGeneric : public handler<dims, HardSphereGeneric<dims, handler> > {
  public:
    //! \brief Typedef for the handler type.
    typedef handler<dims, HardSphereGeneric<dims, handler> > Handler;

    //! \brief Constructor.
    HardSphereGeneric(GFlow *gflow) : Handler(gflow) {};

    //! \brief Constructor, sets repulsion
    HardSphereGeneric(GFlow *gflow, RealType rp) : Handler(gflow), repulsion(rp) {};

    //! \brief The force kernel for hard sphere forces.
    virtual void force(const int id1, const int id2, const RealType R1, const RealType R2, const RealType rsqr, RealType *X, RealType **f) const {
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

  template<int dims> 
  using HardSphereVL = HardSphereGeneric<dims, VerletList>;



  /*
  * \brief Generic class for hard sphere forces with dissipation.
  *
  */
  template<int dims, template<int, class> class handler> class HardSphereDsGeneric : public handler<dims, HardSphereDsGeneric<dims, handler> > {
  public:
    //! \brief Typedef for the handler type.
    typedef handler<dims, HardSphereDsGeneric<dims, handler> > Handler;

    //! \brief Constructor.
    HardSphereDsGeneric(GFlow *gflow) : Handler(gflow) {};

    //! \brief Constructor, sets repulsion
    HardSphereDsGeneric(GFlow *gflow, RealType rp) : Handler(gflow), repulsion(rp) {};

    //! \brief The force kernel for hard sphere forces.
    virtual void force(const int id1, const int id2, const RealType R1, const RealType R2, const RealType rsqr, RealType *X, RealType **f) const {
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
    //! \brief Typedef for the handler type.
    typedef handler<dims, LennardJonesGeneric<dims, handler> > Handler;

    //! \brief Constructor.
    LennardJonesGeneric(GFlow *gflow) : Handler(gflow) {
      setCutoff(2.5);
    };

    //! \brief Constructor, sets strength.
    LennardJonesGeneric(GFlow *gflow, RealType str) : Handler(gflow), strength(str) {
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

    virtual void force(const int id1, const int id2, const RealType R1, const RealType R2, const RealType rsqr, RealType *X, RealType **f) const {
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
    //! \brief Typedef for the handler type.
    typedef handler<dims, CoulombGeneric<dims, handler> > Handler;

    //! \brief Constructor.
    CoulombGeneric(GFlow *gflow) : Handler(gflow) {
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

    virtual void force(const int id1, const int id2, const RealType R1, const RealType R2, const RealType rsqr, RealType *X, RealType **f) const {
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

  //! \brief Define the CoulombVLP class to be Coulomb force using verlet list pairs container.
  template<int dims> 
  using CoulombVLP = CoulombGeneric<dims, VerletListPairs>;



  /*
  * \brief Generic class for Buckingham forces.
  *
  */
  template<int dims, template<int, class> class handler> class BuckinghamGeneric : public handler<dims, BuckinghamGeneric<dims, handler> > {
  public:
    //! \brief Typedef for the handler type.
    typedef handler<dims, BuckinghamGeneric<dims, handler> > Handler;

    //! \brief Constructor.
    BuckinghamGeneric(GFlow *gflow) : Handler(gflow) {
      setCutoff(3.);
    };

    //! \brief Constructor, sets strength.
    BuckinghamGeneric(GFlow *gflow, RealType str) : Handler(gflow), strength(str) {
      setCutoff(3.);
    };

    //! \brief Set the cutoff factor.
    void setCutoff(RealType c) { 
      Interaction::cutoff = c>0 ? c : Interaction::cutoff;
      // Calculate the potential when r is the cutoff radius
      potential_energy_shift = 0; //\todo Calculate this.
    }

    virtual void force(const int id1, const int id2, const RealType R1, const RealType R2, const RealType rsqr, RealType *X, RealType **f) const {
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
      Fn = min(inner_force, Fn);

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



  /*
  * \brief Generic class for Hertz-type forces.
  *
  */
  template<int dims, template<int, class> class handler> class HertzGeneric : public handler<dims, HertzGeneric<dims, handler> > {
  public:
    //! \brief Typedef for the handler type.
    typedef handler<dims, HertzGeneric<dims, handler> > Handler;

    //! \brief Constructor.
    HertzGeneric(GFlow *gflow) : Handler(gflow) {
      // Add the needed data entries and integrator.
      if (dims==2) {
        om_add = Base::simData->requestScalarData("Om");
        tq_add = Base::simData->requestScalarData("Tq");
        // Add an angular integrator.
        Base::gflow->addIntegrator(new AngularVelocityVerlet2d(gflow));
      }
    };

    virtual void interact() const override {
      // Make sure pointers are correct.
      if (-1<om_add && -1<tq_add) {
        om = Base::simData->ScalarData(om_add);
        tq = Base::simData->ScalarData(tq_add);
      }

      // Now do interactions.
      Handler::interact();
    }

    void force(const int id1, const int id2, const RealType R1, const RealType R2, const RealType rsqr, RealType *X, RealType **f) const {
      // Calculate distance, inverse distance.
      RealType r = sqrt(rsqr);
      RealType invr = 1./r;
      // Create a normal vector
      scalar_mult_eq_vec<dims>(X, invr);

      // Necessary quantities.
      RealType delta = R1 + R2 - r;
      RealType R_eff = R1*R2/(R1+R2);
      RealType m1 = 1./Base::simData->Im(id1);
      RealType m2 = 1./Base::simData->Im(id2);
      RealType M_eff = m1*m2/(m1+m2);

      // In case one or both particles are infinitely massive.
      M_eff = isnan(M_eff) ? 1. : M_eff;
      
      // Calculate relative velocity.
      RealType V[dims];
      subtract_vec<dims>(Base::simData->V(id1), Base::simData->V(id2), V);
      // Calculate normal velocity.
      RealType vn = dot_vec<dims>(V, X);
      // Calculate prefactor.
      RealType c1 = sqrt(delta*R_eff);

      // Calculate the magnitude of the force
      RealType Fn = max(RealType(0.), c1 * (K_n*delta - M_eff*gamma_n*vn));

      // Update normal forces.
      sum_eq_vec_scaled<dims>(f[id1], X,  Fn);
      sum_eq_vec_scaled<dims>(f[id2], X, -Fn);

      // Calculate tangential force.
      RealType Vn[dims], Vt[dims];
      // Create the Vn (normal velocity component) vector.
      copy_vec<dims>(X, Vn);
      scalar_mult_eq_vec<dims>(Vn, vn);
      // Subtract, setting Vt = V - Vn.
      subtract_vec<dims>(V, Vn, Vt);

      // Angular velocities.
      if (dims==2 && om!=nullptr && tq!=nullptr) {
        // Tangential normal direction for 2d.
        RealType Nt[dims];
        Nt[0] = -X[1];
        Nt[1] =  X[0];

        RealType vt = dot_vec<dims>(Vt, Nt); //sqrt(dot_vec<dims>(Vt, Vt));        

        // Calculate difference in velocities at the point of intersection.
        vt -= om[id1]*R1 + om[id2]*R2; // Relative surface velocity due to angular velocity.

        // Instead of K_t * 0, it should really be based on an "angular spring" stretching, 
        // see e.g. Luding, Gran. Matter 2008, v 10, p. 235. But I'd have to keep track of the
        // history of particles, even between neighbor list updates, which is too much of a pain 
        // for now.
        RealType Ft = - c1 * (K_t*0 + M_eff*vt);

        RealType maxF = fabs(mu*Fn);
        if (Ft>maxF) Ft = mu*Fn;
        else if (Ft<-maxF) Ft = -maxF;

        // Update tangential forces.
        sum_eq_vec_scaled<dims>(f[id1], Nt, Ft);
        sum_eq_vec_scaled<dims>(f[id2], Nt, -Ft);

        // Update torques.
        tq[id1] += Ft*R1;
        tq[id2] += Ft*R2;
      }

      // Calculate potential
      if (Interaction::do_potential) {
        //Interaction::potential += 0;
      }
      // Calculate virial
      if (Interaction::do_virial) {
        Interaction::virial += Fn*r;
      }      
    }

    void setKn(const RealType kn) { K_n = 0<=kn ? kn : K_n; }
    void setKt(const RealType kt) { K_t = 0<=kt ? kt : K_t; }

    void setGammaN(const RealType gn) { gamma_n = 0<=gn ? gn : gamma_n; }
    void setGammaT(const RealType gt) { gamma_t = 0<=gt ? gt : gamma_t; }

    void setMu(const RealType m) { mu = 0<=m ? m : mu; }

  private:
    //! \brief Normal elastic constant.
    RealType K_n = 5'000.;
    //! \brief Tangential elastic constant.
    RealType K_t = 5'000. * 2./7.;

    //! \brief Normal viscoelastic constant.
    RealType gamma_n = 1000.;
    //! \brief Tangential viscoelastic constant.
    RealType gamma_t = 700.;

    //! \brief The largest possible (absolute value) friction force is mu*Fn.
    RealType mu = 0.5;

    //! \brief The array address for angular velocity.
    int om_add = -1;
    //! \brief The array address for torque.
    int tq_add = -1;

    //! \brief Pointers to necessary data.
    mutable RealType *om = nullptr, *tq = nullptr;
  };

  // --- Define specific forces

  //! \brief Define the HertzVLP class to be Hertz force using verlet list pairs container.
  template<int dims> 
  using HertzVLP = HertzGeneric<dims, VerletListPairs>;



  /*
  * \brief Generic class for Hertz-type forces.
  *
  */
  template<int dims, template<int, class> class handler> class HookeGeneric : public handler<dims, HookeGeneric<dims, handler> > {
  public:
    //! \brief Typedef for the handler type.
    typedef handler<dims, HookeGeneric<dims, handler> > Handler;

    //! \brief Constructor.
    HookeGeneric(GFlow *gflow) : Handler(gflow) {
      // Add the needed data entries and integrator.
      if (dims==2) {
        om_add = Base::simData->requestScalarData("Om");
        tq_add = Base::simData->requestScalarData("Tq");
        // Add an angular integrator.
        Base::gflow->addIntegrator(new AngularVelocityVerlet2d(gflow));
      }
    };

    virtual void interact() const override {
      // Make sure pointers are correct.
      if (-1<om_add && -1<tq_add) {
        om = Base::simData->ScalarData(om_add);
        tq = Base::simData->ScalarData(tq_add);
      }

      // Now do interactions.
      Handler::interact();
    }

    void force(const int id1, const int id2, const RealType R1, const RealType R2, const RealType rsqr, RealType *X, RealType **f) const {
      // Calculate distance, inverse distance.
      RealType r = sqrt(rsqr);
      RealType invr = 1./r;
      // Create a normal vector
      scalar_mult_eq_vec<dims>(X, invr);

      // Necessary quantities.
      RealType delta = R1 + R2 - r;
      RealType m1 = 1./Base::simData->Im(id1);
      RealType m2 = 1./Base::simData->Im(id2);
      RealType M_eff = m1*m2/(m1+m2);

      // In case one or both particles are infinitely massive.
      M_eff = isnan(M_eff) ? 1. : M_eff;
      
      // Calculate relative velocity.
      RealType V[dims];
      subtract_vec<dims>(Base::simData->V(id1), Base::simData->V(id2), V);
      // Calculate normal velocity.
      RealType vn = dot_vec<dims>(V, X);

      // Calculate the magnitude of the force
      RealType Fn = max(RealType(0.), K_n*delta - M_eff*gamma_n*vn);

      // Update normal forces.
      sum_eq_vec_scaled<dims>(f[id1], X,  Fn);
      sum_eq_vec_scaled<dims>(f[id2], X, -Fn);

      // Calculate tangential force.
      RealType Vn[dims], Vt[dims];
      // Create the Vn (normal velocity component) vector.
      copy_vec<dims>(X, Vn);
      scalar_mult_eq_vec<dims>(Vn, vn);
      // Subtract, setting Vt = V - Vn.
      subtract_vec<dims>(V, Vn, Vt);

      // Angular velocities.
      if (dims==2 && om!=nullptr && tq!=nullptr) {
        // Tangential normal direction for 2d.
        RealType Nt[dims];
        Nt[0] = -X[1];
        Nt[1] =  X[0];

        RealType vt = dot_vec<dims>(Vt, Nt); //sqrt(dot_vec<dims>(Vt, Vt));        

        // Calculate difference in velocities at the point of intersection.
        vt -= om[id1]*R1 + om[id2]*R2; // Relative surface velocity due to angular velocity.

        // Instead of K_t * 0, it should really be based on an "angular spring" stretching, 
        // see e.g. Luding, Gran. Matter 2008, v 10, p. 235. But I'd have to keep track of the
        // history of particles, even between neighbor list updates, which is too much of a pain 
        // for now.
        RealType Ft = - (K_t*0 + M_eff*vt);

        RealType maxF = fabs(mu*Fn);
        if (Ft>maxF) Ft = mu*Fn;
        else if (Ft<-maxF) Ft = -maxF;

        // Update tangential forces.
        sum_eq_vec_scaled<dims>(f[id1], Nt, Ft);
        sum_eq_vec_scaled<dims>(f[id2], Nt, -Ft);

        // Update torques.
        tq[id1] += Ft*R1;
        tq[id2] += Ft*R2;
      }

      // Calculate potential
      if (Interaction::do_potential) {
        //Interaction::potential += 0;
      }
      // Calculate virial
      if (Interaction::do_virial) {
        Interaction::virial += Fn*r;
      }      
    }

    void setKn(const RealType kn) { K_n = 0<=kn ? kn : K_n; }
    void setKt(const RealType kt) { K_t = 0<=kt ? kt : K_t; }

    void setGammaN(const RealType gn) { gamma_n = 0<=gn ? gn : gamma_n; }
    void setGammaT(const RealType gt) { gamma_t = 0<=gt ? gt : gamma_t; }

    void setMu(const RealType m) { mu = 0<=m ? m : mu; }

  private:
    //! \brief Normal elastic constant.
    RealType K_n = 5'000.;
    //! \brief Tangential elastic constant.
    RealType K_t = 5'000. * 2./7.;

    //! \brief Normal viscoelastic constant.
    RealType gamma_n = 1000.;
    //! \brief Tangential viscoelastic constant.
    RealType gamma_t = 700.;

    //! \brief The largest possible (absolute value) friction force is mu*Fn.
    RealType mu = 0.5;

    //! \brief The array address for angular velocity.
    int om_add = -1;
    //! \brief The array address for torque.
    int tq_add = -1;

    //! \brief Pointers to necessary data.
    mutable RealType *om = nullptr, *tq = nullptr;
  };

  // --- Define specific forces

  //! \brief Define the HookeVLP class to be Hooke force using verlet list pairs container.
  template<int dims> 
  using HookeVLP = HookeGeneric<dims, VerletListPairs>;

};

#endif // __HARD_SPHERE_GENERIC_HPP__GFLOW__