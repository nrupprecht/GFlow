/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 16, 2017
 *
 */

#ifndef __DATA_CREATOR_HPP__
#define __DATA_CREATOR_HPP__

// Includes
#include "../../include/Utility.hpp"
#include "../../include/DefaultConstants.hpp"
#include "../objects/Particle.hpp"
#include "../objects/Wall.hpp"
#include "../creation/Region.hpp"

namespace GFlow {

  // Forward declaration to region
  struct Region;

  // Data creation class
  struct DataCreator {
    virtual void makeValues(Region&, vector<Particle>&) = 0;
  };

  // Square lattice
  struct Square_Lattice : public DataCreator {
    Square_Lattice(RealType sig) : sigma(sig), packing(1.) {};
    Square_Lattice(RealType sig, RealType pack) : sigma(sig), packing(pack) {};
    virtual void makeValues(Region&, vector<Particle>&);
  private:
    RealType sigma, packing;
  };

  // Radius (and number)
  struct Fixed_Number_Uniform_Radii : public DataCreator {
    Fixed_Number_Uniform_Radii(int num, RealType sig, RealType disp) : number(num), sigma(sig), dispersion(disp) {};
    virtual void makeValues(Region&, vector<Particle>&);
  private:
    int number;
    RealType sigma, dispersion;
  };

  // Radius (and number)
  struct Fixed_Phi_Uniform_Radii : public DataCreator {
  public:
    Fixed_Phi_Uniform_Radii(RealType ph, RealType sig, RealType disp) : phi(ph), sigma(sig), dispersion(disp) {};
    virtual void makeValues(Region&, vector<Particle>&);
  private:
    RealType phi, sigma, dispersion;
  };

  // Radius (and number)
  struct Fixed_Phi_Scale_Free_Radii : public DataCreator {
  public:
    Fixed_Phi_Scale_Free_Radii(RealType ph, RealType smin, RealType smax) : phi(ph), sigMin(smin), sigMax(smax) {};
    virtual void makeValues(Region&, vector<Particle>&);
  private:
    RealType phi, sigMin, sigMax;
  };

  // Position
  struct Uniform_Space_Distribution : public DataCreator {
    virtual void makeValues(Region&, vector<Particle>&);
  };

  // Velocity
  struct Normal_Random_Velocity : public DataCreator {
    Normal_Random_Velocity() : velocity(0.25) {};
    Normal_Random_Velocity(RealType vel) : velocity(vel) {};
    virtual void makeValues(Region&, vector<Particle>&);
  private:
    RealType velocity;
  };

  // Velocity 
  struct Normal_Random_KE : public DataCreator {
    Normal_Random_KE() : vsgma(0.25) {};
    Normal_Random_KE(RealType s) : vsgma(s) {};
    virtual void makeValues(Region&, vector<Particle>&);
  private:
    RealType vsgma;
  };

  // Velocity
  struct Constant_Velocity : public DataCreator {
    Constant_Velocity() : velocity(Zero) {};
    Constant_Velocity(RealType vx, RealType vy) : velocity(vec2(vx,vy)) {};
    virtual void makeValues(Region&, vector<Particle>&);
  private:
    vec2 velocity;
  };

  // Theta
  struct Uniform_Random_Theta : public DataCreator {
    virtual void makeValues(Region&, vector<Particle>&);
  };

  // Theta
  struct Align_Velocity_Theta : public DataCreator {
    virtual void makeValues(Region&, vector<Particle>&);
  };

  // Omega
  struct Normal_Random_Omega : public DataCreator {
    Normal_Random_Omega() : omega(0) {};
    Normal_Random_Omega(RealType om) : omega(om) {};
    virtual void makeValues(Region&, vector<Particle>&);
  private:
    RealType omega;
  };

  // Dissipation
  struct Uniform_Random_Dissipation : public DataCreator {
    Uniform_Random_Dissipation() : value(default_particle_dissipation), dispersion(0) {};
    Uniform_Random_Dissipation(RealType v, RealType dis=0) : value(v), dispersion(dis) {};
    virtual void makeValues(Region&, vector<Particle>&);
  private:
    RealType value, dispersion;
  };

  // Repulsion
  struct Uniform_Random_Repulsion : public DataCreator {
    Uniform_Random_Repulsion() : value(default_particle_repulsion), dispersion(0) {};
    Uniform_Random_Repulsion(RealType v, RealType dis=0) : value(v), dispersion(dis) {};
    virtual void makeValues(Region&, vector<Particle>&);
  private:
    RealType value, dispersion;
  };

  // Coeff
  struct Uniform_Random_Coeff : public DataCreator {
    Uniform_Random_Coeff() : value(default_particle_coeff), dispersion(0) {};
    Uniform_Random_Coeff(RealType v, RealType dis=0) : value(v), dispersion(0) {};
    virtual void makeValues(Region&, vector<Particle>&);
  private:
    RealType value, dispersion;
  };

  // Interaction
  struct Homogeneous_Interaction : public DataCreator {
    Homogeneous_Interaction() : interaction(0) {};
    Homogeneous_Interaction(int inter) : interaction(inter) {};
    virtual void makeValues(Region&, vector<Particle>&);
  private:
    int interaction;
  };

  // Inertia (mass and moment of inertia)
  struct Constant_Density : public DataCreator {
    Constant_Density() : density(1.), dispersion(0) {}
    Constant_Density(RealType den, RealType dis=0) : density(den), dispersion(dis) {};
    virtual void makeValues(Region&, vector<Particle>&);
  private:
    RealType density, dispersion;
  };

}
#endif // __DATA_CREATOR_HPP__
