#include "DataCreator.hpp" // Our Header file
#include "../creation/FileParser.hpp"
#include "../creation/Creator.hpp"

namespace GFlow {

  void Fixed_Number_Uniform_Radii::makeValues(Region& region, vector<Particle>& particles) {
    for (int i=0; i<number; ++i) {
      RealType s = sigma*(1-dispersion*drand48());
      particles.push_back(Particle(0,0,s));
    }
  }

  void Fixed_Phi_Uniform_Radii::makeValues(Region& region, vector<Particle>& particles) {
    // Get parameters
    RealType volume = region.bounds.volume();
    RealType aim = phi*volume, vol = 0;

    // Add more radii until we have the desired packing ratio
    while (vol<aim) {
      RealType s = sigma*(1-dispersion*drand48());
      particles.push_back(Particle(0,0,s));
      vol += PI*sqr(s);
    }
  }

  void Fixed_Phi_Uniform_Radii::makeValues(Region& region, vector<Particle>& particles)  {
    // Get parameters
    RealType volume = region.bounds.volume();
    RealType aim = phi*volume, vol = 0;

    // Add more radii until we have the desired packing ratio
    while (vol<aim) {
      RealType s = sigma*(1-dispersion*drand48());
      particles.push_back(Particle(0,0,s));
      vol += PI*sqr(s);
    }
  }

  void Uniform_Space_Distribution::makeValues(Region& region, vector<Particle>& particles) {
    // Get parameters
    RealType left = region.bounds.left, width = region.bounds.right - left;
    RealType bottom = region.bounds.bottom, height = region.bounds.top - bottom;

    // Distribute positions uniformly
    for (auto& p : particles) p.position = vec2(left+p.sigma+(width-2*p.sigma)*drand48(), bottom+p.sigma+(height-2*p.sigma)*drand48());
  }

  void Normal_Random_Velocity::makeValues(Region& region, vector<Particle>& particles) {
    if (velocity!=0)
      for (auto& p : particles) {
	RealType theta = 2*PI*drand48();
	RealType vel = velocity*normal_dist(generator);
	p.velocity = vec2(vel*cos(theta), vel*sin(theta));
      }
    else 
      for (auto& p : particles) p.velocity = Zero;
  }

  void Normal_Random_KE::makeValues(Region& region, vector<Particle>& particles) {
    if (vsgma!=0)
      for (auto& p : particles) {
        RealType theta = 2*PI*drand48();
        vec2 v(cos(theta), sin(theta));
	double ke = fabs(vsgma*normal_dist(generator));
	double velocity = sqrt(2*p.invMass*ke/127.324);
        p.velocity = velocity*v;
      }
    else
      for (auto& p : particles) p.velocity = Zero;
  }

  void Constant_Velocity::makeValues(Region& region, vector<Particle>& particles) {
    for (auto& p : particles) p.velocity = velocity;
  }

  void Uniform_Random_Theta::makeValues(Region& region, vector<Particle>& particles) {
    for (auto& p : particles) p.theta = 2*PI*drand48();
  }

  void Align_Velocity_Theta::makeValues(Region& region, vector<Particle>& particles) {
    for (auto& p : particles) p.theta = atan2(p.velocity.y, p.velocity.x);
  }

  void Normal_Random_Omega::makeValues(Region& region, vector<Particle>& particles) {
    if (omega!=0)
      for (auto& p : particles) p.omega = omega*normal_dist(generator);
    else
      for (auto& p : particles) p.omega = 0;
  }

  void Uniform_Random_Dissipation::makeValues(Region& region, vector<Particle>& particles) {
    for (auto& p : particles) p.dissipation = value*(1-dispersion*drand48());
  }

  void Uniform_Random_Repulsion::makeValues(Region& region, vector<Particle>& particles) {
    for (auto& p : particles) p.repulsion = value*(1-dispersion*drand48());
  }

  void Uniform_Random_Coeff::makeValues(Region& region, vector<Particle>& particles) {
    for (auto& p : particles) p.coeff = value*(1-dispersion*drand48());
  }

  void Homogeneous_Interaction::makeValues(Region& region, vector<Particle>& particles) {
    for (auto& p : particles) p.interaction = interaction;
  }

  void Constant_Density::makeValues(Region& region, vector<Particle>& particles) {
    for (auto& p : particles) {
      RealType mass = PI*sqr(p.sigma)*density*(1-dispersion*drand48());
      p.invMass = 1./mass;
      p.invII = 1./(0.5*mass*sqr(p.sigma));
    }
  }
  
}
