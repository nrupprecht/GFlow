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
    RealType volume = (region.right - region.left)*(region.top - region.bottom);
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
    RealType left = region.left, width = region.right - left;
    RealType bottom = region.bottom, height = region.top - bottom;

    // Distribute positions uniformly
    for (auto& p : particles) p.position = vec2(left+width*drand48(), bottom+height*drand48());
  }

  void Uniformly_Random_Dissipation::makeValues(Region& region, vector<Particle>& particles) {
    for (auto& p : particles) p.dissipation = value*(1-dispersion*drand48());
  }

  void Uniformly_Random_Repulsion::makeValues(Region& region, vector<Particle>& particles) {
    for (auto& p : particles) p.repulsion = value*(1-dispersion*drand48());
  }

  void Uniformly_Random_Coeff::makeValues(Region& region, vector<Particle>& particles) {
    for (auto& p : particles) p.coeff = value*(1-dispersion*drand48());
  }
  
}
