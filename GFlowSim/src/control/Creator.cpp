#include "Creator.hpp"

namespace GFlow {
  
  SimData* Creator::create() {
    // Create bounds and hand them to a sim data object
    Bounds simBounds(0,4,0,4);
    SimData* simData = new SimData(simBounds, simBounds);

    // Set wrapping
    simData->setWrap(false);

    // Add boundary walls
    simData->addWall(Wall(0,0,0,4)); // Left   wall
    simData->addWall(Wall(0,0,4,0)); // Bottom wall
    simData->addWall(Wall(4,0,4,4)); // Right  wall
    simData->addWall(Wall(0,4,4,4)); // Top    wall

    // Make room for particles
    int domain_size = 100, edge_size = domain_size/10;
    simData->reserve(domain_size, edge_size);
    // Add some particles
    RealType sigma = 0.05;
    for (int i=0; i<domain_size; ++i) {
      // Random particle
      RealType X     = (4-2*sigma)*drand48()+sigma;
      RealType Y     = (4-2*sigma)*drand48()+sigma;
      RealType theta = 2*PI*drand48();
      RealType V     = normal_dist(generator);
      // Set particle values
      Particle P(X, Y, sigma);
      P.velocity = V*vec2(cos(theta), sin(theta));
      P.theta = theta;
      P.dissipation = 0;
      P.coeff = 0;
      // Add the particle
      simData->addParticle(P);
    }

    return simData;
  }
  
}
