#include "Creator.hpp"

namespace GFlow {
  
  SimData* Creator::create() {
    
    // Create bounds and hand them to a sim data object
    Bounds simBounds(0,4,0,4);
    SimData* simData = new SimData(simBounds, simBounds);

    // Add boundary walls
    simData->addWall(Wall(0,0,0,4));
    simData->addWall(Wall(0,0,4,0));
    simData->addWall(Wall(4,0,4,4));
    simData->addWall(Wall(0,4,4,4));

    // Add some particles
    int domain_size = 1000, edge_size = 100;
    simData->reserve(domain_size, edge_size);
    
    RealType sigma = 0.05;
    for (int i=0; i<domain_size; ++i) {
      RealType X     = 4*drand48();
      RealType Y     = 4*drand48();
      RealType theta = 2*PI*drand48();
      RealType V     = normal_dist(generator);
      Particle P(X, Y, sigma);

      P.velocity = V*vec2(cos(theta), sin(theta));
      simData->addParticle(P);
    }

    return simData;
  }
  
}
