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
      vec2 pos(4*drand48(), 4*drand48());
      simData->addParticle(Particle(pos, sigma));
    }

    return simData;
  }
  
}
