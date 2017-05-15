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
    int domain_size = 1018, edge_size = 0;
    simData->reserve(domain_size, edge_size);
    // Add some particles
    RealType sigma = 0.05;
    RealType edge = 0.5;
    for (int i=0; i<domain_size; ++i) {
      // Random particle
      RealType X     = (4-2*sigma-2*edge)*drand48()+sigma+edge;
      RealType Y     = (4-2*sigma-2*edge)*drand48()+sigma+edge;
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

    // Relax, so there is no overlap
    VelocityVerletIntegrator verlet(simData);
    DragForce *largeDrag = new ViscousDrag(1.); // This is large enough
    verlet.addDragForce(largeDrag);
    verlet.initialize(0.25);
    verlet.integrate();
    // Clean up drag force
    delete largeDrag;

    // Give random velocities
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    RealType *ds = simData->getDsPtr();
    RealType *cf = simData->getCfPtr();
    for (int i=0; i<domain_size; ++i) {
      // Random particle
      RealType V     = normal_dist(generator);
      RealType theta = 2*PI*drand48();
      // Set particle values
      vx[i] = V*cos(theta);
      vy[i] = V*sin(theta);
      // Set particle values
      ds[i] = 0;
      cf[i] = 0;
    }

    return simData;
  }
  
}
