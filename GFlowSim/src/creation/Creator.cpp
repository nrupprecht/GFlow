#include "Creator.hpp"

namespace GFlow {
  
  SimData* Creator::create() {
    // Create bounds and hand them to a sim data object
    RealType width = 4, height = 4;
    Bounds simBounds(0,width,0,height);
    SimData* simData = new SimData(simBounds, simBounds);

    // Set wrapping
    simData->setWrap(false);

    RealType edge = 0.;

    // Add boundary walls
    simData->addWall(Wall(edge,edge,edge,height-edge)); // Left   wall
    simData->addWall(Wall(edge,edge,width-edge,edge)); // Bottom wall
    simData->addWall(Wall(width-edge,edge,width-edge,height-edge)); // Right  wall
    simData->addWall(Wall(edge,height-edge,height-edge,height-edge)); // Top    wall

    // Make room for particles
    int domain_size = 1018, edge_size = 0;
    simData->reserve(domain_size, edge_size);
    // Add some particles
    RealType sigma = 0.05;
    RealType baseVelocity = 0.02;
    for (int i=0; i<domain_size; ++i) {
      // Random particle
      RealType X     = (width-2*sigma-2*edge)*drand48()+sigma+edge;
      RealType Y     = (height-2*sigma-2*edge)*drand48()+sigma+edge;
      RealType theta = 2*PI*drand48();
      RealType V     = baseVelocity*normal_dist(generator);
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
    verlet.addExternalForce(new ViscousDrag(1.)); // This is large enough
    verlet.initialize(0.25);
    verlet.integrate();
    // Remove drag force
    simData->clearExternalForces();

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

  bool Creator::createRegion(Region& region, SimData* simData) {
    if (region.left<simData->getSimBounds().left || simData->getSimBounds().right<region.right || region.bottom<simData->getSimBounds().bottom || region.top<simData->getSimBounds().top) return false;

    // List of particles we are creating
    vector<Particle> particles;

    // First create radii list (this also allocates the correct number of particles)
    region.sigma->makeValues(region, particles);

    // Assign interaction
    region.interaction->makeValues(region, particles);

    // Assign positions and velocities
    region.position->makeValues(region, particles);
    region.theta->makeValues(region, particles);
    region.velocity->makeValues(region, particles);
    region.omega->makeValues(region, particles);

    // Assign inertias (masses and moments of inertia)
    region.inertia->makeValues(region, particles);

    // Assign repulsion, dissipation, coefficient of friction
    region.repulsion->makeValues(region, particles);
    region.dissipation->makeValues(region, particles);
    region.coeff->makeValues(region, particles);


    return true;
  }
  
}
