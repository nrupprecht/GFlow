#include "Creator.hpp"

namespace GFlow {
  
  void Creator::create(SimData *& simData, Integrator *& integrator) {
    // Create bounds and hand them to a sim data object
    RealType width = 4, height = 4;
    Bounds simBounds(0,width,0,height);
    simData = new SimData(simBounds, simBounds);

    // Set wrapping
    simData->setWrap(false);

    RealType edge = 0.;

    // Add boundary walls
    Wall w(edge,edge,edge,height-edge);  // Left wall
    w.dissipation = 0;
    w.coeff = 0;
    simData->addWall(w);

    w = Wall(edge,edge,width-edge,edge); // Bottom wall
    w.dissipation = 0;
    w.coeff = 0;
    simData->addWall(w);

    w = Wall(width-edge,edge,width-edge,height-edge); // Right wall
    w.dissipation = 0;
    w.coeff = 0;
    simData->addWall(w);

    w = Wall(edge,height-edge,height-edge,height-edge); // Top wall
    w.dissipation = 0;
    w.coeff = 0;
    simData->addWall(w);

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

    // Create integrator
    integrator = new VelocityVerletIntegrator(simData);
  }

  void Creator::createBuoyancy(SimData *& simData, Integrator *& integrator, RealType radius, RealType density, vec2 velocity) {
    Bounds sb = simData->getSimBounds();
    // Change the bounds for the addition of the new particle
    RealType height = 1; //** Height of the tops of the balls
    Bounds nb(sb.left, sb.right, sb.bottom, sb.top+2*radius);
    // Set new bounds
    simData->setSimBounds(nb);
    simData->setBounds(nb);
    // Get rid of the old walls
    simData->getWalls().clear();
    // Create new walls
    simData->addWall(Wall(nb.left, nb.bottom, nb.right, nb.bottom)); // Bottom
    simData->addWall(Wall(nb.left, nb.bottom, nb.left, nb.top));     // Left
    simData->addWall(Wall(nb.right, nb.bottom, nb.right, nb.top));   // Right
    simData->addWall(Wall(nb.left, nb.top, nb.right, nb.top));       // Top

    // Insert the intruding particle
    if (radius>0) {
      Particle P(0.5*(nb.right-nb.left), height + radius, radius);
      P.setDensity(density);
      P.velocity = velocity;
    }

    // Make sure we use a small time step
    reinterpret_cast<VelocityVerletIntegrator*>(integrator)->setMinTimeStep(1e-4);
  }

  bool Creator::createRegion(Region& region, SimData* simData) {
    // Make sure the region is within the bounds
    if (!simData->simBounds.contains(region.bounds)) return false;

    // List of particles we are creating
    vector<Particle> particles;

    // First create radii list (this also allocates the correct number of particles)
    if (region.sigma) region.sigma->makeValues(region, particles);
    else return false; // We require a sigma function

    // Assign interaction
    if (region.interaction) region.interaction->makeValues(region, particles);
    else Homogeneous_Interaction().makeValues(region, particles);

    // Assign positions
    if (region.position) region.position->makeValues(region, particles);
    else Uniform_Space_Distribution().makeValues(region, particles);
    
    // Assign inertias (masses and moments of inertia)
    if (region.inertia) region.inertia->makeValues(region, particles);
    else Constant_Density().makeValues(region, particles);

    // Assign repulsion
    if (region.repulsion) region.repulsion->makeValues(region, particles);
    else Uniform_Random_Repulsion().makeValues(region, particles);

    // Assign dissipation
    if (region.dissipation) region.dissipation->makeValues(region, particles);
    else Uniform_Random_Dissipation().makeValues(region, particles);

    // Assign coefficient of friction
    if (region.coeff) region.coeff->makeValues(region, particles);
    else Uniform_Random_Coeff().makeValues(region, particles);

    // Let the simulation relax, using heavy viscous drag, so particles do not overlap with things
    SimData *relax = new SimData(region.bounds, region.bounds);
    relax->setWrap(false);
    relax->reserve(particles.size());
    relax->addParticle(particles);
    relax->addWall(region.bounds);
    VelocityVerletIntegrator verlet(relax);
    verlet.setAdjustTimeStep(false);
    verlet.addExternalForce(new ViscousDrag(1.)); // This is large enough
    verlet.initialize(0.25);
    verlet.integrate();
    // Remove drag force
    particles = relax->getParticles();
    delete relax;
    
    // Set velocities and omegas
    // Assign velocity
    if (region.velocity) region.velocity->makeValues(region, particles);
    else Normal_Random_Velocity().makeValues(region, particles);

    // Assign theta
    if (region.theta) region.theta->makeValues(region, particles);
    else Uniform_Random_Theta().makeValues(region, particles);

    // Assign omega
    if (region.omega) region.omega->makeValues(region, particles);
    else Normal_Random_Omega().makeValues(region, particles);    
    
    // Reserve more room for the new particles
    simData->reserveAdditional(particles.size(), 0);
    // Add the particles to simData
    simData->addParticle(particles);

    // Return success
    return true;
  }
  
}
